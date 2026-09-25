#include "sphsolver.h"
#include "RigidBody.h"
#include <cmath>
#include <algorithm>

SPHSolver::SPHSolver(Parameters params) : params_(params) {}

double SPHSolver::computeSoundSpeed(double u, double gamma) const {
    return std::sqrt(gamma * (gamma - 1.0) * std::max(u, 1e-6));
}

void SPHSolver::updateSmoothingLengths(std::vector<std::shared_ptr<Particle>>& fluidParticles) const {
    for (auto& p : fluidParticles) {
        double rho = std::max(p->density, 1e-5);
        p->h = params_.eta * std::sqrt(p->mass / rho);
    }
}

void SPHSolver::computeDensityAndPressure(
    std::vector<std::shared_ptr<Particle>>& fluidParticles,
    const std::vector<std::shared_ptr<Particle>>& internalRigidParticles,
    const std::vector<std::shared_ptr<Particle>>& ghostParticles,
    const SpatialGrid& grid,
    const Kernel& kernel) const 
{
    for (auto& p_i : fluidParticles) {
        p_i->density = 0.0;
        double h_i = p_i->h;

        // 1. Fluid-Fluid Density Contributions
        auto neighborIndices = grid.getNeighborIndices(p_i->pos);
        for (size_t idx : neighborIndices) {
            if (idx >= fluidParticles.size()) continue;
            const auto& p_j = fluidParticles[idx];
            Vector2D r_ij = p_i->pos - p_j->pos;
            p_i->density += p_j->mass * kernel.value(r_ij, h_i);
        }

        // 2. Internal Rigid Body Particle Contributions
        for (const auto& rp : internalRigidParticles) {
            Vector2D r_ij = p_i->pos - rp->pos;
            p_i->density += rp->mass * kernel.value(r_ij, h_i);
        }

        // 3. External Ghost Boundary Contributions
        for (const auto& ghost : ghostParticles) {
            Vector2D r_ij = p_i->pos - ghost->pos;
            p_i->density += ghost->mass * kernel.value(r_ij, h_i);
        }

        if (p_i->density < 1e-5) p_i->density = 1e-5;
    }
}

void SPHSolver::computeDerivatives(
    std::vector<std::shared_ptr<Particle>>& fluidParticles,
    const std::vector<std::shared_ptr<Particle>>& internalRigidParticles,
    const std::vector<std::shared_ptr<Particle>>& ghostParticles,
    const SpatialGrid& grid,
    const Kernel& kernel) const
{
    size_t n = fluidParticles.size();

    for (auto& p : fluidParticles) {
        p->accel = params_.enableGravity ? params_.gravity : Vector2D(0.0, 0.0);
        p->dudt = 0.0;
    }

    // 1. Fluid-Fluid Interactions
    for (size_t i = 0; i < n; ++i) {
        auto& p_i = fluidParticles[i];
        double rho_i = p_i->density;
        double u_i   = p_i->u;
        double P_i   = (params_.gamma - 1.0) * rho_i * u_i;
        double h_i   = p_i->h;
        double c_i   = computeSoundSpeed(u_i, params_.gamma);

        auto neighborIndices = grid.getNeighborIndices(p_i->pos);
        for (size_t idx : neighborIndices) {
            if (idx <= i || idx >= n) continue;
            const auto& p_j = fluidParticles[idx];

            Vector2D r_ij = p_i->pos - p_j->pos;
            double r2 = r_ij.normSq();
            double h_j  = p_j->h;
            double h_ij = 0.5 * (h_i + h_j);
            Vector2D gradW_ij = 0.5 * (kernel.gradient(r_ij, h_i) + kernel.gradient(r_ij, h_j));

            if (gradW_ij.normSq() < 1e-24) continue;

            Vector2D v_ij = p_i->vel - p_j->vel;
            double rho_j  = p_j->density;
            double u_j    = p_j->u;
            double P_j    = (params_.gamma - 1.0) * rho_j * u_j;
            double c_j    = computeSoundSpeed(u_j, params_.gamma);

            double pi_ij = 0.0;
            double v_dot_r = v_ij.dot(r_ij);
            if (v_dot_r < 0.0) {
                double rho_ij = 0.5 * (rho_i + rho_j);
                double c_ij   = 0.5 * (c_i + c_j);
                double mu_ij  = (h_ij * v_dot_r) / (r2 + 0.01 * h_ij * h_ij);
                pi_ij = (-params_.alpha * c_ij * mu_ij + params_.beta * mu_ij * mu_ij) / rho_ij;
            }

            double p_term = (P_i / (rho_i * rho_i)) + (P_j / (rho_j * rho_j)) + pi_ij;
            Vector2D force_ij = gradW_ij * (-p_i->mass * p_j->mass * p_term);
            double thermal_ij = 0.5 * p_term * v_ij.dot(gradW_ij);

            p_i->accel += force_ij / p_i->mass;
            p_j->accel -= force_ij / p_j->mass;

            p_i->dudt += p_j->mass * thermal_ij;
            p_j->dudt += p_i->mass * thermal_ij;
        }

        // 2. Fluid-Internal Rigid Particles Interaction
        for (const auto& rp : internalRigidParticles) {
            Vector2D r_ij = p_i->pos - rp->pos;
            double r2 = r_ij.normSq();
            double h_j  = rp->h;
            double h_ij = 0.5 * (h_i + h_j);
            Vector2D gradW_ij = 0.5 * (kernel.gradient(r_ij, h_i) + kernel.gradient(r_ij, h_j));

            if (gradW_ij.normSq() < 1e-24) continue;

            Vector2D v_ij = p_i->vel - rp->vel;
            double rho_j  = p_i->density;
            double P_j    = P_i; // Mirrored pressure to prevent boundary suction

            double pi_ij = 0.0;
            double v_dot_r = v_ij.dot(r_ij);
            if (v_dot_r < 0.0) {
                double rho_ij = rho_i;
                double mu_ij  = (h_ij * v_dot_r) / (r2 + 0.01 * h_ij * h_ij);
                pi_ij = (-params_.alpha * c_i * mu_ij + params_.beta * mu_ij * mu_ij) / rho_ij;
            }

            double p_term = (P_i / (rho_i * rho_i)) + (P_j / (rho_j * rho_j)) + pi_ij;
            Vector2D force_fluid_on_rp = gradW_ij * (p_i->mass * rp->mass * p_term);

            // Apply reactive force onto fluid
            p_i->accel -= force_fluid_on_rp / p_i->mass;
            p_i->dudt += 0.5 * rp->mass * p_term * v_ij.dot(gradW_ij);

            // Route equal-opposite force to the owning RigidObject at this specific location
            if (rp->parentBody) {
                rp->parentBody->addForceAtPosition(force_fluid_on_rp, rp->pos);
            }
        }

        // 3. Fluid-External Wall Ghost Interactions
        for (const auto& ghost : ghostParticles) {
            Vector2D r_ij = p_i->pos - ghost->pos;
            double r2 = r_ij.normSq();
            double h_j  = ghost->h;
            double h_ij = 0.5 * (h_i + h_j);
            Vector2D gradW_ij = 0.5 * (kernel.gradient(r_ij, h_i) + kernel.gradient(r_ij, h_j));

            if (gradW_ij.normSq() < 1e-24) continue;

            Vector2D v_ij = p_i->vel - ghost->vel;
            double P_j = P_i;

            double pi_ij = 0.0;
            double v_dot_r = v_ij.dot(r_ij);
            if (v_dot_r < 0.0) {
                double rho_ij = rho_i;
                double mu_ij  = (h_ij * v_dot_r) / (r2 + 0.01 * h_ij * h_ij);
                pi_ij = (-params_.alpha * c_i * mu_ij + params_.beta * mu_ij * mu_ij) / rho_ij;
            }

            double p_term = (P_i / (rho_i * rho_i)) + (P_j / (P_j * P_j)) + pi_ij;
            Vector2D f_ghost = gradW_ij * (-p_i->mass * ghost->mass * p_term);

            p_i->accel += f_ghost / p_i->mass;
            p_i->dudt += 0.5 * ghost->mass * p_term * v_ij.dot(gradW_ij);
        }
    }
}