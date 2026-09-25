#include "sphsolver.h"
#include <cmath>
#include <algorithm>

SPHSolver::SPHSolver(Parameters params) : params_(params) {}

double SPHSolver::computeSoundSpeed(double u, double gamma) {
    return std::sqrt(gamma * (gamma - 1.0) * std::max(u, 1e-6));
}

void SPHSolver::updateSmoothingLengths(std::vector<std::shared_ptr<Particle>>& particles) const {
    for (auto& p : particles) {
        double rho = std::max(p->density, 1e-5);
        p->h = params_.eta * std::sqrt(p->mass / rho);
    }
}

void SPHSolver::computeDensityAndPressure(
    std::vector<std::shared_ptr<Particle>>& particles,
    const std::vector<std::shared_ptr<Particle>>& ghosts,
    const SpatialGrid& grid,
    const Kernel& kernel) const 
{
    for (auto& p_i : particles) {
        p_i->density = 0.0;
        double h_i = p_i->h;

        auto neighborIndices = grid.getNeighborIndices(p_i->pos);
        for (size_t idx : neighborIndices) {
            if (idx >= particles.size()) continue;
            const auto& p_j = particles[idx];
            Vector2D r_ij = p_i->pos - p_j->pos;
            p_i->density += p_j->mass * kernel.value(r_ij, h_i);
        }

        for (const auto& ghost : ghosts) {
            Vector2D r_ij = p_i->pos - ghost->pos;
            p_i->density += ghost->mass * kernel.value(r_ij, h_i);
        }

        if (p_i->density < 1e-5) p_i->density = 1e-5;
    }
}

void SPHSolver::computeDerivatives(
    std::vector<std::shared_ptr<Particle>>& particles,
    const std::vector<std::shared_ptr<Particle>>& ghosts,
    const SpatialGrid& grid,
    const Kernel& kernel,
    Vector2D& netBoundaryForce) const
{
    netBoundaryForce = Vector2D(0.0, 0.0);
    size_t n = particles.size();

    // Reset accelerations for ALL particles (fluid and rigid/solid boundary)
    for (auto& p : particles) {
        p->accel = Vector2D(0.0, 0.0);
        p->dudt = 0.0;
    }

    for (size_t i = 0; i < n; ++i) {
        auto& p_i = particles[i];
        if (!p_i->isDynamic()) continue; // Dynamic fluid particles drive physical interactions

        double rho_i = p_i->density;
        double u_i   = p_i->u;
        double P_i   = (params_.gamma - 1.0) * rho_i * u_i;
        double c_i   = computeSoundSpeed(u_i, params_.gamma);
        double h_i   = p_i->h;

        auto neighborIndices = grid.getNeighborIndices(p_i->pos);
        for (size_t idx : neighborIndices) {
            if (idx >= n || idx == i) continue;
            const auto& p_j = particles[idx];

            Vector2D r_ij = p_i->pos - p_j->pos;
            double r2 = r_ij.normSq();
            
            // Pairwise symmetrized smoothing length and kernel gradient
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
            
            // Hydrodynamic force vector exerted on particle i by particle j
            Vector2D force_ij = gradW_ij * (-p_i->mass * p_j->mass * p_term);

            // Apply force to fluid particle i
            p_i->accel += force_ij / p_i->mass;
            p_i->dudt += 0.5 * p_j->mass * p_term * v_ij.dot(gradW_ij);

            // Equal & Opposite force applied to particle j (read by RigidObject::accumulateForces)
            if (!p_j->isFluid()) {
                p_j->accel -= force_ij / p_j->mass;
            }
        }

        // Ghost Boundary interactions
        for (const auto& ghost : ghosts) {
            Vector2D r_ij = p_i->pos - ghost->pos;
            double h_ij = 0.5 * (h_i + ghost->h);
            Vector2D gradW_ij = 0.5 * (kernel.gradient(r_ij, h_i) + kernel.gradient(r_ij, ghost->h));

            if (gradW_ij.normSq() < 1e-24) continue;

            Vector2D v_ij = p_i->vel - ghost->vel;
            double rho_j  = ghost->density;
            double u_j    = ghost->u;
            double P_j    = (params_.gamma - 1.0) * rho_j * u_j;
            double c_j    = computeSoundSpeed(u_j, params_.gamma);

            double pi_ij = 0.0;
            double v_dot_r = v_ij.dot(r_ij);
            if (v_dot_r < 0.0) {
                double rho_ij = 0.5 * (rho_i + rho_j);
                double c_ij   = 0.5 * (c_i + c_j);
                double mu_ij  = (h_ij * v_dot_r) / (r_ij.normSq() + 0.01 * h_ij * h_ij);
                pi_ij = (-params_.alpha * c_ij * mu_ij + params_.beta * mu_ij * mu_ij) / rho_ij;
            }

            double p_term = (P_i / (rho_i * rho_i)) + (P_j / (rho_j * rho_j)) + pi_ij;
            Vector2D f_ghost = gradW_ij * (-p_i->mass * ghost->mass * p_term);

            p_i->accel += f_ghost / p_i->mass;
            p_i->dudt += 0.5 * ghost->mass * p_term * v_ij.dot(gradW_ij);

            netBoundaryForce += f_ghost;
        }
    }
}