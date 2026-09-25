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
    // Step 1: Compute densities for all particles
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

    // Step 2: Compute pressures (Fluid gets ideal gas, Wall gets average neighbor fluid pressure)
    for (auto& p_i : particles) {
        if (p_i->isFluid()) {
            p_i->pressure = (params_.gamma - 1.0) * p_i->density * p_i->u;
        } else {
            // Wall / Boundary: copy local average fluid pressure to avoid vacuum suction
            double sumPressure = 0.0;
            int count = 0;
            auto neighborIndices = grid.getNeighborIndices(p_i->pos);
            for (size_t idx : neighborIndices) {
                if (idx < particles.size() && particles[idx]->isFluid()) {
                    sumPressure += (params_.gamma - 1.0) * particles[idx]->density * particles[idx]->u;
                    count++;
                }
            }
            p_i->pressure = (count > 0) ? (sumPressure / count) : 0.0;
        }
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

    for (auto& p : particles) {
        p->accel = Vector2D(0.0, -9.81); // Apply constant gravity
        p->dudt = 0.0;
    }

    for (size_t i = 0; i < n; ++i) {
        auto& p_i = particles[i];
        if (!p_i->isDynamic()) continue; // Only dynamic particles update acceleration

        double rho_i = p_i->density;
        double P_i   = p_i->isFluid() ? (params_.gamma - 1.0) * rho_i * p_i->u : p_i->pressure;
        double c_i   = computeSoundSpeed(p_i->u, params_.gamma);
        double h_i   = p_i->h;

        auto neighborIndices = grid.getNeighborIndices(p_i->pos);
        for (size_t idx : neighborIndices) {
            if (idx >= n || idx == i) continue;
            const auto& p_j = particles[idx];

            Vector2D r_ij = p_i->pos - p_j->pos;
            double r2 = r_ij.normSq();
            
            double h_j  = p_j->h;
            double h_ij = 0.5 * (h_i + h_j);
            Vector2D gradW_ij = 0.5 * (kernel.gradient(r_ij, h_i) + kernel.gradient(r_ij, h_j));

            if (gradW_ij.normSq() < 1e-24) continue;

            Vector2D v_ij = p_i->vel - p_j->vel;
            double rho_j  = p_j->density;
            double P_j    = p_j->isFluid() ? (params_.gamma - 1.0) * rho_j * p_j->u : p_j->pressure;
            double c_j    = computeSoundSpeed(p_j->u, params_.gamma);

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

            p_i->accel += force_ij / p_i->mass;
            
            if (p_i->isFluid()) {
                p_i->dudt += 0.5 * p_j->mass * p_term * v_ij.dot(gradW_ij);
            }

            if (!p_j->isFluid() && p_j->isDynamic()) {
                p_j->accel -= force_ij / p_j->mass;
            }
        }
    }
}