#include "sphsolver.h"
#include "RigidBody.h"
#include <cmath>
#include <algorithm>
#include <atomic> // Required for std::atomic_ref (C++20)

SPHSolver::SPHSolver(Parameters params, std::shared_ptr<EquationOfState> eos) 
    :params_(params), eos_(std::move(eos))
{
}

void SPHSolver::addThreadPool(std::shared_ptr<ThreadPool> pool)
{
    threadPool_ = std::move(pool);
    threadPool_->start();
}

//double SPHSolver::computeSoundSpeed(double u, double gamma) const {
//    return std::sqrt(gamma * (gamma - 1.0) * std::max(u, 1e-6));
//}

void SPHSolver::updateSmoothingLengths(std::vector<std::shared_ptr<Particle>>& fluidParticles) const {
    //can mulitthread
    threadPool_->parallel_for(0, fluidParticles.size(), [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            auto& p = fluidParticles[i];
            double rho = std::max(p->density, 1e-5);
            p->h = params_.eta * std::sqrt(p->mass / rho);
        }
    });
}

void SPHSolver::computeDensityAndPressure(
    std::vector<std::shared_ptr<Particle>>& fluidParticles,
    const std::vector<std::shared_ptr<Particle>>& internalRigidParticles,
    const std::vector<std::shared_ptr<Particle>>& ghostParticles,
    const SpatialGrid& grid,
    const Kernel& kernel) const 
{

    threadPool_->parallel_for(0, fluidParticles.size(), [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            auto& p_i = fluidParticles[i];
            p_i->density = 0.0;
            double h_i = p_i->h;

            auto neighborIndices = grid.getNeighborIndices(p_i->pos);
            for (size_t idx : neighborIndices) {
                if (idx >= fluidParticles.size()) continue;
                const auto& p_j = fluidParticles[idx];
                Vector2D r_ij = p_i->pos - p_j->pos;
                p_i->density += p_j->mass * kernel.value(r_ij, h_i);
            }
        }
    });

    auto applyBoundaryDensity = [&](const std::vector<std::shared_ptr<Particle>>& boundaries) {
        threadPool_->parallel_for(0, boundaries.size(), [&](size_t start, size_t end) {
            for (size_t b = start; b < end; ++b) {
                const auto& bp = boundaries[b];
                auto neighbors = grid.getNeighborIndices(bp->pos);

                for (size_t idx : neighbors) {
                    if (idx >= fluidParticles.size()) continue;
                    auto& p_i = fluidParticles[idx];

                    Vector2D r_ij = p_i->pos - bp->pos;
                    double mass_contrib = bp->mass * kernel.value(r_ij, p_i->h);

                    if (mass_contrib > 0.0) {
                        std::atomic_ref<double> density_ref(p_i->density);
                        density_ref.fetch_add(mass_contrib, std::memory_order_relaxed);
                    }
                }
            }
        });
    };


    applyBoundaryDensity(internalRigidParticles);
    applyBoundaryDensity(ghostParticles);
    //threadPool_->waitFinished();

    threadPool_->parallel_for(0, fluidParticles.size(), [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            if (fluidParticles[i]->density < 1e-5) {
                fluidParticles[i]->density = 1e-5;
            }
        }
    });

    //can multithread 
    //for (auto& p_i : fluidParticles) {
    //    threadPool_->QueueJob([&](){
    //        p_i->density = 0.0;
    //        double h_i = p_i->h;
//
    //        // 1. Fluid-Fluid Density Contributions
    //        auto neighborIndices = grid.getNeighborIndices(p_i->pos);
    //        for (size_t idx : neighborIndices) {
    //            if (idx >= fluidParticles.size()) continue;
    //            const auto& p_j = fluidParticles[idx];
    //            Vector2D r_ij = p_i->pos - p_j->pos;
    //            p_i->density += p_j->mass * kernel.value(r_ij, h_i);
    //        }
//
    //        // 2. Internal Rigid Body Particle Contributions
    //        for (const auto& rp : internalRigidParticles) {
    //            Vector2D r_ij = p_i->pos - rp->pos;
    //            p_i->density += rp->mass * kernel.value(r_ij, h_i);
    //        }
//
    //        // 3. External Ghost Boundary Contributions
    //        //this is wrong
    //        for (const auto& ghost : ghostParticles) {
    //            Vector2D r_ij = p_i->pos - ghost->pos;
    //            p_i->density += ghost->mass * kernel.value(r_ij, h_i);
    //        }
//
    //        if (p_i->density < 1e-5) p_i->density = 1e-5;
    //    });
    //}
    //threadPool_->waitFinished();
}

void SPHSolver::computeDerivatives(
    std::vector<std::shared_ptr<Particle>>& fluidParticles,
    const std::vector<std::shared_ptr<Particle>>& internalRigidParticles,
    const std::vector<std::shared_ptr<Particle>>& ghostParticles,
    const SpatialGrid& grid,
    const Kernel& kernel) const
{
    size_t n = fluidParticles.size();

    threadPool_->parallel_for(0, n, [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            fluidParticles[i]->accel = params_.enableGravity ? params_.gravity : Vector2D(0.0, 0.0);
            fluidParticles[i]->dudt = 0.0;
        }
    });

    double alpha_u = params_.thermal;
    
    threadPool_->parallel_for(0, n, [&, alpha_u](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            auto& p_i = fluidParticles[i];
            double rho_i = p_i->density;
            double u_i   = p_i->u;
            double P_i   = eos_->computePressure(rho_i, u_i);
            double h_i   = p_i->h;
            double c_i   = eos_->computeSoundSpeed(rho_i, u_i);

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
                double P_j    = eos_->computePressure(rho_j, u_j);
                double c_j    = eos_->computeSoundSpeed(rho_j, u_j);

                double pi_ij = 0.0;
                double v_dot_r = v_ij.dot(r_ij);
                if (v_dot_r < 0.0) {
                    double rho_ij = 0.5 * (rho_i + rho_j);
                    double c_ij   = 0.5 * (c_i + c_j);
                    double mu_ij  = (h_ij * v_dot_r) / (r2 + 0.01 * h_ij * h_ij);
                    pi_ij = (-params_.alpha * c_ij * mu_ij + params_.beta * mu_ij * mu_ij) / rho_ij;
                }

                double p_term = ((P_i + P_j) / (rho_i * rho_j)) + pi_ij;
                Vector2D force_ij = gradW_ij * (-p_i->mass * p_j->mass * p_term);
                double thermal_ij = 0.5 * p_term * v_ij.dot(gradW_ij);

                double rho_ij = 0.5 * (rho_i + rho_j);
                double v_u_ij = std::sqrt(std::abs(P_i - P_j) / rho_ij);
                double cond_term = (alpha_u * v_u_ij / rho_ij) * (u_i - u_j) * 
                                   (r_ij.dot(gradW_ij) / (r2 + 0.01 * h_ij * h_ij));

                std::atomic_ref<double> ax_i(p_i->accel.x), ay_i(p_i->accel.y);
                std::atomic_ref<double> ax_j(p_j->accel.x), ay_j(p_j->accel.y);
                std::atomic_ref<double> du_i(p_i->dudt), du_j(p_j->dudt);

                ax_i.fetch_add(force_ij.x / p_i->mass, std::memory_order_relaxed);
                ay_i.fetch_add(force_ij.y / p_i->mass, std::memory_order_relaxed);
                ax_j.fetch_sub(force_ij.x / p_j->mass, std::memory_order_relaxed);
                ay_j.fetch_sub(force_ij.y / p_j->mass, std::memory_order_relaxed);

                du_i.fetch_add(p_j->mass * (thermal_ij + cond_term), std::memory_order_relaxed);
                du_j.fetch_add(p_i->mass * (thermal_ij - cond_term), std::memory_order_relaxed);
            }
        }
    });

    // Helper for boundary forces
    auto applyBoundaryForces = [&](const std::vector<std::shared_ptr<Particle>>& boundaries, bool isRigid) {
        threadPool_->parallel_for(0, boundaries.size(), [&](size_t start, size_t end) {
            for (size_t b = start; b < end; ++b) {
                const auto& bp = boundaries[b];
                auto neighbors = grid.getNeighborIndices(bp->pos);

                for (size_t idx : neighbors) {
                    if (idx >= n) continue;
                    auto& p_i = fluidParticles[idx];

                    Vector2D r_ij = p_i->pos - bp->pos;
                    double r2 = r_ij.normSq();
                    double h_i  = p_i->h;
                    double h_j  = bp->h;
                    double h_ij = 0.5 * (h_i + h_j);
                    Vector2D gradW_ij = 0.5 * (kernel.gradient(r_ij, h_i) + kernel.gradient(r_ij, h_j));

                    if (gradW_ij.normSq() < 1e-24) continue;

                    double rho_i = p_i->density;
                    double u_i   = p_i->u;
                    double P_i   = eos_->computePressure(rho_i, u_i);
                    double c_i   = eos_->computeSoundSpeed(rho_i, u_i);

                    Vector2D v_ij = p_i->vel - bp->vel;
                    double rho_j = rho_i; 
                    double P_j   = P_i; 
                    double u_j   = bp->u;  

                    double pi_ij = 0.0;
                    double v_dot_r = v_ij.dot(r_ij);
                    if (v_dot_r < 0.0) {
                        double rho_ij = rho_i;
                        double c_ij   = c_i;
                        double mu_ij  = (h_ij * v_dot_r) / (r2 + 0.01 * h_ij * h_ij);
                        pi_ij = (-params_.alpha * c_ij * mu_ij + params_.beta * mu_ij * mu_ij) / rho_ij;
                    }

                    double p_term = ((P_i + P_j) / (rho_i * rho_j)) + pi_ij;
                    Vector2D force = gradW_ij * (p_i->mass * bp->mass * p_term);
                    Vector2D force_on_fluid = -force;

                    // Thermal conduction with boundary
                    double rho_ij = 0.5 * (rho_i + rho_j);
                    double v_u_ij = std::sqrt(std::abs(P_i - P_j) / rho_ij);
                    double alpha_u = 0.5;
                    double cond_term = (alpha_u * v_u_ij / rho_ij) * (u_i - u_j) * 
                                       (r_ij.dot(gradW_ij) / (r2 + 0.01 * h_ij * h_ij));
                    
                    std::atomic_ref<double> ax_i(p_i->accel.x), ay_i(p_i->accel.y);
                    std::atomic_ref<double> du_i(p_i->dudt);

                    ax_i.fetch_add(force_on_fluid.x / p_i->mass, std::memory_order_relaxed);
                    ay_i.fetch_add(force_on_fluid.y / p_i->mass, std::memory_order_relaxed);

                    double thermal_work = 0.5 * p_term * v_ij.dot(gradW_ij);
                    du_i.fetch_add(bp->mass * (thermal_work + cond_term), std::memory_order_relaxed);

                    if (isRigid && bp->parentBody) {
                        bp->parentBody->addForceAtPosition(force, bp->pos);
                    }
                }
            }
        });
    };

    // 2. Fluid-Rigid Interactions
    applyBoundaryForces(internalRigidParticles, true);

    // 3. Fluid-Ghost Interactions
    applyBoundaryForces(ghostParticles, false);

}