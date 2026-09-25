#include <algorithm>
#include <cmath>
#include "VelocityVerlet.h"

void VelocityVerletIntegrator::kickFirstHalf(
    std::vector<std::shared_ptr<Particle>>& particles, 
    double dt, 
    ThreadPool& pool) const 
{
    int n = static_cast<int>(particles.size());
    int numThreads = pool.getMaxThreads();
    int chunkSize = (n + numThreads - 1) / numThreads;
    const double halfDt = 0.5 * dt;

    for (int t = 0; t < numThreads; ++t) {
        int start = t * chunkSize;
        int end = std::min(start + chunkSize, n);
        if (start >= n) break;

        pool.QueueJob([&particles, halfDt, start, end]() {
            for (int i = start; i < end; ++i) {
                auto& p = particles[i];
                if (!p->isDynamic()) continue;

                p->vel += p->accel * halfDt;
                p->u   += p->dudt  * halfDt;

                if (p->u < 1e-5) p->u = 1e-5;
            }
        }, t);
    }
    pool.waitFinished();
}

void VelocityVerletIntegrator::drift(
    std::vector<std::shared_ptr<Particle>>& particles, 
    double dt, 
    ThreadPool& pool) const 
{
    int n = static_cast<int>(particles.size());
    int numThreads = pool.getMaxThreads();
    int chunkSize = (n + numThreads - 1) / numThreads;

    for (int t = 0; t < numThreads; ++t) {
        int start = t * chunkSize;
        int end = std::min(start + chunkSize, n);
        if (start >= n) break;

        pool.QueueJob([&particles, dt, start, end]() {
            for (int i = start; i < end; ++i) {
                auto& p = particles[i];
                if (!p->isDynamic()) continue;

                p->pos += p->vel * dt;
            }
        }, t);
    }
    pool.waitFinished();
}

void VelocityVerletIntegrator::kickSecondHalf(
    std::vector<std::shared_ptr<Particle>>& particles, 
    double dt, 
    ThreadPool& pool) const 
{
    int n = static_cast<int>(particles.size());
    int numThreads = pool.getMaxThreads();
    int chunkSize = (n + numThreads - 1) / numThreads;
    const double halfDt = 0.5 * dt;

    for (int t = 0; t < numThreads; ++t) {
        int start = t * chunkSize;
        int end = std::min(start + chunkSize, n);
        if (start >= n) break;

        pool.QueueJob([&particles, halfDt, start, end]() {
            for (int i = start; i < end; ++i) {
                auto& p = particles[i];
                if (!p->isDynamic()) continue;

                p->vel += p->accel * halfDt;
                p->u   += p->dudt  * halfDt;

                if (p->u < 1e-5) p->u = 1e-5;
            }
        }, t);
    }
    pool.waitFinished();
}

double VelocityVerletIntegrator::computeAdaptiveTimestep(
    const std::vector<std::shared_ptr<Particle>>& particles,
    double currentDt) const 
{
    // Evaluates strided sample (sampleSize = 128) sequentially
    // Sequential execution is faster than thread dispatch overhead for 128 samples
    VolatilityMetrics metrics = SystemVolatility::evaluate(
        particles, h_, gamma_, sampleSize_
    );

    double dtTarget = maxDt_;
    if (metrics.compositeVolatility > 1e-8) {
        dtTarget = cflSafety_ / metrics.compositeVolatility;
    }

    if (currentDt > 0.0) {
        dtTarget = std::min(dtTarget, 1.2 * currentDt);
    }

    return std::clamp(dtTarget, minDt_, maxDt_);
}