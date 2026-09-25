#pragma once

#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <filesystem>

#include "Particle.h"
#include "Vector2D.h"

namespace fs = std::filesystem;
struct SystemMetrics {
    double time = 0.0;
    int step = 0;

    // Conservation quantities
    Vector2D totalMomentum{0.0, 0.0};
    double kineticEnergy = 0.0;
    double internalEnergy = 0.0;
    double totalEnergy = 0.0;
    double totalMass = 0.0;

    Vector2D fluidMomentum{0.0, 0.0};
    Vector2D rigidBodyMomentum{0.0, 0.0};
    double fluidKineticEnergy = 0.0;
    double rigidBodyKineticEnergy = 0.0;

    // Thermodynamic state & stability bounds
    double avgDensity = 0.0;
    double minDensity = 1e9;
    double maxDensity = -1e9;
    double minPressure = 1e9;
    double maxPressure = -1e9;
    double maxSpeed = 0.0;

    size_t particleCount = 0;
    size_t rigidBodyCount = 0;
};

class SystemAggregator {
private:
    std::ofstream metricsLogFile_;
    std::string frameOutputDir_ = "output/frames";
    bool metricsLoggingEnabled_ = false;

public:
    SystemAggregator() = default;
    explicit SystemAggregator(const std::string& logFilePath, 
                             const std::string& frameOutputDir = "output/frames");
    ~SystemAggregator();

    void setFrameDirectory(const std::string& frameOutputDir);
    bool openLogFile(const std::string& logFilePath);
    static void printSummary(const SystemMetrics& metrics);

    // --- Templated Methods (Header-defined for generic container instantiation) ---

    template <typename ParticleContainer>
    SystemMetrics compute(const ParticleContainer& particles, double time = 0.0, int step = 0) const {
        std::vector<int> emptyRigidObjects; // Dummy fallback
        return compute(particles, emptyRigidObjects, time, step);
    }

    // Compute combined metrics for particles AND rigid bodies
    template <typename ParticleContainer, typename RigidContainer>
    SystemMetrics compute(const ParticleContainer& particles, 
                          const RigidContainer& rigidObjects, 
                          double time = 0.0, 
                          int step = 0) const 
    {
        SystemMetrics metrics;
        metrics.time = time;
        metrics.step = step;

        double densitySum = 0.0;

        // 1. Process Fluid / General Particles
        for (const auto& p : particles) {
            Vector2D momentum = p->getMomentum();
            double kinE = p->getKineticEnergy();
            double intE = p->getInternalEnergy();

            metrics.fluidMomentum += momentum;
            metrics.fluidKineticEnergy += kinE;
            metrics.internalEnergy += intE;
            metrics.totalMass += p->mass;

            // Density bounds
            double d = p->density;
            densitySum += d;
            metrics.minDensity = std::min(metrics.minDensity, d);
            metrics.maxDensity = std::max(metrics.maxDensity, d);

            // Pressure bounds
            //double press = p->pressure;
            //metrics.minPressure = std::min(metrics.minPressure, press);
            //metrics.maxPressure = std::max(metrics.maxPressure, press);
//
            // Maximum speed tracking
            double speed = p->vel.length();
            metrics.maxSpeed = std::max(metrics.maxSpeed, speed);

            metrics.particleCount++;
        }

        if (metrics.particleCount > 0) {
            metrics.avgDensity = densitySum / static_cast<double>(metrics.particleCount);
        } else {
            metrics.minDensity = 0.0;
            metrics.maxDensity = 0.0;
            metrics.minPressure = 0.0;
            metrics.maxPressure = 0.0;
        }

        // 2. Process Rigid Bodies
        for (const auto& obj : rigidObjects) {
            metrics.rigidBodyCount++;
            for (const auto& p : obj->getParticles()) {
                metrics.rigidBodyMomentum += p->getMomentum();
                metrics.rigidBodyKineticEnergy += p->getKineticEnergy();
                metrics.totalMass += p->mass;
            }
        }

        // Combined Totals
        metrics.totalMomentum = metrics.fluidMomentum + metrics.rigidBodyMomentum;
        metrics.kineticEnergy = metrics.fluidKineticEnergy + metrics.rigidBodyKineticEnergy;
        metrics.totalEnergy = metrics.kineticEnergy + metrics.internalEnergy;

        return metrics;
    }

    template <typename ParticleContainer, typename RigidContainer>
    SystemMetrics processAndLog(const ParticleContainer& particles, 
                                const RigidContainer& rigidObjects, 
                                double time, 
                                int step) 
    {
        SystemMetrics metrics = compute(particles, rigidObjects, time, step);

        if (metricsLoggingEnabled_ && metricsLogFile_.is_open()) {
            metricsLogFile_ << std::scientific << std::setprecision(8)
                            << metrics.time << ","
                            << metrics.step << ","
                            << metrics.totalMomentum.x << ","
                            << metrics.totalMomentum.y << ","
                            << metrics.kineticEnergy << ","
                            << metrics.internalEnergy << ","
                            << metrics.totalEnergy << ","
                            << metrics.totalMass << ","
                            << metrics.avgDensity << ","
                            << metrics.minDensity << ","
                            << metrics.maxDensity << ","
                            << metrics.maxSpeed << ","
                            << metrics.particleCount << ","
                            << metrics.rigidBodyCount << "\n";
        }

        return metrics;
    }

    // Process and log particle-only metrics overload
    template <typename ParticleContainer>
    SystemMetrics processAndLog(const ParticleContainer& particles, double time, int step) {
        std::vector<int> emptyRigidObjects;
        return processAndLog(particles, emptyRigidObjects, time, step);
    }

    // Frame CSV export including both fluid particles and rigid body particles
    template <typename ParticleContainer>
    bool exportFrameCSV(size_t frameIndex, 
                        const ParticleContainer& particles) const 
    {
        if (frameOutputDir_.empty()) return false;

        std::filesystem::create_directories(frameOutputDir_);

        std::ostringstream filename;
        filename << frameOutputDir_ << "/frame_" 
                 << std::setfill('0') << std::setw(5) << frameIndex << ".csv";

        std::ofstream file(filename.str());
        if (!file.is_open()) return false;

        file << "id,x,y,vx,vy,density,internal_energy,mass,type\n";

        for (const auto& p : particles) {
            file << p->id << ","
                 << p->pos.x << "," << p->pos.y << ","
                 << p->vel.x << "," << p->vel.y << ","
                 << p->density << ","
                 << p->u << ","
                 << p->mass << ","
                 << (p->isDynamic() ? 0 : 1) << "\n";
        }

        return true;
    }
};