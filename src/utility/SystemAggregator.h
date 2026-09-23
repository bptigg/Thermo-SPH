#pragma once

#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "Particle.h"
#include "Vector2D.h"

// Snapshot container for aggregated global metrics
struct SystemMetrics {
    double time = 0.0;
    int step = 0;
    Vector2D totalMomentum{0.0, 0.0};
    double kineticEnergy = 0.0;
    double internalEnergy = 0.0;
    double totalEnergy = 0.0;
    double totalMass = 0.0;
    double avgDensity = 0.0;
    size_t particleCount = 0;
};

class SystemAggregator {
private:
    std::ofstream logFile_;
    bool loggingEnabled_ = false;

public:
    SystemAggregator() = default;

    // Optional constructor that opens a CSV file for time-series logging
    explicit SystemAggregator(const std::string& logFilePath) {
        openLogFile(logFilePath);
    }

    ~SystemAggregator() {
        if (logFile_.is_open()) {
            logFile_.close();
        }
    }

    // Opens log file and writes CSV header
    bool openLogFile(const std::string& logFilePath) {
        logFile_.open(logFilePath);
        if (logFile_.is_open()) {
            loggingEnabled_ = true;
            logFile_ << "time,step,px,py,e_kin,e_int,e_tot,total_mass,avg_density\n";
            return true;
        }
        loggingEnabled_ = false;
        return false;
    }

    // Computes system sums over std::unique_ptr<Particle> containers
    SystemMetrics compute(const std::vector<std::unique_ptr<Particle>>& particles, 
                          double time = 0.0, 
                          int step = 0) const {
        SystemMetrics metrics;
        metrics.time = time;
        metrics.step = step;

        double densitySum = 0.0;

        for (const auto& p : particles) {
            metrics.totalMomentum += p->getMomentum();
            metrics.kineticEnergy += p->getKineticEnergy();
            metrics.internalEnergy += p->getInternalEnergy();
            metrics.totalMass += p->mass;
            densitySum += p->density;
            metrics.particleCount++;
        }

        metrics.totalEnergy = metrics.kineticEnergy + metrics.internalEnergy;

        if (metrics.particleCount > 0) {
            metrics.avgDensity = densitySum / static_cast<double>(metrics.particleCount);
        }

        return metrics;
    }

    // Overload for std::shared_ptr<Particle> containers (e.g. RigidObjects/Pistons)
    SystemMetrics compute(const std::vector<std::shared_ptr<Particle>>& particles, 
                          double time = 0.0, 
                          int step = 0) const {
        SystemMetrics metrics;
        metrics.time = time;
        metrics.step = step;

        double densitySum = 0.0;

        for (const auto& p : particles) {
            metrics.totalMomentum += p->getMomentum();
            metrics.kineticEnergy += p->getKineticEnergy();
            metrics.internalEnergy += p->getInternalEnergy();
            metrics.totalMass += p->mass;
            densitySum += p->density;
            metrics.particleCount++;
        }

        metrics.totalEnergy = metrics.kineticEnergy + metrics.internalEnergy;

        if (metrics.particleCount > 0) {
            metrics.avgDensity = densitySum / static_cast<double>(metrics.particleCount);
        }

        return metrics;
    }

    // Calculates current system metrics and writes a line to the log file
    SystemMetrics processAndLog(const std::vector<std::unique_ptr<Particle>>& particles, 
                                double time, 
                                int step) {
        SystemMetrics metrics = compute(particles, time, step);

        if (loggingEnabled_ && logFile_.is_open()) {
            logFile_ << std::scientific << std::setprecision(8)
                     << metrics.time << ","
                     << metrics.step << ","
                     << metrics.totalMomentum.x << ","
                     << metrics.totalMomentum.y << ","
                     << metrics.kineticEnergy << ","
                     << metrics.internalEnergy << ","
                     << metrics.totalEnergy << ","
                     << metrics.totalMass << ","
                     << metrics.avgDensity << "\n";
        }

        return metrics;
    }

    // Formatted terminal output for step monitoring
    static void printSummary(const SystemMetrics& metrics) {
        std::cout << std::fixed << std::setprecision(5)
                  << "--- [Step " << metrics.step << " | Time: " << metrics.time << "s] ---\n"
                  << "  Total Energy:    " << metrics.totalEnergy 
                  << " (Kin: " << metrics.kineticEnergy << " | Int: " << metrics.internalEnergy << ")\n"
                  << "  Total Momentum:  (" << metrics.totalMomentum.x << ", " << metrics.totalMomentum.y << ")\n"
                  << "  Total Mass:      " << metrics.totalMass << "\n"
                  << "  Avg Density:     " << metrics.avgDensity << "\n"
                  << "----------------------------------------\n";
    }
};