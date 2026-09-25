#include "SystemAggregator.h"
#include <iostream>
#include <iomanip>

namespace fs = std::filesystem;

SystemAggregator::SystemAggregator(const std::string& logFilePath, const std::string& frameOutputDir)
    : frameOutputDir_(frameOutputDir) 
{
    openLogFile(logFilePath);
    setFrameDirectory(frameOutputDir);
}

SystemAggregator::~SystemAggregator() {
    if (metricsLogFile_.is_open()) {
        metricsLogFile_.close();
    }
}

void SystemAggregator::setFrameDirectory(const std::string& frameOutputDir) {
    frameOutputDir_ = frameOutputDir;
    if (!frameOutputDir_.empty()) {
        fs::create_directories(frameOutputDir_);
    }
}

bool SystemAggregator::openLogFile(const std::string& logFilePath) {
    fs::path p(logFilePath);
    if (p.has_parent_path()) {
        fs::create_directories(p.parent_path());
    }

    metricsLogFile_.open(logFilePath);
    if (metricsLogFile_.is_open()) {
        metricsLoggingEnabled_ = true;
        metricsLogFile_ << "time,step,px,py,e_kin,e_int,e_tot,total_mass,"
                        << "avg_density,min_density,max_density,max_speed,n_particles,n_rigid_bodies\n";
        return true;
    }
    metricsLoggingEnabled_ = false;
    return false;
}

void SystemAggregator::printSummary(const SystemMetrics& metrics) {
    std::cout << std::fixed << std::setprecision(5)
              << "--- [Step " << metrics.step << " | Time: " << metrics.time << "s] ---\n"
              << "  Total Energy:        " << metrics.totalEnergy 
              << " (Kin: " << metrics.kineticEnergy 
              << " [Fluid: " << metrics.fluidKineticEnergy << " | Rigid: " << metrics.rigidBodyKineticEnergy << "]"
              << " | Int: " << metrics.internalEnergy << ")\n"
              << "  Total Momentum:      (" << metrics.totalMomentum.x << ", " << metrics.totalMomentum.y << ")\n"
              << "  Total Mass:          " << metrics.totalMass << "\n"
              << "  Fluid Density:       [" << metrics.minDensity << " - " << metrics.maxDensity << "] (Avg: " << metrics.avgDensity << ")\n"
              << "  Active Rigid Bodies: " << metrics.rigidBodyCount << "\n"
              << "----------------------------------------\n";
}