//#include <iostream>
//
//#include "Vector2D.h"
//#include "Particle.h"
//#include "kernel.h"
//#include "SpatialGrid.h"
//#include "SystemVolatility.h"
//#include "VelocityVerlet.h"
//#include "GhostManager.h"
//#include "Planer.h"
//#include "ThreadPool.h"

#include <iostream>
#include <memory>
#include <iomanip>

#include "simulation/simulation.h"
#include "simulation/sphsolver.h"
#include "phys/Integrator/kernel.h"

#include "simulation/initial_conditions/dam_break.h"
#include "simulation/initial_conditions/piston_compression.h"
#include "simulation/initial_conditions/explosion.h"
#include "simulation/initial_conditions/sod_shock_tube.h"
#include "simulation/initial_conditions/thermal_equilibrium.h"

#include "utility/SystemAggregator.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>
#include <thread>
#include <algorithm>

//inline double computeSoundSpeed(double u, double gamma = 1.4) {
//    return std::sqrt(gamma * (gamma - 1.0) * std::max(u, 1e-6));
//}
//
//void computeDensityAndPressure(
//    std::vector<std::unique_ptr<Particle>>& particles,
//    const std::vector<std::unique_ptr<Particle>>& ghosts,
//    const SpatialGrid& grid,
//    const Kernel& kernel,
//    double h,
//    double gamma = 1.4) 
//{
//    // Combined particle view for neighbor interactions
//    std::vector<const Particle*> allParticles;
//    allParticles.reserve(particles.size() + ghosts.size());
//    for (const auto& p : particles) allParticles.push_back(p.get());
//    for (const auto& g : ghosts) allParticles.push_back(g.get());
//
//    for (auto& p_i : particles) {
//        p_i->density = 0.0;
//
//        // Query spatial grid cell neighbors
//        auto neighborIndices = grid.getNeighborIndices(p_i->pos);
//        for (size_t idx : neighborIndices) {
//            if (idx >= particles.size()) continue;
//            const auto& p_j = particles[idx];
//            Vector2D r_ij = p_i->pos - p_j->pos;
//            p_i->density += p_j->mass * kernel.value(r_ij, h);
//        }
//
//        // Add contributions from ghost particles
//        for (const auto& ghost : ghosts) {
//            Vector2D r_ij = p_i->pos - ghost->pos;
//            p_i->density += ghost->mass * kernel.value(r_ij, h);
//        }
//
//        // Prevent density underflow
//        if (p_i->density < 1e-5) p_i->density = 1e-5;
//    }
//}
//
//void computeDerivatives(
//    std::vector<std::unique_ptr<Particle>>& particles,
//    const std::vector<std::unique_ptr<Particle>>& ghosts,
//    const SpatialGrid& grid,
//    const Kernel& kernel,
//    double h,
//    double gamma = 1.4,
//    double alpha = 1.0,
//    double beta = 2.0)
//{
//    size_t n = particles.size();
//
//    for (size_t i = 0; i < n; ++i) {
//        auto& p_i = particles[i];
//        p_i->accel = Vector2D(0.0, 0.0);
//        p_i->dudt = 0.0;
//
//        if (!p_i->isDynamic()) continue;
//
//        double rho_i = p_i->density;
//        double u_i   = p_i->u;
//        double P_i   = (gamma - 1.0) * rho_i * u_i;
//        double c_i   = computeSoundSpeed(u_i, gamma);
//
//        // Interaction with fluid particles
//        auto neighborIndices = grid.getNeighborIndices(p_i->pos);
//        for (size_t idx : neighborIndices) {
//            if (idx >= n || idx == i) continue;
//            const auto& p_j = particles[idx];
//
//            Vector2D r_ij = p_i->pos - p_j->pos;
//            Vector2D gradW = kernel.gradient(r_ij, h);
//            if (gradW.normSq() < 1e-24) continue;
//
//            Vector2D v_ij = p_i->vel - p_j->vel;
//            double rho_j  = p_j->density;
//            double u_j    = p_j->u;
//            double P_j    = (gamma - 1.0) * rho_j * u_j;
//            double c_j    = computeSoundSpeed(u_j, gamma);
//
//            // Monaghan Artificial Viscosity
//            double pi_ij = 0.0;
//            double v_dot_r = v_ij.dot(r_ij);
//            if (v_dot_r < 0.0) {
//                double rho_ij = 0.5 * (rho_i + rho_j);
//                double c_ij   = 0.5 * (c_i + c_j);
//                double mu_ij  = (h * v_dot_r) / (r_ij.normSq() + 0.01 * h * h);
//                pi_ij = (-alpha * c_ij * mu_ij + beta * mu_ij * mu_ij) / rho_ij;
//            }
//
//            double p_term = (P_i / (rho_i * rho_i)) + (P_j / (rho_j * rho_j)) + pi_ij;
//
//            // Momentum: dv/dt = - sum m_j * P_term * gradW
//            p_i->accel -= gradW * (p_j->mass * p_term);
//
//            // Thermal Work: du/dt = 0.5 * sum m_j * P_term * (v_ij . gradW)
//            p_i->dudt += 0.5 * p_j->mass * p_term * v_ij.dot(gradW);
//        }
//
//        // Interaction with ghost boundary particles
//        for (const auto& ghost : ghosts) {
//            Vector2D r_ij = p_i->pos - ghost->pos;
//            Vector2D gradW = kernel.gradient(r_ij, h);
//            if (gradW.normSq() < 1e-24) continue;
//
//            Vector2D v_ij = p_i->vel - ghost->vel;
//            double rho_j  = ghost->density;
//            double u_j    = ghost->u;
//            double P_j    = (gamma - 1.0) * rho_j * u_j;
//            double c_j    = computeSoundSpeed(u_j, gamma);
//
//            double pi_ij = 0.0;
//            double v_dot_r = v_ij.dot(r_ij);
//            if (v_dot_r < 0.0) {
//                double rho_ij = 0.5 * (rho_i + rho_j);
//                double c_ij   = 0.5 * (c_i + c_j);
//                double mu_ij  = (h * v_dot_r) / (r_ij.normSq() + 0.01 * h * h);
//                pi_ij = (-alpha * c_ij * mu_ij + beta * mu_ij * mu_ij) / rho_ij;
//            }
//
//            double p_term = (P_i / (rho_i * rho_i)) + (P_j / (rho_j * rho_j)) + pi_ij;
//
//            p_i->accel -= gradW * (ghost->mass * p_term);
//            p_i->dudt += 0.5 * ghost->mass * p_term * v_ij.dot(gradW);
//        }
//    }
//}
//
//void exportResults(const std::vector<std::unique_ptr<Particle>>& particles, const std::string& filename, double gamma = 1.4) {
//    std::ofstream file(filename);
//    file << "x,y,vx,vy,rho,p,u\n";
//    for (const auto& p : particles) {
//        if (!p->isFluid()) continue;
//        double pressure = (gamma - 1.0) * p->density * p->u;
//        file << p->pos.x << "," << p->pos.y << ","
//             << p->vel.x << "," << p->vel.y << ","
//             << p->density << "," << pressure << "," << p->u << "\n";
//    }
//}
//
//int main2() {
//    std::cout << "=========================================================\n";
//    std::cout << "        SPH 1D Sod Shock Tube Baseline Verification      \n";
//    std::cout << "=========================================================\n";
//
//    // 1. Problem Setup Parameters
//    const double gamma = 1.4;
//    const double domainLength = 1.0;  // Domain x in [-0.5, 0.5]
//    const double xMin = -0.5;
//    const double xMax = 0.5;
//    const double tFinal = 0.05;       // Standard Sod analytical evaluation time
//    
//    // Sod Initial States
//    const double rhoL = 1.0,   PL = 1.0;
//    const double rhoR = 0.125, PR = 0.1;
//    const double uL = PL / ((gamma - 1.0) * rhoL); // uL = 2.5
//    const double uR = PR / ((gamma - 1.0) * rhoR); // uR = 2.0
//
//    const int numParticlesLeft = 400;
//    const double dxL = 0.5 / numParticlesLeft;
//    const double particleMass = rhoL * dxL;
//    const double dxR = particleMass / rhoR;
//    const int numParticlesRight = static_cast<int>(0.5 / dxR);
//
//    const double h = 0.015; // Smoothing length
//    int particleID = 0;
//
//    std::vector<std::unique_ptr<Particle>> particles;
//
//    // 2. Generate Left Fluid Domain
//    for (int i = 0; i < numParticlesLeft; ++i) {
//        double x = xMin + (i + 0.5) * dxL;
//        particles.push_back(std::make_unique<FluidParticle>(
//            particleID++, Vector2D(x, 0.0), Vector2D(0.0, 0.0), particleMass, rhoL, uL
//        ));
//    }
//
//    // 3. Generate Right Fluid Domain
//    for (int i = 0; i < numParticlesRight; ++i) {
//        double x = (i + 0.5) * dxR;
//        if (x >= xMax) break;
//        particles.push_back(std::make_unique<FluidParticle>(
//            particleID++, Vector2D(x, 0.0), Vector2D(0.0, 0.0), particleMass, rhoR, uR
//        ));
//    }
//
//    // 4. Attach Boundary Walls at End Points
//    auto leftWall = std::make_unique<SolidParticle>(
//        particleID++, Vector2D(xMin, 0.0), Vector2D(0.0, 0.0), particleMass, rhoL, uL, MotionType::STATIC
//    );
//    leftWall->setBoundaryComponent(std::make_unique<PlanarBoundary>(Vector2D(1.0, 0.0), /*isNoSlip=*/false));
//    particles.push_back(std::move(leftWall));
//
//    auto rightWall = std::make_unique<SolidParticle>(
//        particleID++, Vector2D(xMax, 0.0), Vector2D(0.0, 0.0), particleMass, rhoR, uR, MotionType::STATIC
//    );
//    rightWall->setBoundaryComponent(std::make_unique<PlanarBoundary>(Vector2D(-1.0, 0.0), /*isNoSlip=*/false));
//    particles.push_back(std::move(rightWall));
//
//    std::cout << "Initialized " << particles.size() << " total particles (Fluid + Boundaries).\n";
//
//    // 5. Initialize Subsystems
//    unsigned int hardwareThreads = std::thread::hardware_concurrency();
//    int numThreads = (hardwareThreads > 0) ? static_cast<int>(hardwareThreads) : 4;
//    ThreadPool threadPool(numThreads);
//    threadPool.start();
//
//    CubicSplineKernel kernel;
//    VelocityVerletIntegrator integrator(h, gamma, /*cflSafety=*/0.25, /*maxDt=*/0.001, /*minDt=*/1e-6, /*sampleSize=*/128);
//
//    SystemAggregator aggregator("energy_conservation.csv");
//
//    double dt = 0.0002;
//    double t = 0.0;
//    int step = 0;
//
//    // Initial grid construction & density/derivative pass
//    integrator.getSpatialGrid().build(particles);
//    auto ghosts = GhostManager::generateAllGhosts(particles, 2.0 * h);
//    computeDensityAndPressure(particles, ghosts, integrator.getSpatialGrid(), kernel, h, gamma);
//    computeDerivatives(particles, ghosts, integrator.getSpatialGrid(), kernel, h, gamma);
//
//    aggregator.processAndLog(particles, t, step);
//
//    std::cout << "Starting simulation loop up to t = " << tFinal << "s...\n";
//
//    // 6. Time Integration Loop
//    while (t < tFinal) {
//        // Kick 1 & Drift (Parallelized across ThreadPool)
//        integrator.kickFirstHalf(particles, dt, threadPool);
//        integrator.drift(particles, dt, threadPool);
//
//        // Rebuild Spatial Grid & Update Ghosts
//        integrator.getSpatialGrid().build(particles);
//        ghosts = GhostManager::generateAllGhosts(particles, 2.0 * h);
//
//        // Field Property Updates & Derivatives
//        computeDensityAndPressure(particles, ghosts, integrator.getSpatialGrid(), kernel, h, gamma);
//        computeDerivatives(particles, ghosts, integrator.getSpatialGrid(), kernel, h, gamma);
//
//        // Kick 2
//        integrator.kickSecondHalf(particles, dt, threadPool);
//
//        // Evaluate Volatility on Strided Sample & Adapt Timestep
//        dt = integrator.computeAdaptiveTimestep(particles, dt);
//
//        t += dt;
//        step++;
//
//        if (step % 10 == 0) {
//            SystemMetrics metrics = aggregator.processAndLog(particles, t, step);
//
//            if (step % 100 == 0) {
//                std::cout << "Step: " << step << " | Time: " << t << " / " << tFinal << "s | Next dt: " << dt << "s\n";
//                SystemAggregator::printSummary(metrics);
//            }
//        }
//    }
//
//    threadPool.Stop();
//
//    // 7. Output Final Results
//    exportResults(particles, "sod_shock_tube.csv", gamma);
//    std::cout << "Simulation complete. Output written to 'sod_shock_tube.csv'.\n";
//
//    return 0;
//}

int main() {
    std::cout << "========================================\n";
    std::cout << "  SPH FSI & Thermodynamics Simulator    \n";
    std::cout << "========================================\n\n";

    // 1. Configure SPH Solver Parameters
    Parameters solverParams;
    solverParams.gamma = 1.4; // Ratio of specific heats for ideal gas
    solverParams.alpha = 0.5; // Artificial viscosity alpha coefficient
    solverParams.beta  = 1.0; // Artificial viscosity beta coefficient
    solverParams.eta   = 1.3; // Scale factor for adaptive smoothing length (h_i)
    solverParams.thermal = 2.0; // Thermal conduction coefficient

    // 2. Instantiate Concrete Kernel (passed via shared_ptr for polymorphic Kernel base)
    auto kernel = std::make_shared<CubicSplineKernel>();

    //auto ic = std::make_unique<DamBreakIC>(damParams);
    //ThermalEquilibriumIC::Parameters TEp;
    //auto ic = std::make_unique<ThermalEquilibriumIC>(TEp);

    //DamBreakIC::Parameters damParams;
    //damParams.enableDamWallPivot = true;
    //damParams.gravity = Vector2D(0.0, -0.981);
    //damParams.damWallRestoringStiffness = 0.01;
    //damParams.damWallRestoringDamping = 0.002;
    //auto ic = std::make_unique<DamBreakIC>(damParams);

    enum class Scenario { PistonCompression, Explosion, PistonEngine, FlywheelEngine };
    const Scenario scenario = Scenario::FlywheelEngine; // Change this to select different scenarios

    PistonCompressionIC::Parameters damParams;
    damParams.driveMode = PistonCompressionIC::DriveMode::ConstantForce;
    ExplosionIC::Parameters blastParams;

    std::unique_ptr<InitialConditions> ic;
    double gasGamma = damParams.gamma;
    bool gravityEnabled = damParams.enableGravity;
    Vector2D gravityVector = damParams.gravity;
    Vector2D domainMin = damParams.domainMin, domainMax = damParams.domainMax;
    double tFinal = 1000.0;
    double initialDt = 1e-4;
    double frameRate = 12000.0; // exported frames per simulated second
    if (scenario == Scenario::Explosion) {
        ic = std::make_unique<ExplosionIC>(blastParams);
        gasGamma = blastParams.gamma;
        gravityEnabled = blastParams.enableGravity;
        gravityVector = blastParams.gravity;
        domainMin = blastParams.domainMin;
        domainMax = blastParams.domainMax;
        tFinal = 0.1;
        initialDt = 1e-8;   // the first step already sees a 2000x energy jump
        frameRate = 3000.0; // the blast crosses the box in a few ms
        solverParams.thermal = 0.5;
    } else if (scenario == Scenario::PistonEngine) {
        // A heavy load squeezes and heats the gas; the gas then pushes the piston back up.
        // The power take-off siphons work out of the piston on every stroke.
        damParams.constantForceFactor = 1.5;
        damParams.pistonCeilingAboveStart = 0.25;
        damParams.powerTakeOff = 0.15;
        ic = std::make_unique<PistonCompressionIC>(damParams);
        tFinal = 10.0;
        frameRate = 120.0;
        solverParams.thermal = 0.5;
    } else if (scenario == Scenario::FlywheelEngine) {
        // Slider-crank engine: the flywheel compresses the gas, heat added near the bottom of the stroke
        // raises the pressure, and the gas drives the flywheel round again. The shaft generator siphons work.
        damParams.driveMode = PistonCompressionIC::DriveMode::Flywheel;
        damParams.crankRadius = 0.12;
        damParams.flywheelInertia = 0.02;
        damParams.flywheelInitialOmega = 3.0;
        damParams.shaftDamping = 0.000;
        damParams.loadFactor =0.0;
        damParams.constantForceFactor = 0.0;
        ic = std::make_unique<PistonCompressionIC>(damParams);
        tFinal = 20.0;
        frameRate = 120.0;
        solverParams.thermal = 0.5;
    } else {
        ic = std::make_unique<PistonCompressionIC>(damParams);
    }

    solverParams.enableGravity = gravityEnabled;
    solverParams.gravity = gravityVector;

    size_t threadCount = 40;
    auto pool = std::make_shared<ThreadPool>(threadCount);
    //auto liquidEOS = std::make_shared<TaitEOS>(1.0 /*rho0*/, 20.0 /*soundSpeed0*/, 7.0 /*gamma*/);
    auto GAS = std::make_shared<IdealGasEOS>(gasGamma);
    SPHSolver solver(solverParams, GAS);
    solver.addThreadPool(pool);

    // 2. Configure Simulation Engine with Gravity
    SimulationEngine::Configuration config;
    //config.enableGravity = false;
    //config.gravity = Vector2D(0.0,0.0);
    //config.domainMin = TEp.domainMin;
    //config.domainMax = TEp.domainMax;
    config.enableGravity = gravityEnabled;
    config.gravity = gravityVector;
    config.domainMin = domainMin;
    config.domainMax = domainMax;
    config.tFinal = tFinal;
    config.initialDt = initialDt;

    // Heat bath touching the walls; colder than the gas so it removes heat and damps the piston bounce
    const bool useThermalBath = scenario != Scenario::Explosion;
    const double bathInternalEnergy = scenario == Scenario::PistonEngine ? damParams.gasInternalEnergy : 0.3;
    config.bath.enabled = useThermalBath;
    config.bath.targetU = bathInternalEnergy;
    config.bath.tau = 0.7;
    config.bath.wallThickness = 0.15;
    if (scenario == Scenario::FlywheelEngine) {
        // Cold sink most of the cycle, heat source only while the piston is near the bottom of its stroke
        config.bath.targetU = 0.4;
        config.bath.tau = 0.1;
        config.bath.hotU = 1.5;
        config.bath.hotBelowY = 0.6;//0.65;// - 2.0 * damParams.crankRadius + 0.06;
        std::cout<< "Bath hot region below y = " << config.bath.hotBelowY << "\n";
        config.bath.mode = BathMode::HOT_ONLY;
    }



    // 3. Create Engine
    SimulationEngine engine(config, solver, kernel, std::move(ic), threadCount);

    std::cout << "Simulation Initialized:\n";
    std::cout << "  - Fluid Particles : " << engine.getParticles().size() << "\n";
    std::cout << "  - Rigid Bodies    : " << engine.getRigidObjects().size() << "\n";
    std::cout << "  - Worker Threads  : " << threadCount << "\n\n";

    // 6. Time Integration Step Loop
    double logInterval = 0.05;  // Print output every 0.05 physical seconds
    double nextLogTime = 0.0;

    std::string runName = "piston_engine";
    if (scenario == Scenario::Explosion) runName = "explosion";
    else if (scenario == Scenario::FlywheelEngine) runName = "flywheel_engine";

    std::string metricsPath = "output/" + runName + "/metrics.csv";
    std::string framesDir = "output/" + runName + "/frames";
    std::string bathPath = "output/" + runName + "/bath.csv";

    SystemAggregator logger(metricsPath, framesDir);
    {
        // Read by animate.py to draw the bath region
       std::ofstream bathFile(bathPath);
        bathFile << "enabled,target_u,tau,wall_thickness,xmin,ymin,xmax,ymax,hot_u,hot_below_y,mode\n"
                 << config.bath.enabled << "," 
                 << config.bath.targetU << "," 
                 << config.bath.tau << ","
                 << config.bath.wallThickness << "," 
                 << config.domainMin.x << "," 
                 << config.domainMin.y << ","
                 << config.domainMax.x << "," 
                 << config.domainMax.y << ","
                 << config.bath.hotU << "," 
                 << config.bath.hotBelowY << "," 
                 << static_cast<int>(config.bath.mode) << "\n";
    }
    double nextFrameTime = 0.0;
    double frameInterval = 1.0 / frameRate;
    size_t frameIdx = 0;

    engine.getThreadPool().start();
    SystemMetrics metrics = logger.processAndLog(engine.getParticles(), engine.getRigidObjects(), 0.0, 0);
    logger.exportFrameCSV(frameIdx++, engine.getParticles(), engine.getRigidObjects());
    while (engine.getCurrentTime() < config.tFinal) {
        engine.step();

        double currentTime = engine.getCurrentTime();
        int currentStep = engine.getCurrentStep();

        // 1. Process and write CSV time-series metrics on every integration step
        metrics = logger.processAndLog(engine.getParticles(), engine.getRigidObjects(), currentTime, currentStep);

        // 2. Export Visualization Frame (Unnested: runs every step to catch precise 60 FPS target times)
        if (currentTime >= nextFrameTime) {
            logger.exportFrameCSV(frameIdx, engine.getParticles(), engine.getRigidObjects());
            logger.exportActiveGridCSV(frameIdx, engine.getSpatialGrid()); // Export C++ ground-truth grid
            logger.exportBathFrameCSV(frameIdx, currentTime, engine.getThermalBath());
            frameIdx++;
            nextFrameTime += frameInterval;
        }

        // 3. Console Terminal Output & Summary Print
        if (currentTime >= nextLogTime) {
            std::cout << "[Step " << std::setw(6) << currentStep << "] "
                      << "t = " << std::fixed << std::setprecision(4) << currentTime << "s"
                      << " | bath heat = " << std::setprecision(5) << engine.getBathHeat()
                      << " | bath active u = " << std::setprecision(5) << engine.getCurrentActiveU();
            for (const auto& body : engine.getRigidObjects()) {
                if (auto piston = std::dynamic_pointer_cast<Piston>(body)) {
                    std::cout << " | piston y = " << std::setprecision(3) << piston->getCenterOfMass().y
                              << " | harvested = " << std::setprecision(5) << piston->getHarvestedEnergy();
                    if (piston->hasFlywheel()) {
                        std::cout << " | omega = " << std::setprecision(3) << piston->getFlywheelOmega()
                                  << " | flywheel E = " << std::setprecision(5) << piston->getFlywheelEnergy();
                    }
                }
            }

            // Track dynamic rigid body Center of Mass if present in scene
            //if (!engine.getRigidObjects().empty()) {
            //    const auto& body = engine.getRigidObjects()[4];
            //    Vector2D com = body->getCenterOfMass();
            //    std::cout << " | Rigid Body CoM: (" 
            //              << std::setprecision(3) << com.x << ", " << com.y << ")";
            //}
            std::cout << "\n";

            if (currentStep % 100 == 0) {
                SystemAggregator::printSummary(metrics);
            }

            nextLogTime += logInterval;
        }
    }
    engine.getThreadPool().Stop();

    std::cout << "\nSimulation Complete! Final time = " 
              << engine.getCurrentTime() << "s\n";

    return 0;
}