# SPH FSI & Thermodynamics Simulator

## Project Overview

This is a **Smoothed Particle Hydrodynamics (SPH)** simulation engine designed for **Fluid-Structure Interaction (FSI)** and **thermodynamic simulations**. The project is written in C++ and implements a mesh-free Lagrangian method for simulating compressible fluid dynamics coupled with rigid body mechanics.

---

## Architecture

### Directory Structure

```
├── src/
│   ├── main.cpp                    # Entry point with scenario configuration
│   ├── objects/                    # Particle and rigid body definitions
│   │   ├── Particle.h              # Fluid/Solid particle classes
│   │   ├── RigidBody.h             # Rigid object base class
│   │   ├── Piston.h/.cpp           # Specialized piston rigid body
│   │   └── ...
│   ├── phys/                       # Physics components
│   │   ├── EquationOfState.h       # EOS interfaces (Ideal Gas, Tait)
│   │   ├── Vector2D.h/.cpp         # 2D vector math
│   │   ├── IBoundary.h             # Boundary condition interface
│   │   ├── GhostManager.h          # Ghost particle generation
│   │   ├── ThermalZone.h           # Thermal region definitions
│   │   ├── Boundaries/             # Boundary implementations
│   │   └── Integrator/             # Time integration schemes
│   │       ├── kernel.h            # SPH kernel functions
│   │       └── VelocityVerlet.h    # Verlet integrator
│   ├── simulation/                 # Core simulation engine
│   │   ├── simulation.h/.cpp       # Main simulation engine
│   │   ├── sphsolver.h/.cpp        # SPH solver implementation
│   │   ├── rigidbodysolver.h/.cpp  # Rigid body dynamics solver
│   │   ├── initialconditions.h     # IC base class
│   │   ├── initial_conditions/     # Scenario implementations
│   │   │   ├── dam_break.h
│   │   │   ├── piston_compression.h
│   │   │   ├── explosion.h
│   │   │   ├── sod_shock_tube.h
│   │   │   └── thermal_equilibrium.h
│   │   └── thermalbath.h           # Berendsen-style heat bath
│   └── utility/                    # Utilities
│       ├── ThreadPool.h/.cpp       # Parallel computation
│       ├── SystemAggregator.h/.cpp # Metrics collection & export
│       ├── SpatialGrid.h           # Neighbor search acceleration
│       └── parser.h/.cpp           # Input file parsing
├── tests/                          # Unit tests
└── output/                         # Simulation results (CSV frames, metrics)
```

---

## Core Components

### 1. Particle System (`src/objects/Particle.h`)

Particles are the fundamental simulation elements with properties:
- **Position, Velocity, Acceleration** (`Vector2D`)
- **Mass, Density, Pressure**
- **Specific Internal Energy** (`u`) and its rate (`dudt`)
- **Smoothing Length** (`h`)

**Particle Types:**
| Type | MotionType | Description |
|------|------------|-------------|
| `FluidParticle` | DYNAMIC | Compressible fluid elements |
| `SolidParticle` | STATIC/NONSTATIC/DYNAMIC | Boundary walls, obstacles |

### 2. Equation of State (`src/phys/EquationOfState.h`)

Abstract interface for computing pressure from density and internal energy:

```cpp
class EquationOfState {
    virtual double computePressure(double density, double internalEnergy) const = 0;
    virtual double computeSoundSpeed(double density, double internalEnergy) const = 0;
    virtual double computeTemperature(double internalEnergy) const = 0;
};
```

**Implementations:**
- **`IdealGasEOS`**: `P = (γ-1)·ρ·u` for gases (γ = 1.4 default)
- **`TaitEOS`**: `P = B·((ρ/ρ₀)^γ - 1)` for liquids like water

### 3. SPH Solver (`src/simulation/sphsolver.h`)

Computes fluid dynamics using SPH formulations:

**Key Methods:**
- `computeDensityAndPressure()` - Density summation with ghost particles
- `computeDerivatives()` - Acceleration and internal energy rate via SPH interpolation
- `updateSmoothingLengths()` - Adaptive smoothing length (h_i = η·(m_i/ρ_i)^(1/d))

**SPH Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `gamma` | 1.4 | Ratio of specific heats |
| `alpha` | 1.0 | Artificial viscosity (linear) |
| `beta` | 1.0 | Artificial viscosity (quadratic) |
| `eta` | 1.3 | Smoothing length scale factor |
| `thermal` | 0.5 | Thermal conduction coefficient |

### 4. Rigid Bodies (`src/objects/RigidBody.h`)

Rigid objects composed of constituent particles:

**Features:**
- Linear and angular dynamics
- Pivot points with angular springs
- Constraint locking (X, Y, rotation)
- AABB collision bounds
- Force accumulation from fluid interactions

**Specialized: `Piston`**
- Direction locks (e.g., upward-only for gas-driven pistons)
- Constant load (weight/external force)
- Power take-off (energy harvesting via damping)
- Flywheel/slider-crank mechanism
- Travel limits and seals

### 5. Time Integration (`src/phys/Integrator/VelocityVerlet.h`)

**Velocity Verlet** scheme with adaptive timestepping:

```
1. Kick (first half):  v ← v + a·(dt/2)
2. Drift:              x ← x + v·dt
3. Compute forces
4. Kick (second half): v ← v + a·(dt/2)
```

**Adaptive Timestep:** Based on system volatility sampling (128-particle stride)

**CFL Condition:** `dt ≤ cflSafety · h / c_s` (sound speed limit)

### 6. SPH Kernels (`src/phys/Integrator/kernel.h`)

Smoothing kernels for particle interpolation:

- **`CubicSplineKernel`**: Standard C² continuous kernel (cutoff = 2h)
- **`WendlandC2Kernel`**: Compact support, positive definite

### 7. Thermal Bath (`src/simulation/thermalbath.h`)

**Berendsen-style thermostat** for temperature control:

```cpp
du/dt += (targetU - u) / τ
```

**Features:**
- Wall-proximity coupling (only particles near boundaries)
- Mode switching: `DEFAULT`, `HOT_ONLY`, `COLD_ONLY`
- Piston-position-dependent heating (engine cycles)

### 8. Simulation Engine (`src/simulation/simulation.h`)

Orchestrates the simulation loop:

```cpp
while (currentTime < tFinal) {
    engine.step();  // One integration step
    logger.processAndLog();  // Metrics
    logger.exportFrameCSV(); // Visualization frames
}
```

### 9. Thread Pool (`src/utility/ThreadPool.h`)

Parallel computation for:
- Particle integration (kick/drift)
- Density/derivative computation
- Spatial grid operations

---

## Simulation Scenarios

Configurable via `main.cpp` scenario enum:

| Scenario | Description |
|----------|-------------|
| `PistonCompression` | Gas compression by piston (constant velocity/force) |
| `Explosion` | High-energy blast wave propagation |
| `PistonEngine` | Single-acting heat engine with power take-off |
| `FlywheelEngine` | Slider-crank engine with flywheel energy storage |

**Initial Condition Modules:**
- `DamBreakIC` - Gravity-driven fluid collapse
- `SodShockTubeIC` - 1D shock tube Riemann problem
- `ThermalEquilibriumIC` - Temperature relaxation test
- `ExplosionIC` - Central high-pressure region
- `PistonCompressionIC` - Piston-cylinder gas dynamics

---

## Output & Visualization

### Metrics Logging (`output/metrics.csv`)

| Column | Description |
|--------|-------------|
| `time`, `step` | Simulation timestamp |
| `totalMomentum_x`, `totalMomentum_y` | Conservation check |
| `kineticEnergy`, `internalEnergy`, `totalEnergy` | Energy tracking |
| `avgDensity`, `minDensity`, `maxDensity` | Density bounds |
| `maxSpeed` | CFL timestep driver |
| `particleCount`, `rigidBodyCount` | System size |

### Frame Export (`output/frames/frame_XXXXX.csv`)

Per-particle snapshot for visualization:
```csv
id,x,y,vx,vy,density,internal_energy,mass,type
```
- `type=0`: Fluid particles
- `type=1`: Rigid body particles

### Bath Configuration (`output/bath.csv`)

Thermal bath parameters for animation overlay.

---

## Key Physics Features

### Artificial Viscosity (Monaghan)

```cpp
π_ij = (-α·c_ij·μ_ij + β·μ_ij²) / ρ_ij
```
Where `μ_ij = (h·v_ij·r_ij) / (|r_ij|² + ε·h²)`

### Ghost Particles

Boundary conditions enforced via mirror particles generated by `GhostManager`:
- No-slip/slip walls
- Pressure reflection
- Thermal boundary conditions

### Fluid-Structure Coupling

1. Fluid particles exert pressure forces on rigid body constituent particles
2. Rigid bodies move according to accumulated forces/torques
3. Ghost particles regenerated from updated boundary positions
4. Fluid particles constrained from penetrating solid boundaries

---

## Usage Example

```cpp
// Configure solver
Parameters solverParams;
solverParams.gamma = 1.4;
solverParams.alpha = 0.5;
solverParams.thermal = 0.5;

auto kernel = std::make_shared<CubicSplineKernel>();
auto EOS = std::make_shared<IdealGasEOS>(1.4);
SPHSolver solver(solverParams, EOS);

// Configure scenario
PistonCompressionIC::Parameters params;
params.driveMode = PistonCompressionIC::DriveMode::Flywheel;
params.crankRadius = 0.12;
params.flywheelInertia = 0.02;
auto ic = std::make_unique<PistonCompressionIC>(params);

// Configure thermal bath
config.bath.enabled = true;
config.bath.targetU = 0.4;
config.bath.mode = BathMode::HOT_ONLY;

// Run simulation
SimulationEngine engine(config, solver, kernel, std::move(ic), 40);
engine.run();
```

---

## Testing

Unit tests in `tests/` directory:
- `test_Vector2D.cpp` - Math primitives
- `test_RigidBody.cpp` - Rigid body dynamics
- `test_Piston.cpp` - Piston mechanics
- `test_Explosion.cpp` - Blast wave validation

---

## Dependencies

- **C++17** (filesystem, memory, threading)
- **Standard Library** only (no external dependencies)

---

## Performance Considerations

1. **Spatial Hashing**: O(1) neighbor lookup via `SpatialGrid`
2. **Thread Pool**: Configurable worker count (default: 40 threads)
3. **Adaptive Timestep**: Prevents instability during high-velocity events
4. **Strided Sampling**: Volatility computed on 128-particle sample for efficiency

---

## Future Extensions (Commented Code in main.cpp)

- GPU acceleration (KernelType dispatcher)
- Multi-phase flows (TaitEOS for liquids)
- Advanced boundary conditions (PlanarBoundary, etc.)
- Input file parsing (`parser.h`, `test.inp`)

---

## Summary

This is a **research-grade SPH simulation framework** for:
- ✅ Compressible gas dynamics
- ✅ Fluid-structure interaction
- ✅ Thermodynamic cycles (engines, compressors)
- ✅ Shock wave propagation
- ✅ Energy conservation tracking
- ✅ Parallel computation

The codebase is well-structured with clear separation between physics (SPH, EOS), numerics (integrators, kernels), and application logic (scenarios, I/O).

</content>