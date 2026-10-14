#include "test_Explosion.h"
#include "explosion.h"
#include "SpatialGrid.h"
#include "ThreadPool.h"
#include "GhostManager.h"
#include <algorithm>
#include <cmath>

static bool testBlastRegionAndEnergy() {
    ExplosionIC::Parameters params;
    ExplosionIC ic(params);
    auto particles = ic.generateParticles();

    size_t hot = 0, ambient = 0;
    double maxU = 0.0;
    for (const auto& p : particles) {
        if (!p->isFluid()) continue;
        maxU = std::max(maxU, p->u);
        const double r = std::hypot(p->pos.x - params.blastCentre.x, p->pos.y - params.blastCentre.y);
        if (p->u > params.ambientInternalEnergy) {
            ++hot;
            TEST_ASSERT(r <= params.blastRadius + 1e-9, "Only particles inside the blast radius should be hot");
        } else {
            ++ambient;
        }
    }
    TEST_ASSERT(hot > 10, "The blast region should contain particles");
    TEST_ASSERT(ambient > 10 * hot, "Most of the gas should be ambient");
    TEST_ASSERT_NEAR(maxU, params.ambientInternalEnergy * params.blastEnergyRatio, 1e-9, "Blast energy ratio not applied");

    const double mass = params.gasDensity * params.spacing * params.spacing;
    const double expected = hot * mass * params.ambientInternalEnergy * (params.blastEnergyRatio - 1.0);
    TEST_ASSERT_NEAR(ic.getBlastEnergy(), expected, 1e-9, "Released blast energy should match the hot particle count");
    return true;
}

static bool testClosedBoxGhosts() {
    ExplosionIC ic;
    auto all = ic.generateParticles();
    std::vector<std::shared_ptr<Particle>> fluid, walls;
    for (const auto& p : all) {
        if (p->isFluid()) fluid.push_back(p);
        else if (p->hasBoundary()) walls.push_back(p);
    }

    // Place a gas particle next to each corner so every wall and corner host has something to mirror
    const double h = ic.getParameters().spacing;
    const Vector2D min = ic.getParameters().domainMin, max = ic.getParameters().domainMax;
    std::vector<std::shared_ptr<Particle>> probes;
    for (const Vector2D& pos : {Vector2D(min.x + h, min.y + h), Vector2D(max.x - h, min.y + h),
                                Vector2D(min.x + h, max.y - h), Vector2D(max.x - h, max.y - h)}) {
        probes.push_back(std::make_shared<FluidParticle>(-1, pos, Vector2D(0.0, 0.0), 1.0, 1.0, 1.0));
    }

    SpatialGrid grid(0.2);
    grid.build(probes);
    ThreadPool pool(2);
    pool.start();
    auto ghosts = GhostManager::generateAllGhosts(walls, probes, 0.2, grid, pool);
    pool.Stop();

    // Each corner probe sees two wall images and one diagonal corner image
    TEST_ASSERT(ghosts.size() == 12, "Each corner particle should have 3 ghosts (2 walls + 1 corner image)");
    return true;
}

void registerExplosionTests(TestSuite& suite) {
    suite.startSection("Explosion Module");
    suite.runTest("Blast Region and Energy", testBlastRegionAndEnergy);
    suite.runTest("Closed Box Ghosts", testClosedBoxGhosts);
}
