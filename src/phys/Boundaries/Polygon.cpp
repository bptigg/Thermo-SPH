#include "Polygon.h"
#include <algorithm>
#include <cmath>
#include "Particle.h"

PolygonBoundary::PolygonBoundary(const std::vector<Vector2D>& localVertices, bool isNoSlip)
    : localVertices_(localVertices), isNoSlip_(isNoSlip) {}

std::vector<std::shared_ptr<Particle>> PolygonBoundary::generateGhosts(
    const Vector2D& hostPos,
    const Vector2D& hostVel,
    double hostEnergy,
    const std::vector<std::shared_ptr<Particle>>& fluidParticles,
    double supportRadius) const 
{
    std::vector<std::shared_ptr<Particle>> ghosts;
    int ghostID = -1;

    size_t numVerts = localVertices_.size();
    if (numVerts < 2) return ghosts;

    for (const auto& fluid : fluidParticles) {
        if (!fluid->isFluid()) continue;

        double minDistSq = supportRadius * supportRadius;
        Vector2D bestNormal(0.0, 0.0);
        double bestDist = -1.0;
        bool foundEdge = false;

        // Iterate over all edges defined by adjacent vertices
        for (size_t i = 0; i < numVerts; ++i) {
            Vector2D pA = hostPos + localVertices_[i];
            Vector2D pB = hostPos + localVertices_[(i + 1) % numVerts];

            Vector2D edge = pB - pA;
            double edgeLenSq = edge.normSq();
            if (edgeLenSq < 1e-12) continue;

            // Project fluid particle onto the edge segment parameter t in [0, 1]
            double t = std::clamp((fluid->pos - pA).dot(edge) / edgeLenSq, 0.0, 1.0);
            Vector2D closestPoint = pA + edge * t;

            Vector2D diff = fluid->pos - closestPoint;
            double distSq = diff.normSq();

            if (distSq < minDistSq) {
                minDistSq = distSq;
                bestDist = std::sqrt(distSq);

                // Outward normal for counter-clockwise vertices: (edge.y, -edge.x)
                bestNormal = Vector2D(edge.y, -edge.x) / std::sqrt(edgeLenSq);
                foundEdge = true;
            }
        }

        // Generate mirror ghost if within interaction distance
        if (foundEdge && bestDist > 0.0 && bestDist < supportRadius) {
            Vector2D ghostPos = fluid->pos - bestNormal * (2.0 * bestDist);
            Vector2D ghostVel;

            if (isNoSlip_) {
                // v_ghost = 2 * v_wall - v_fluid
                ghostVel = hostVel * 2.0 - fluid->vel;
            } else {
                // Free-slip reflection along outward edge normal
                Vector2D vRel = fluid->vel - hostVel;
                double vRelNorm = vRel.dot(bestNormal);
                ghostVel = fluid->vel - bestNormal * (2.0 * vRelNorm);
            }

            // Instantiate dynamic ghost as a FluidParticle
            ghosts.push_back(std::make_unique<FluidParticle>(
                ghostID--, ghostPos, ghostVel, fluid->mass, fluid->density, hostEnergy
            ));
        }
    }

    return ghosts;
}