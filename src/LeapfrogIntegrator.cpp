//
// Created by Andrew on 2026-07-31.
//

#include "LeapfrogIntegrator.h"

#include <cmath>

#include "GPUOperators.h"

#pragma omp begin declare target
namespace LeapfrogIntegrator {
    float densityForParticle(int partIndex, int nodeIndex, float *xyzh, TreeData *tree, SimConfigDevice *config) {
        float density = 0.0f;

        GPU::NeighbourList neighbours;
        float newH = xyzh[4 * partIndex + 3];
        getNeighbours(tree, nodeIndex, config, newH, &neighbours);

        for (int i = 0; i < neighbours.count; i++) {
            int partB = neighbours.indices[i];
            float dist = GPU::distBetween(config, xyzh, partIndex, partB);
            float tarH = xyzh[4 * partIndex + 3];
            density += GPU::valueAt(dist / tarH) * config->mass / (newH * newH * newH);
        }
        return density;
    }

    float omegaForParticle(int partIndex, int nodeIndex, float density, float *xyzh, TreeData *tree, SimConfigDevice *config) {
        float omega = 0.0f;
        float newH = xyzh[4 * partIndex + 3];
        float grad = -1 * (newH / (3 * density));

        GPU::NeighbourList neighbours;
        getNeighbours(tree, nodeIndex, config, newH, &neighbours);

        for (int i = 0; i < neighbours.count; i++) {
            int partB = neighbours.indices[i];
            float dist = GPU::distBetween(config, xyzh, partIndex, partB);
            omega += config->mass * GPU::dWdhAt(dist / newH);
        }

        return 1 - grad * omega / (newH * newH * newH * newH);
    }

    float pressureForParticle(int partIndex, float density, float *vxyzu, TreeData *tree, SimConfigDevice *simConfig) {
        return (GAMMA - 1) * density * vxyzu[4 * partIndex + 3];
    }

    junipermath::Point3f accForParticle(int partIndex, int nodeIndex, float *xyzh, float *vxyzu, float *density, float* omega, TreeData *tree, SimConfigDevice *config) {
        float ax = 0, ay = 0, az = 0;
        float partDensity = density[partIndex];
        float partOmega = omega[partIndex];
        float partPressure = pressureForParticle(partIndex, partDensity, vxyzu, tree, config);

        GPU::NeighbourList neighbours;
        getNeighbours(tree, nodeIndex, config, xyzh[4 * partIndex + 3], &neighbours);
        for (int i = 0; i < neighbours.count; i++) {
            int partNeighbour = neighbours.indices[i];

            float partQAB = qabForParticle(partIndex, partNeighbour, xyzh, vxyzu, partDensity, tree, config);
            float partRatio = (partPressure + partQAB) / (partDensity * partDensity * partOmega);

            float neighbourDensity = density[partNeighbour];
            float neighbourOmega = omega[partNeighbour];
            float neighbourQAB = qabForParticle(partNeighbour, partIndex, xyzh, vxyzu, neighbourDensity, tree, config);
            float neighbourPressure = pressureForParticle(partNeighbour, neighbourDensity, vxyzu, tree, config);
            float neighbourRatio = (neighbourPressure + neighbourQAB) / (neighbourDensity * neighbourDensity * neighbourOmega);

            junipermath::Point3f disp = displacementBetween(partIndex, partNeighbour, xyzh, config);
            junipermath::Point3f dispNorm = norm(disp);
            float distance = GPU::distBetween(config, xyzh, partIndex, partNeighbour);
            float hPart = xyzh[4 * partIndex + 3], hNeighbour = xyzh[4 * partNeighbour + 3];

            float partGradient = GPU::gradientAt(distance / hPart) / (hPart * hPart * hPart * hPart);
            float neighbourGradient = GPU::gradientAt(distance / hNeighbour) / (hNeighbour * hNeighbour * hNeighbour * hNeighbour);

            ax -= config->mass * (partRatio * dispNorm.x * partGradient + neighbourRatio * dispNorm.x * neighbourGradient);
            ay -= config->mass * (partRatio * dispNorm.y * partGradient + neighbourRatio * dispNorm.y * neighbourGradient);
            az -= config->mass * (partRatio * dispNorm.z * partGradient + neighbourRatio * dispNorm.z * neighbourGradient);
        }

        return junipermath::Point3f{ax, ay, az};
    }

    float energyChangeForParticle(int partIndex, int nodeIndex, float *xyzh, float *vxyzu, float *density, float* omega, TreeData *tree, SimConfigDevice *config) {
        float uChange = 0;

        float partDensity = density[partIndex];
        float partOmega = omega[partIndex];
        float partPressure = pressureForParticle(partIndex, partDensity, vxyzu, tree, config);
        float partRatio = partPressure / (partDensity * partDensity * partOmega);
        float hPart = xyzh[4 * partIndex + 3];

        GPU::NeighbourList neighbours;
        getNeighbours(tree, nodeIndex, config, hPart, &neighbours);
        for (int i = 0; i < neighbours.count; i++) {
            int partNeighbour = neighbours.indices[i];

            junipermath::Point3f vDiff = velocityDiffBetween(partIndex, partNeighbour, vxyzu);
            junipermath::Point3f disp = displacementBetween(partIndex, partNeighbour, xyzh, config);
            junipermath::Point3f dispNorm = norm(disp);

            float distance = GPU::distBetween(config, xyzh, partIndex, partNeighbour);
            float partGradient = GPU::gradientAt(distance / hPart) / (hPart * hPart * hPart * hPart);

            float xComp = vDiff.x * dispNorm.x * partGradient;
            float yComp = vDiff.y * dispNorm.y * partGradient;
            float zComp = vDiff.z * dispNorm.z * partGradient;

            uChange += config->mass * (xComp + yComp + zComp);
        }

        return uChange * partRatio;
    }

    float qabForParticle(int partIndex, int otherPart, float *xyzh, float *vxyzu, float densityA, TreeData *tree, SimConfigDevice *config) {
        const junipermath::Point3f velocityDiff = velocityDiffBetween(partIndex, otherPart, vxyzu);
        const junipermath::Point3f displacement = displacementBetween(partIndex, otherPart, xyzh, config);
        const junipermath::Point3f dispNorm = norm(displacement);

        const float losDot = dot(velocityDiff, dispNorm);
        if (losDot >= 0) {
            return 0;
        }

        const float pressureA = pressureForParticle(partIndex, densityA, vxyzu, tree, config);
        const float soundSpeed = sqrt(GAMMA * pressureA / densityA);
        const float signalSpeed = 1 * soundSpeed + 2 * std::abs(losDot);

        return -0.5 * densityA * signalSpeed * losDot;
    }

    junipermath::Point3f norm(junipermath::Point3f p) {
        const float mag = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);

        if (mag <= 0) {
            return junipermath::Point3f(0, 0, 0);
        }
        return {p.x / mag, p.y / mag, p.z / mag};
    }

    float dot(const junipermath::Point3f p1, const junipermath::Point3f p2) {
        return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
    }

    junipermath::Point3f velocityDiffBetween(int target, int part, float* vxyzu) {
        return {vxyzu[4 * target + 0] - vxyzu[4 * part + 0],
                    vxyzu[4 * target + 1] - vxyzu[4 * part + 1],
                    vxyzu[4 * target + 2] - vxyzu[4 * part + 2]};
    }

    junipermath::Point3f displacementBetween(int target, int part, float* xyzh, SimConfigDevice *config) {
        return displacementBetween(xyzh[4 * target + 0], xyzh[4 * part + 0],
                    xyzh[4 * target + 1], xyzh[4 * part + 1],
                    xyzh[4 * target + 2], xyzh[4 * part + 2], config);
    }

    junipermath::Point3f displacementBetween(float x1, float x2, float y1, float y2, float z1, float z2, SimConfigDevice* config) {
        float dx = x1 - x2;
        float dy = y1 - y2;
        float dz = z1 - z2;

        float lx = config->boundingBox.x2 - config->boundingBox.x1;
        float ly = config->boundingBox.y2 - config->boundingBox.y1;
        float lz = config->boundingBox.z2 - config->boundingBox.z1;

        if (dx > lx / 2) {
            dx -= lx;
        } else if (dx < -lx / 2) {
            dx += lx;
        }

        if (dy > ly / 2) {
            dy -= ly;
        } else if (dy < -ly / 2) {
            dy += ly;
        }

        if (dz > lz / 2) {
            dz -= lz;
        } else if (dz < -lz / 2) {
            dz += lz;
        }

        return {dx, dy, dz};
    }
}

#pragma omp end declare target
