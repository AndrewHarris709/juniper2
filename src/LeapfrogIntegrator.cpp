//
// Created by Andrew on 2026-07-31.
//

#include "LeapfrogIntegrator.h"

#include <cassert>
#include <cmath>
#include <omp.h>

#include "GPUOperators.h"
#include "SimData.h"

namespace LeapfrogIntegrator {

#pragma omp declare target
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

    float periodicWrap(float coord, float lower, float upper) {
        const float length = upper - lower;

        if (coord >= upper) {
            return coord - length;
        }
        if (coord < lower) {
            return coord + length;
        }
        return coord;
    }
#pragma omp end declare target

    constexpr int THREAD_LIMIT = 512;

    // Builds the neighbour tree from the host copy of xyzh (CPU), then computes the acceleration and energy change of
    // every particle (GPU). Expects xyzh, vxyzu, density, omega, accs and energies to already be mapped on the device.
    static void computeAccelerationAndEnergy(SimData& data) {
        Tree tree(data.getParticleCount());
        tree.build(data);

        TreeNode* contents = tree.data->contents;
        int* mapping = tree.data->mapping;
        const int nodeCount = tree.getNodeCount();
        const int partCount = tree.data->partCount;
        float* xyzh = data.xyzh;
        float* vxyzu = data.vxyzu;
        float* density = data.density;
        float* omega = data.omega;
        float* accs = data.accs;
        float* energies = data.energies;
        SimConfigDevice* config = data.deviceConfig;

        // The acceleration of a particle depends on the density & omega of its neighbours, so these must be complete
        // for every particle (hence a separate kernel) before any acceleration is computed.
        #pragma omp target teams distribute parallel for  map(to: contents[0:nodeCount]) \
                                                          map(to: mapping[0:partCount]) \
                                                          thread_limit(THREAD_LIMIT)
        for (int nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
            if (contents[nodeIndex].leftChild != -1) {
                // non-leaf node
                continue;
            }

            TreeData td{contents, mapping, nodeCount, partCount};
            NodeRange partIndices = GPU::getPartsFromNode(&td, nodeIndex);
            for (int i = 0; i < partIndices.size; i++) {
                int partIndex = partIndices.data[i];
                density[partIndex] = densityForParticle(partIndex, nodeIndex, xyzh, &td, config);
                omega[partIndex] = omegaForParticle(partIndex, nodeIndex, density[partIndex], xyzh, &td, config);
            }
        }

        #pragma omp target teams distribute parallel for  map(to: contents[0:nodeCount]) \
                                                          map(to: mapping[0:partCount]) \
                                                          thread_limit(THREAD_LIMIT)
        for (int nodeIndex = 0; nodeIndex < nodeCount; nodeIndex++) {
            if (contents[nodeIndex].leftChild != -1) {
                // non-leaf node
                continue;
            }

            TreeData td{contents, mapping, nodeCount, partCount};
            NodeRange partIndices = GPU::getPartsFromNode(&td, nodeIndex);
            for (int i = 0; i < partIndices.size; i++) {
                int partIndex = partIndices.data[i];
                junipermath::Point3f acc = accForParticle(partIndex, nodeIndex, xyzh, vxyzu, density, omega, &td, config);
                accs[3 * partIndex + 0] = acc.x;
                accs[3 * partIndex + 1] = acc.y;
                accs[3 * partIndex + 2] = acc.z;
                energies[partIndex] = energyChangeForParticle(partIndex, nodeIndex, xyzh, vxyzu, density, omega, &td, config);
            }
        }
    }

    void step(SimData& data, const float timestep) {
        const int n = data.getParticleCount();
        const float halfStep = 0.5f * timestep;
        float* xyzh = data.xyzh;
        float* vxyzu = data.vxyzu;
        float* density = data.density;
        float* omega = data.omega;
        float* accs = data.accs;
        float* energies = data.energies;
        SimConfigDevice* config = data.deviceConfig;

        // v^{n+1/2} and u^{n+1/2}, set aside while vxyzu holds the predicted values used by the second force evaluation
        auto* vxyzuHalf = new float[n * 4];

        assert(omp_target_is_present(xyzh, omp_get_default_device()) && "xyzh not mapped on device");
        assert(omp_target_is_present(config, omp_get_default_device()) && "config not mapped on device");

        #pragma omp target data map(tofrom: vxyzu[0:n * 4]) \
                                map(from: density[0:n], accs[0:n * 3], energies[0:n]) \
                                map(alloc: omega[0:n], vxyzuHalf[0:n * 4])
        {
            // a^n and du/dt^n from x^n, v^n, u^n
            computeAccelerationAndEnergy(data);

            // Kick by a half step, then drift by a full step:
            //     v^{n+1/2} = v^n + dt/2 * a^n,    x^{n+1} = x^n + dt * v^{n+1/2}
            // The forces depend on velocity (artificial viscosity, du/dt), so they are next evaluated at the predicted
            // v* = v^{n+1/2} + dt/2 * a^n rather than v^{n+1/2}, as in Phantom (Price et al. 2018). Likewise for u.
            #pragma omp target teams distribute parallel for
            for (int i = 0; i < n; i++) {
                for (int dim = 0; dim < 3; dim++) {
                    float vHalf = vxyzu[4 * i + dim] + halfStep * accs[3 * i + dim];
                    vxyzuHalf[4 * i + dim] = vHalf;
                    vxyzu[4 * i + dim] = vHalf + halfStep * accs[3 * i + dim];
                    xyzh[4 * i + dim] += timestep * vHalf;
                }

                float uHalf = vxyzu[4 * i + 3] + halfStep * energies[i];
                vxyzuHalf[4 * i + 3] = uHalf;
                vxyzu[4 * i + 3] = uHalf + halfStep * energies[i];

                xyzh[4 * i + 0] = periodicWrap(xyzh[4 * i + 0], config->boundingBox.x1, config->boundingBox.x2);
                xyzh[4 * i + 1] = periodicWrap(xyzh[4 * i + 1], config->boundingBox.y1, config->boundingBox.y2);
                xyzh[4 * i + 2] = periodicWrap(xyzh[4 * i + 2], config->boundingBox.z1, config->boundingBox.z2);
            }

            // The tree is built on the host, so it needs the drifted positions
            data.pullGPU();

            // a^{n+1} and du/dt^{n+1} from x^{n+1}, v*, u*
            computeAccelerationAndEnergy(data);

            // Kick by the remaining half step: v^{n+1} = v^{n+1/2} + dt/2 * a^{n+1}
            #pragma omp target teams distribute parallel for
            for (int i = 0; i < n; i++) {
                for (int dim = 0; dim < 3; dim++) {
                    vxyzu[4 * i + dim] = vxyzuHalf[4 * i + dim] + halfStep * accs[3 * i + dim];
                }
                vxyzu[4 * i + 3] = vxyzuHalf[4 * i + 3] + halfStep * energies[i];
            }
        }

        delete[] vxyzuHalf;
        data.time += timestep;
    }
}
