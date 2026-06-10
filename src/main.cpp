#include <iostream>
#include <chrono>
#include <cmath>
#include <thread>
#include <omp.h>
#include <mutex>
#include "kernel.h"
#include "SimData.h"

#pragma omp declare target
float distBetween(float x1, float x2, float y1, float y2, float z1, float z2) {
    // Periodic Boundary-aware distance calculation
    float dx = std::abs(x1 - x2);
    float dy = std::abs(y1 - y2);
    float dz = std::abs(z1 - z2);

    float temp_bounds[6]{0, 1, 0, 1, 0, 1};

    if (dx > (temp_bounds[1] - temp_bounds[0]) / 2) {
        dx = (temp_bounds[1] - temp_bounds[0]) - dx;
    }
    if (dy > (temp_bounds[3] - temp_bounds[2]) / 2) {
        dy = (temp_bounds[3] - temp_bounds[2]) - dy;
    }
    if (dz > (temp_bounds[5] - temp_bounds[4]) / 2) {
        dz = (temp_bounds[5] - temp_bounds[4]) - dz;
    }

    if (x1 > x2)
        dx *= -1;
    if (y1 > y2)
        dy *= -1;
    if (z1 > z2)
        dz *= -1;

    return sqrtf(dx * dx + dy * dy + dz * dz);
}
#pragma omp end declare target

float densityAt(SimData& simData, int part) {
    float density = 0.0;
    Kernel kernel = Kernel();
    int pCount = simData.getParticleCount();

    // You cannot offload parts of an unsafe class, even if those parameters are safe by themselves.
    float  m     = simData.m;
    float* xyzh  = simData.xyzh;
    int    count = pCount;

    // On OpenMP 5.0+: #pragma omp target teams distribute parallel for reduction(+:density) map(to: m) map(to:pCount) map(to:kernel) map(to:part) map(tofrom:density) map(present: xyzh[0:pcount*4]) defaultmap(none)

    #pragma omp target enter data map(to: xyzh[0:pCount*4])
    #pragma omp target teams distribute parallel for reduction(+:density) \
                                                     map(to: m, pCount, kernel, part) \
                                                     map(tofrom:density) \
                                                     map(present: xyzh[0:pcount*4]) \
                                                     defaultmap(firstprivate)
    for (int i = 0; i < pCount; i++) {
        float dist = distBetween(
            xyzh[4 * i], xyzh[4 * part],
            xyzh[4 * i + 1], xyzh[4 * part + 1],
            xyzh[4 * i + 2], xyzh[4 * part + 2]);

        float tarH = xyzh[4 * part+3];
        density += kernel.valueAt(dist / tarH) * m;
    }
    #pragma omp target exit data map(from: xyzh[0:pCount*4])

    return density;
}

int main()
{
    int num_devices = omp_get_num_devices();
    std::cout << "OpenMP target devices available: " << num_devices << std::endl;
    if (num_devices == 0) {
        std::cerr << "No GPU device found — offloading will fall back to host." << std::endl;
    }

    std::cout << "Loading Data" << std::endl;

    SimData simData = SimData(".//test.csv");

    std::cout << "Data Loaded, starting density calculation." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    float sum = densityAt(simData, 35);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    float expected = 5.27347e-05;
    std::cout << "Density Calculated: " << sum << std::endl;
    std::cout << "Expected Value    : " << expected << std::endl;
    std::cout << "  Percent Error   : " << 100.0 * std::abs(sum - expected) / expected << std::endl;
    std::cout << "        Runtime   : " << (duration / 1e6) << " seconds." << std::endl;

    return 0;
}
