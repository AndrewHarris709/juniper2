#include <iostream>
#include <chrono>
#include <cmath>
#include <thread>
#include <omp.h>
#include <mutex>
#include "kernel.h"

const int NUM_POINTS = 1000000;

// We need a separate random number generator that can operate on the GPU.
#pragma omp declare target
inline double rand(unsigned int &state)
{
    // Basic Xorshift32 algorithm
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return static_cast<double>(state) / static_cast<double>(UINT_MAX);
}
#pragma omp end declare target

#pragma omp declare target
bool is_inside_circle(double x, double y)
{
    return x * x + y * y <= 1.0;
}
#pragma omp end declare target

float sample_kernel(int npoints)
{
    float total_area;

    Kernel kernel = Kernel();

    #pragma omp target teams distribute parallel for reduction(+:total_area) map(to:npoints) map(to:kernel) map(tofrom:total_area) defaultmap(none)
    for (int i = 0; i < npoints; i++)
    {
        total_area += kernel.valueAt(0);
    }

    std::cout << "Sampled Value: " << total_area << std::endl;

    return total_area / static_cast<float>(npoints);
}

int main()
{
    int num_devices = omp_get_num_devices();
    std::cout << "OpenMP target devices available: " << num_devices << std::endl;
    if (num_devices == 0) {
        std::cerr << "No GPU device found — offloading will fall back to host." << std::endl;
    }

    auto start = std::chrono::high_resolution_clock::now();

    float sum = sample_kernel(NUM_POINTS);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    std::cout << "Expected Value : " << M_1_PI * NUM_POINTS << std::endl;
    std::cout << "  Percent Error: " << 100.0 * std::abs(sum - (M_1_PI * NUM_POINTS)) / (M_1_PI * NUM_POINTS) << std::endl;
    std::cout << "        Runtime: " << (duration / 1e6) << " seconds." << std::endl;

    return 0;
}
