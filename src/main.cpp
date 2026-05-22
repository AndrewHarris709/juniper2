#include <iostream>
#include <chrono>
#include <cmath>
#include <thread>
#include <omp.h>
#include <mutex>

const int NUM_POINTS = 10000000000;
const int NUM_THREADS = 24;

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

int sample_random_points(int npoints)
{
    int points_inside = 0;

    #pragma omp target teams distribute parallel for reduction(+:points_inside) map(to:npoints) map(tofrom:points_inside) defaultmap(none)
    for (int i = 0; i < npoints; i++)
    {
        unsigned int state = (unsigned int)(i+1) * 100;

        rand(state);
        double x = rand(state);
        double y = rand(state);

        if (is_inside_circle(x, y)) {
            points_inside++;
        }
    }

    return points_inside;
}

int main()
{
    int num_devices = omp_get_num_devices();
    std::cout << "OpenMP target devices available: " << num_devices << std::endl;
    if (num_devices == 0) {
        std::cerr << "No GPU device found — offloading will fall back to host." << std::endl;
    }

    auto start = std::chrono::high_resolution_clock::now();

    int sum = sample_random_points(NUM_POINTS);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    double pi_estimate = 4.0 * static_cast<double>(sum) / NUM_POINTS;
    std::cout << "Estimate for pi: " << pi_estimate << std::endl;
    std::cout << "  Percent Error: " << 100.0 * std::abs(pi_estimate - M_PI) / M_PI << std::endl;
    std::cout << "        Runtime: " << (duration / 1e6) << " seconds." << std::endl;

    return 0;
}
