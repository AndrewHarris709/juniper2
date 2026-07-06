#include <iostream>
#include <chrono>
#include <omp.h>

#include "SimData.h"

int main()
{
    int num_devices = omp_get_num_devices();
    std::cout << "OpenMP target devices available: " << num_devices << std::endl;
    if (num_devices == 0) {
        std::cerr << "No GPU device found — offloading will fall back to host." << std::endl;
    }

    std::cout << "Loading Data" << std::endl;

    SimData simData = SimData(".//test.csv");

    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "Beginning Density Iteration" << std::endl;
    simData.densityIterate(simData);
    std::cout << "Done Density Iteration" << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    std::cout << "        Runtime   : " << (duration / 1e6) << " seconds." << std::endl;

    return 0;
}
