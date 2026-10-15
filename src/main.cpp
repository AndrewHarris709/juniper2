#include <iostream>
#include <chrono>
#include <filesystem>
#include <format>
#include <omp.h>

#include "LeapfrogIntegrator.h"
#include "SimData.h"

constexpr float TIMESTEP = 1e-3;
constexpr int STEP_COUNT = 500;

int main()
{
    int num_devices = omp_get_num_devices();
    std::cout << "OpenMP target devices available: " << num_devices << std::endl;
    if (num_devices == 0) {
        std::cerr << "No GPU device found — offloading will fall back to host." << std::endl;
    }

    std::cout << "Loading Data" << std::endl;

    SimData simData = SimData(".//kh_00000.csv", ".//kh_00000.toml");
    std::filesystem::create_directories(".//output");

    auto start = std::chrono::high_resolution_clock::now();

    simData.generateReports();

    std::cout << "Beginning Density Iteration" << std::endl;
    simData.densityIterate();
    std::cout << "Done Density Iteration" << std::endl;

    std::cout << "Beginning Leapfrog Steps (dt = " << TIMESTEP << ")" << std::endl;
    for (int stepIndex = 1; stepIndex <= STEP_COUNT; stepIndex++) {
        LeapfrogIntegrator::step(simData, TIMESTEP);

        std::string filename = std::format(".//output//kh_{:05d}.csv", stepIndex);
        simData.toCSV(filename);
        std::cout << "    Step " << stepIndex << "/" << STEP_COUNT << " (t = " << simData.time << "), wrote " << filename << std::endl;
    }
    std::cout << "Done Leapfrog Steps" << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    std::cout << "        Runtime   : " << (duration / 1e6) << " seconds." << std::endl;

    return 0;
}
