#include <filesystem>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <span>
#include "SimData.h"

std::vector<float> toVector(float* arr, int size) {
    return std::vector<float>(arr, arr + size);
}

void assertTestSim(SimData& simData) {
    ASSERT_THAT(toVector(simData.xyzh, simData.getParticleCount() * 4), testing::ElementsAre(3, 5, 2, 0.5, 1, 5, 3, 2, 2, 20, 1, 0.7, 5, 10, 2, 0.1));
    ASSERT_THAT(toVector(simData.vxyzu, simData.getParticleCount() * 4), testing::ElementsAre(1.5, 2.5, 0.25, 0.3, 4, 1, -0.5, 0.3, 0.5, 0.5, 0.5, 0.3, 0.5, 0.5, 0.5, 0.3));
    ASSERT_THAT(toVector(simData.fxyz, simData.getParticleCount() * 3), testing::ElementsAre(1, 2, 3, 3, 2.3, 7.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5));
}

TEST(IOTest, BasicReaderTest) {
    SimData data = SimData(".//files//test.csv", ".//files//test.toml");

    assertTestSim(data);
}

TEST(IOTest, BasicWriterTest) {
    SimData simData = SimData(".//files//test.csv", ".//files//test.toml");
    simData.toCSV(".//files//temp_writing_test.csv");
    SimData writtenSim = SimData(".//files//temp_writing_test.csv", ".//files//test.toml");

    try {
        std::filesystem::remove(".//files//temp_writing_test.csv");
    } catch (const std::filesystem::filesystem_error& err) {
        std::cout << "Error deleting temporary file!: " << err.what() << "\n";
    }

    assertTestSim(writtenSim);
}

