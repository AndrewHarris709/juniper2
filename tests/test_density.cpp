#include <random>
#include <gtest/gtest.h>
#include <stack>
#include "SimData.h"

TEST(DensityTest, PhantomDensityTest) {
    SimData data = SimData("test.csv");
    data.m = 3.0517578125e-05;

    // Even if we randomly perturb the h-values of some particles, we should still return to the same state.
    std::vector<float> oldxyzh(data.xyzh, data.xyzh + 4 * data.getParticleCount());
    std::default_random_engine el(15);
    std::uniform_real_distribution<float> distribution(0.9, 1.1);
    for (int i = 0; i < data.getParticleCount(); i++) {
        data.h(i) = data.h(i) * distribution(el);
    }

    data.densityIterate(data);
    for (int i = 0; i < 100; i++) {
        EXPECT_NEAR(data.h(i), oldxyzh[4 * i + 3], 1e-4) << "Mismatch at particle index i = " << i;
    }
}
