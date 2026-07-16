#include <random>
#include <gtest/gtest.h>
#include <stack>
#include "SimData.h"

TEST(DensityTest, PhantomDensityTest) {
    SimData data = SimData(".//files//hydro32.csv", ".//files//hydro32.toml");

    // Even if we randomly perturb the h-values of some particles, we should still return to the same state.
    std::vector<float> oldxyzh(data.xyzh, data.xyzh + 4 * data.getParticleCount());
    std::default_random_engine el(15);
    std::uniform_real_distribution<float> distribution(0.9, 1.1);
    for (int i = 0; i < data.getParticleCount(); i++) {
        data.h(i) = data.h(i) * distribution(el);
    }

    data.densityIterate();
    for (int i = 0; i < data.getParticleCount(); i++) {
        EXPECT_NEAR(data.h(i), oldxyzh[4 * i + 3], 1e-4) << "Mismatch at particle index i = " << i;
    }
}
