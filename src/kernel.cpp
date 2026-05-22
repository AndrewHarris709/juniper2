#include <cmath>
#include "kernel.h"

constexpr float NORM = 1 / M_PI;

float Kernel::getRadius() const
{
    return 2;
}

float Kernel::valueAt(const float q) const
{
    float aq = std::abs(q);

    if (aq >= 2) 
    {
        return 0;
    }

    if (aq >= 1)
    {
        return NORM * 0.25 * (2 - aq) * (2 - aq) * (2 - aq);
    }

    return NORM * (1 - 1.5 * aq * aq + 0.75 * aq * aq * aq);
}

float Kernel::gradientAt(const float q) const {
    float aq = std::abs(q);

    if (aq >= 2) {
        return 0;
    }

    if (aq >= 1) {
        return NORM * -0.75 * (aq - 2) * (aq - 2);
    }

    return NORM * aq * (2.25 * aq - 3);
}

float Kernel::dWdhAt(const float q) const {
    return -1.0 * NORM * (3 * this->valueAt(q) + q * this->gradientAt(q));
}
