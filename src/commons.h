#include <cmath>

#ifndef COMMONS_H
#define COMMONS_H

namespace junipermath {
    struct Point3f {
        float x, y, z;
    };

    struct Box3f {
        float x1, y1, z1, x2, y2, z2;
    };

    inline float distBetween(float x1, float x2, float y1, float y2, float z1, float z2) {
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
}

#endif //COMMONS_H
