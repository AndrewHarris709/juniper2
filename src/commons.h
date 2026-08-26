#ifndef COMMONS_H
#define COMMONS_H
#include <cmath>

namespace junipermath {
    struct Point3f {
        float x, y, z;
    };

    struct Box3f {
        float x1, y1, z1, x2, y2, z2;
    };

    inline float distBetween(Box3f bounds, float x1, float x2, float y1, float y2, float z1, float z2) {
        // Periodic Boundary-aware distance calculation
        float dx = std::abs(x1 - x2);
        float dy = std::abs(y1 - y2);
        float dz = std::abs(z1 - z2);

        if (dx > (bounds.x2 - bounds.x1) / 2) {
            dx = (bounds.x2 - bounds.x1) - dx;
        }
        if (dy > (bounds.y2 - bounds.y1) / 2) {
            dy = (bounds.y2 - bounds.y1) - dy;
        }
        if (dz > (bounds.z2 - bounds.z1) / 2) {
            dz = (bounds.z2 - bounds.z1) - dz;
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
