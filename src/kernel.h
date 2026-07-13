#pragma once

#pragma omp declare target
class Kernel
{
    public:
        float valueAt(const float) const;
        float gradientAt(const float) const;
        float dWdhAt(const float) const;
        float getRadius() const;
};
#pragma omp end declare target
