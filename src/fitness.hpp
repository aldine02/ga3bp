#ifndef FITNESS_HPP
#define FITNESS_HPP
#include <cmath>
#include <algorithm>
#include "physics.hpp"

inline double euclidian(const State& initial, const State& final) {
    double d1 = std::pow(final.b1.x - initial.b1.x, 2) + std::pow(final.b1.y - initial.b1.y, 2);
    double d2 = std::pow(final.b2.x - initial.b2.x, 2) + std::pow(final.b2.y - initial.b2.y, 2);
    double d3 = std::pow(final.b3.x - initial.b3.x, 2) + std::pow(final.b3.y - initial.b3.y, 2);
    
    return std::sqrt(d1 + d2 + d3);
}

inline double fitness(const State& initial, const State& final) {
    // ejection Check
    double max_dist = std::max({
        std::sqrt(final.b1.x*final.b1.x + final.b1.y*final.b1.y),
        std::sqrt(final.b2.x*final.b2.x + final.b2.y*final.b2.y),
        std::sqrt(final.b3.x*final.b3.x + final.b3.y*final.b3.y)
    });
    if (max_dist > 4.0) return -9999.0;

    // proximity Check
    double d12 = std::sqrt(std::pow(final.b2.x - final.b1.x, 2) + std::pow(final.b2.y - final.b1.y, 2));
    double d13 = std::sqrt(std::pow(final.b3.x - final.b1.x, 2) + std::pow(final.b3.y - final.b1.y, 2));
    double d23 = std::sqrt(std::pow(final.b3.x - final.b2.x, 2) + std::pow(final.b3.y - final.b2.y, 2));
    if (d12 < 0.2 || d13 < 0.2 || d23 < 0.2) return -8888.0;

    double error = euclidian(initial, final);

    if (std::isnan(error) || std::isinf(error)) return -7777.0;

    return -error; 
}

#endif