#include <iostream>
#include <vector>
#include <cmath>
#include "Galgo.hpp"
#include "physics.hpp"
#include "fitness.hpp"

const double PI = std::acos(-1.0);

template <typename T>
class ThreeBody {
public:
    static std::vector<T> Objective(const std::vector<T>& x) {
        State s;

        double a120 = 2.0 * PI / 3.0;
        double a240 = 4.0 * PI / 3.0;

        s.b1.x = x[0]; s.b1.y = x[1];
        s.b1.vx = x[2]; s.b1.vy = x[3];

        // rotate b1 by 120°
        s.b2.x  = s.b1.x * cos(a120) - s.b1.y * sin(a120);
        s.b2.y  = s.b1.x * sin(a120) + s.b1.y * cos(a120);
        s.b2.vx = s.b1.vx * cos(a120) - s.b1.vy * sin(a120);
        s.b2.vy = s.b1.vx * sin(a120) + s.b1.vy * cos(a120);

        // rotate b1 by 240°
        s.b3.x  = s.b1.x * cos(a240) - s.b1.y * sin(a240);
        s.b3.y  = s.b1.x * sin(a240) + s.b1.y * cos(a240);
        s.b3.vx = s.b1.vx * cos(a240) - s.b1.vy * sin(a240);
        s.b3.vy = s.b1.vx * sin(a240) + s.b1.vy * cos(a240);

        State initial = s;

        verlet(s, 10, 0.001); 
        double error = fitness(initial, s);
        return {error};
}
};

int main() {

    galgo::Parameter<double> b1_x({-2.0, 2.0}), b1_y({-2.0, 2.0});
    galgo::Parameter<double> b1_vx({-1.0, 1.0}),  b1_vy({-1.0, 1.0});

    // initiliazing genetic algorithm
    galgo::GeneticAlgorithm<double> ga(ThreeBody<double>::Objective, 500, 5000, true, 
                                       b1_x, b1_y, b1_vx, b1_vy);

    ga.mutrate = 0.2;
    ga.run();

    return 0;
}