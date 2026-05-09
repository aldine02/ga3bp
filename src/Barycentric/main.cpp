#include <iostream>
#include <vector>
#include "Galgo.hpp"
#include "physics.hpp"
#include "fitness.hpp"

template <typename T>
class ThreeBody {
public:
    static std::vector<T> Objective(const std::vector<T>& x) {
        State s;

        s.b1.x = x[0]; s.b1.y = x[1];
        s.b1.vx = x[2]; s.b1.vy = x[3];

        s.b2.x = x[4]; s.b2.y = x[5];
        s.b2.vx = x[6]; s.b2.vy = x[7];

        s.b3.x = -(s.b1.x + s.b2.x);
        s.b3.y = -(s.b1.y + s.b2.y);
        s.b3.vx = -(s.b1.vx + s.b2.vx);
        s.b3.vy = -(s.b1.vy + s.b2.vy);

        State initial = s;

        verlet(s, 10, 0.001); 
        double error = fitness(initial, s);
        return {error};
}
};

int main() {

    galgo::Parameter<double> b1_x({-2.0, 2.0}), b1_y({-2.0, 2.0});
    galgo::Parameter<double> b1_vx({-1.0, 1.0}),  b1_vy({-1.0, 1.0});

    galgo::Parameter<double> b2_x({-2.0, 2.0}), b2_y({-2.0, 2.0});
    galgo::Parameter<double> b2_vx({-1.0, 1.0}),  b2_vy({-1.0, 1.0});

    // initiliazing genetic algorithm
    galgo::GeneticAlgorithm<double> ga(ThreeBody<double>::Objective, 500, 5000, true, 
                                       b1_x, b1_y, b1_vx, b1_vy, 
                                       b2_x, b2_y, b2_vx, b2_vy);

    ga.mutrate = 0.2;
    ga.run();

    return 0;
}