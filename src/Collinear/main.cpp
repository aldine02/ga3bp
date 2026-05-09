#include <iostream>
#include <vector>
#include "Galgo.hpp"
#include "physics.hpp"
#include "fitness.hpp"

template <typename T>
class ThreeBody {
public:
static std::vector<double> Objective(const std::vector<double>& x) {

    State s;

    s.b1.x = x[0];  
    s.b1.y = 0.0;
    s.b1.vx = x[1];
    s.b1.vy = x[2];

    s.b2.x = -s.b1.x;
    s.b2.y = 0.0;
    s.b2.vx = s.b1.vx;
    s.b2.vy = s.b1.vy;

    s.b3.x = 0.0;      
    s.b3.y = 0.0;
    s.b3.vx = -2.0 * s.b1.vx; 
    s.b3.vy = -2.0 * s.b1.vy;

    State initial = s;

    verlet(s, 10, 0.001); 
    double error = fitness(initial, s);
    return {error};
}
};

int main() {

    galgo::Parameter<double> b1_x({-2.0, 2.0}), b1_vx({-2.0, 2.0}), b1_vy({0.0, 4.0});

    // initiliazing genetic algorithm
    galgo::GeneticAlgorithm<double> ga(ThreeBody<double>::Objective, 500, 5000, true, b1_x, b1_vx, b1_vy);

    ga.mutrate = 0.2;
    ga.run();

    return 0;
}