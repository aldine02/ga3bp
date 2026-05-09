#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "physics.hpp"

int main(int argc, char* argv[]) {
    
    double X1, X2, X3;

    X1 = atof(argv[1]); 
    X2 = atof(argv[2]);
    X3 = atof(argv[3]);
    
    State s;
    s.b1.x = X1; s.b1.y = 0.0;
    s.b1.vx = X2; s.b1.vy = X3;

    s.b2.x = -X1; s.b2.y = 0.0;
    s.b2.vx = X2; s.b2.vy = X3;

    s.b3.x = 0.0;      
    s.b3.y = 0.0;

    // v3 = -(v1 + v2)
    s.b3.vx = -2.0 * X2; 
    s.b3.vy = -2.0 * X3;

    auto logger = [](const State& s, double t) {
    std::cout << std::fixed << std::setprecision(6)
       << t << "," << s.b1.x << "," << s.b1.y << "," 
       << s.b2.x << "," << s.b2.y << "," 
       << s.b3.x << "," << s.b3.y << "\n";
    };

    verlet(s, 10.0, 0.001, logger);

    return 0;
}