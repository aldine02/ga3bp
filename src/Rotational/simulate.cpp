#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include "physics.hpp"

int main(int argc, char* argv[]) {
    
    const double PI = std::acos(-1.0);
    double a120 = 2.0 * PI / 3.0;
    double a240 = 4.0 * PI / 3.0;

    double X1, X2, X3, X4;

    X1 = atof(argv[1]); X2 = atof(argv[2]);
    X3 = atof(argv[3]); X4 = atof(argv[4]);

    State s;

    s.b1.x = X1; s.b1.y = X2;
    s.b1.vx = X3; s.b1.vy = X4;

    // Star 2 (Rotate Star 1 by 120 degrees)
    s.b2.x  = s.b1.x * cos(a120) - s.b1.y * sin(a120);
    s.b2.y  = s.b1.x * sin(a120) + s.b1.y * cos(a120);
    s.b2.vx = s.b1.vx * cos(a120) - s.b1.vy * sin(a120);
    s.b2.vy = s.b1.vx * sin(a120) + s.b1.vy * cos(a120);

    // Star 3 (Rotate Star 1 by 240 degrees)
    s.b3.x  = s.b1.x * cos(a240) - s.b1.y * sin(a240);
    s.b3.y  = s.b1.x * sin(a240) + s.b1.y * cos(a240);
    s.b3.vx = s.b1.vx * cos(a240) - s.b1.vy * sin(a240);
    s.b3.vy = s.b1.vx * sin(a240) + s.b1.vy * cos(a240);


    auto logger = [](const State& s, double t) {
    std::cout << std::fixed << std::setprecision(6)
       << t << "," 
       << s.b1.x << "," << s.b1.y << "," 
       << s.b2.x << "," << s.b2.y << "," 
       << s.b3.x << "," << s.b3.y << "\n";
    };

    verlet(s, 10.0, 0.001, logger);

    return 0;
}