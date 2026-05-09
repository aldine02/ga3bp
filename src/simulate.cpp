#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "physics.hpp"

int main(int argc, char* argv[]) {
    
    double X1, X2, X3, X4;
    double X5, X6, X7, X8;
    double X9, X10, X11, X12;

    X1 = atof(argv[1]); X2 = atof(argv[2]);
    X3 = atof(argv[3]); X4 = atof(argv[4]);
    X5 = atof(argv[5]); X6 = atof(argv[6]);
    X7 = atof(argv[7]); X8 = atof(argv[8]);
    X9 = atof(argv[9]); X10 = atof(argv[10]);
    X11 = atof(argv[11]); X12 = atof(argv[12]);

    State s;
    s.b1.x = X1; s.b1.y = X2;
    s.b1.vx = X3; s.b1.vy = X4;

    s.b2.x = X5; s.b2.y = X6;
    s.b2.vx = X7; s.b2.vy = X8;

    s.b3.x = X9; s.b3.y = X10;
    s.b3.vx = X11; s.b3.vy = X12;

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