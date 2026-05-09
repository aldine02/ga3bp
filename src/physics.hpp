#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include <vector>
#include <cmath>
#include <functional>

struct Body {
    double x, y, vx, vy, ax, ay;
    double mass = 1.0;
};

struct State {
    Body b1, b2, b3;
};

inline void derivs(State& s) {
    
    s.b1.ax = 0; s.b1.ay = 0;
    s.b2.ax = 0; s.b2.ay = 0;
    s.b3.ax = 0; s.b3.ay = 0;

    const double epsilon = 0.0001;

    // b1 and b2
    double dx12 = s.b2.x - s.b1.x;
    double dy12 = s.b2.y - s.b1.y;
    double r2_12 = (dx12 * dx12) + (dy12 * dy12) + epsilon;
    double r_12 = std::sqrt(r2_12);

    s.b1.ax += (s.b2.mass / r2_12) * (dx12 / r_12);
    s.b1.ay += (s.b2.mass / r2_12) * (dy12 / r_12);

    s.b2.ax -= (s.b1.mass / r2_12) * (dx12 / r_12);
    s.b2.ay -= (s.b1.mass / r2_12) * (dy12 / r_12);

    // b1 and b3
    double dx13 = s.b3.x - s.b1.x;
    double dy13 = s.b3.y - s.b1.y;
    double r2_13 = (dx13 * dx13) + (dy13 * dy13) + epsilon;
    double r_13 = std::sqrt(r2_13);

    s.b1.ax += (s.b3.mass / r2_13) * (dx13 / r_13);
    s.b1.ay += (s.b3.mass / r2_13) * (dy13 / r_13);
    
    s.b3.ax -= (s.b1.mass / r2_13) * (dx13 / r_13);
    s.b3.ay -= (s.b1.mass / r2_13) * (dy13 / r_13);

    // b2 and b3
    double dx23 = s.b3.x - s.b2.x;
    double dy23 = s.b3.y - s.b2.y;
    double r2_23 = (dx23 * dx23) + (dy23 * dy23) + epsilon;
    double r_23 = std::sqrt(r2_23);

    s.b2.ax += (s.b3.mass / r2_23) * (dx23 / r_23);
    s.b2.ay += (s.b3.mass / r2_23) * (dy23 / r_23);
    
    s.b3.ax -= (s.b2.mass / r2_23) * (dx23 / r_23);
    s.b3.ay -= (s.b2.mass / r2_23) * (dy23 / r_23);
}

inline void verlet(State& s, double T, double dt, std::function<void(const State&, double)> callback = nullptr) {
    
    derivs(s);

    double dt2 = dt * dt;

    for (double t = 0; t < T; t += dt) {

        if (callback) callback(s, t);

        s.b1.x += s.b1.vx * dt + 0.5 * s.b1.ax * dt2;
        s.b1.y += s.b1.vy * dt + 0.5 * s.b1.ay * dt2;
        s.b2.x += s.b2.vx * dt + 0.5 * s.b2.ax * dt2;
        s.b2.y += s.b2.vy * dt + 0.5 * s.b2.ay * dt2;
        s.b3.x += s.b3.vx * dt + 0.5 * s.b3.ax * dt2;
        s.b3.y += s.b3.vy * dt + 0.5 * s.b3.ay * dt2;

        State old = s;

        // Derive new acceleration
        derivs(s);
        
        // Solve velocity
        s.b1.vx += 0.5 * (old.b1.ax + s.b1.ax) * dt;
        s.b1.vy += 0.5 * (old.b1.ax + s.b1.ay) * dt;
        s.b2.vx += 0.5 * (old.b1.ax + s.b2.ax) * dt;
        s.b2.vy += 0.5 * (old.b1.ax + s.b2.ay) * dt;
        s.b3.vx += 0.5 * (old.b1.ax + s.b3.ax) * dt;
        s.b3.vy += 0.5 * (old.b1.ax + s.b3.ay) * dt;
        
    }
}

#endif