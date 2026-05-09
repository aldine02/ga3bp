# Build

## Compilation
If you want to build the simulation, simply run the makefile

```bash
make
```

By default, it will make the unconstrained GA and the executable to simulate the orbit. Other constraints can be run with

```bash
make barycentric
```
```bash
make collinear
```
```bash
make rotational
```

It will output `evolve` and `simulate`.

## Running GA
Run the GA using the command

```bash
./evolve
```

Output should look like this
```
 Running Genetic Algorithm...
 ----------------------------
 Generation = 0 | X1 = -0.43577 | X2 = -1.32030 | X3 = 0.28243 | X4 = 0.26412 | X5 = 0.67173 | X6 = 1.83618 | X7 = 0.16512 | X8 = 0.37859 | F(x) = -0.53467
 Generation = 10 | X1 = 0.15244 | X2 = -0.25290 | X3 = 0.31591 | X4 = -0.18315 | X5 = 0.63999 | X6 = 1.86938 | X7 = 0.15734 | X8 = 0.44500 | F(x) = -0.29465

 ...
 ...
 ...

 Generation = 4990 | X1 = 0.38926 | X2 = -1.69921 | X3 = 0.21001 | X4 = 0.15920 | X5 = -0.29245 | X6 = 0.19162 | X7 = -0.08397 | X8 = -0.74597 | F(x) = -0.00732
 Generation = 5000 | X1 = 0.38926 | X2 = -1.69921 | X3 = 0.21001 | X4 = 0.15920 | X5 = -0.29245 | X6 = 0.19162 | X7 = -0.08397 | X8 = -0.74597 | F(x) = -0.00732
```
`X` is the value that the GA generate and it represents what the main code programmed it to. For example, the collinear code

```c++
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
```

means that `X1` is the $x$ position, `X2` is the $x$ velocity, and `X3` is the $y$ velocity of the first body.


Depending on the constraint, the number of `X` will vary. The `F(x)` is the error value, the closer to 0 the better.

To simulate the orbit given the GA output, use the `simulate`

```bash
./simulate X1 X2 X3
```

manually input the `X`. I know, not really reliable but it gives us the option to use it with any number of choices too for testing.

The output will be printed to stdout by default, so best to save it to a csv

```bash
./simulate X1 X2 X3 X4 > output.csv
```

Or visualize it by piping it to the python script
```bash
./simulate X1 X2 X3 X4 | python3 scripts/visualize.py
```

More on visualization [here](../scripts/README.md)

Finally, clean the executables
```
make clean
``` 


## Customizing Constraint

Inside the [`main.cpp`](main.cpp) file, rewrite the state
```c++
template <typename T>
class ThreeBody {
public:
    static std::vector<T> Objective(const std::vector<T>& x) {
        State s;

        /* Rewrite here
        s.b1.x = x[0]; s.b1.y = x[1];
        s.b1.vx = x[2]; s.b1.vy = x[3];

        s.b2.x = x[4]; s.b2.y = x[5];
        s.b2.vx = x[6]; s.b2.vy = x[7];

        s.b3.x = x[8]; s.b3.y = x[9];
        s.b3.vx = x[10]; s.b3.vy = x[11];
        */

        State initial = s;

        verlet(s, 10, 0.001); 
        double error = fitness(initial, s);
        return {error};
}
};

```

We can also vary the simulation time per initial state evaluation, `verlet(s, 10, 0.001)` means $T=10$ and $\Delta{t} = 0.001$.

The main function
```c++
int main() {

    /* Change the parameter here
    galgo::Parameter<double> b1_x({-2.0, 2.0}), b1_y({-2.0, 2.0});
    galgo::Parameter<double> b1_vx({-1.0, 1.0}),  b1_vy({-1.0, 1.0});

    galgo::Parameter<double> b2_x({-2.0, 2.0}), b2_y({-2.0, 2.0});
    galgo::Parameter<double> b2_vx({-1.0, 1.0}),  b2_vy({-1.0, 1.0});

    galgo::Parameter<double> b3_x({-2.0, 2.0}), b3_y({-2.0, 2.0});
    galgo::Parameter<double> b3_vx({-1.0, 1.0}),  b3_vy({-1.0, 1.0});
    */

    // initiliazing genetic algorithm
    galgo::GeneticAlgorithm<double> ga(ThreeBody<double>::Objective, 500, 5000, true, 
                                       /* And also here
                                       b1_x, b1_y, b1_vx, b1_vy, 
                                       b2_x, b2_y, b2_vx, b2_vy,
                                       b2_x, b2_y, b2_vx, b2_vy */);

    ga.mutrate = 0.2;
    ga.run();

    return 0;
}
```
we should change the number of parameter of the GA following the constraint on the `Objective(const std::vector<T>& x)`.

`b1_x({-2.0, 2.0})` means the GA guesses a random number between -2.0 and 2.0 for the initial $x$ position of body 1.

`ga.mutrate = 0.2` means 20% mutation rate. For localized search, a good rate is around 0.1% to 1.0%.


The [`simulate.cpp`](simulate.cpp) would also needs a change
```c++
int main(int argc, char* argv[]) {
    
    /* Change here
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

    */

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
```
Done. 