# Genetic Algorithm for Periodic Solutions in the Three-Body Problem

## Intro
I was reading Liu Cixin's The Three-Body Problem when I came across the chapter where the character Wei Cheng talks about using the nature of random numbers to solve the three-body problem. 

>*it shows how, mathematically, random brute force can overcome precise logic ... This is my strategy for solving the three-body problem ... I treat each combination like a life form ... The computation proceeds by eliminating the disadvantaged and preserving the advantaged.*
>
>"It's an evolutionary algorithm," Wang Miao said.

Although the implication is different, it made me think about the application of the Genetic Algorithm (GA) to search for periodic orbits of a three-star system.

I used GA for my thesis during my undergraduate studies in mechanical engineering. Although I used it to optimize the configuration for a low-head hydroturbine to maximize efficiency, the same approach could be used.


## The Approach
I will be using C++ due to its vectorization and parallel ability with C speed, since the computation would take millions of iteration steps. Python will be used to plot using Matplotlib. Build instructions are available [here](src/README.md).


To simplify the search space, I will constrain the system to a 2D plane and assume equal mass for all bodies. The model is structured by: 

```c++
struct Body {
    double x, y, vx, vy, ax, ay;
    double mass = 1.0;
};

struct State {
    Body b1, b2, b3;
};
```

### Gravity
We will use Newton's gravity equation to simulate the physics:

$$F = ma$$

$$\vec{F}_{ij} = G \frac{m_i m_j}{\|\vec{r}_j - \vec{r}_i\|^2} \hat{r}_{ij}$$


For a three-body system, and G is normalized for simplicity, the acceleration can be written as:

$$\vec{a}_i = \sum_{j \neq i} \frac{m_j (\vec{r}_j - \vec{r}_i)}{\|\vec{r}_j - \vec{r}_i\|^3}$$


Then, we decompose the vectors into scalars and add a softening factor ($\epsilon$) to avoid infinite acceleration at near zero distance encounters.

$$r_{ij}^2 = (x_j - x_i)^2 + (y_j - y_i)^2 + \epsilon$$

$$r_{ij} = \sqrt{r_{ij}^2}$$

Therefore, acceleration of body $i$ caused by body $j$ is:

$$\vec{a_i} = \vec{a_i} + \left( \frac{m_j}{r_{ij}^2} \right) \left( \frac{\vec{r_j} - \vec{x_i}}{r_{ij}} \right)$$

Conversely, the acceleration of body $j$ caused by body $i$ is the opposite due to Newton's third law of motion:

$$F_{ij} = -F_{ji}$$

$$\vec{a_j} = \vec{a_j} - \left( \frac{m_j}{r_{ij}^2} \right) \left( \frac{\vec{r_j} - \vec{x_i}}{r_{ij}} \right)$$


### Numerical Method
When it comes to simulating orbits, we need a numerical solver that can preserve the Hamiltonian for a long period of time. The Velocity Verlet is a perfect choice for this as it is a symplectic integrator, which means it preserves energy for long-term stability.

$$\vec{r}(t + \Delta t) = \vec{r}(t) + \vec{v}(t)\Delta t + \frac{1}{2}\vec{a}(t)\Delta t^2$$

$$\vec{v}(t + \frac{1}{2}\Delta t) = \vec{v}(t) + \frac{1}{2}\vec{a}(t)\Delta t$$

$$\vec{a}(t + \Delta t) = \vec{A}(\vec{r}_1, \vec{r}_2, \vec{r}_3)$$

$$\vec{v}(t + \Delta t) = \vec{v}(t + \frac{1}{2}\Delta t) + \frac{1}{2}\vec{a}(t + \Delta t)\Delta t$$

Since the fourth step to calculate the velocity requires the half-step of the velocity (which is the second step of the integration), we simply merge them, so we get:

$$\vec{r} = \vec{r} + \vec{v}\Delta t + \frac{1}{2}\vec{a_{old}}\Delta t^2$$

$$\vec{a_{new}} = \vec{A}(\vec{r}_1, \vec{r}_2, \vec{r}_3)$$

$$\vec{v} = \vec{v} + \frac{1}{2}(\vec{a_{old}} + \vec{a_{new}}) \Delta t$$


### Genetic Algorithm
I will be using a pre-made library [GALGO 2.0](https://github.com/olmallet81/GALGO-2.0) by *Olivier Mallet* to handle the GA. It is a powerful GA library with available parallelization using OpenMP.

For simplicity, the default GA configuration from the library would be used except for the mutation rate, which I would vary between 0.5-20%. Population size is $500$, and the number of generations is $5000$ to give more time for the search.


### Fitness Function
The algorithm will start by generating a bunch of random initial conditions for each generation based on the population size. For each individual, the physics simulator will be used to calculate their orbit. After a certain time period ($T$), the fitness function will calculate the Euclidean distance given by:

$$d = \sqrt{\sum_{i=1}^{3} (x_T - x_0)^2 + (y_T - y_0)^2}\$$

to calculate how much it deviates from the origin. This way, the GA will look for initial conditions of a three-body system where each body returns to its initial position after the given time period. However, this alone has a problem, a body can get ejected out of orbit and be at the position of $\vec{r_0} = \vec{r_T}$, and the GA will think it is a stable orbit. To avoid that, we give a maximum value ($R_{max}$) of the euclidian distance allowed:

$$\max(|\vec{r}_1(T)|, |\vec{r}2(T)|, |\vec{r}3(T)|) ≯  R_{max} $$

Finally, we have to make sure that the GA doesn't just make the bodies orbit each other at an extremely close distance. We give a minimum value ($R_{min}$) of how close the permitted distance between each bodies: 

$$d_{12}, d_{13}, d_{23} ≮  R_{min} $$

The error we return should be negative since the GA library we use maximizes by default.

## Results 
### Unconstrained Search
First, let's run the GA without any constraint, meaning the GA has 12 degrees of freedom for all three bodies' starting position and velocity. 

```c++
s.b1.x = x[0]; s.b1.y = x[1];
s.b1.vx = x[2]; s.b1.vy = x[3];

s.b2.x = x[4]; s.b2.y = x[5];
s.b2.vx = x[6]; s.b2.vy = x[7];

s.b3.x = x[8]; s.b3.y = x[9];
s.b3.vx = x[10]; s.b3.vy = x[11];
```

Output:
```
 Running Genetic Algorithm...
 ----------------------------
 Generation =    0 | X1 =   0.27219 | X2 =   0.46226 | X3 =   0.41640 | X4 =   0.22618 | X5 =  -1.13457 | X6 =  -0.20938 | X7 =  -0.44057 | X8 =  -0.64181 | X9 =  -0.31339 | X10 =  -0.67344 | X11 =  -0.19136 | X12 =   0.24514 | F(x) =     -2.28541
 
 ...

 Generation = 5000 | X1 =  -0.47306 | X2 =  -0.94975 | X3 =   0.00591 | X4 =   0.65853 | X5 =   1.72198 | X6 =  -0.23465 | X7 =   0.99457 | X8 =  -0.03362 | X9 =   1.98297 | X10 =  -0.33689 | X11 =  -0.99750 | X12 =  -0.59973 | F(x) =     -0.29280

```
It converged to -0.29280 error distance to the initial position, which is quite high.

The result is a rough, drifting hierarchical system.

![](results/unconstrained.gif)

Notice that the initial line does not merge with the final line, and is guaranteed to drift if we simulate any further. However, it is satisfying to see that without any constraint, a truly random search gives out the most efficient orbit on which most three-star systems are organized.

### Barycentric Reduction
Most three-star systems are organized in a hierarchical system due to their nature to preserve momentum and the center of mass. Let's look at it further by applying a barycentric reduction and defining a constraint for the GA to the system's center of mass:

$$\frac{m_1\vec{r_1} + m_2\vec{r_2} + m_3\vec{r_3}}{\vec{r_1}+\vec{r_2}+\vec{r_3}} = 0$$

since all bodies have the same mass,

$$\vec{r_1} + \vec{r_2} + \vec{r_3} = 0$$

$$\vec{r_3} = -(\vec{r_1} + \vec{r_2})$$

Then, we conserve the momentum:

$$m_1\vec{v_1} + m_2\vec{v_2} + m_3\vec{v_3} = 0$$

$$\vec{v_1} + \vec{v_2} + \vec{v_3} = 0$$

$$\vec{v_3} = -(\vec{v_1} + \vec{v_2})$$

By constraining $r_3$ and $v_3$, we reduce the GA degree of freedom to 8.

```c++
s.b1.x = x[0]; s.b1.y = x[1];
s.b1.vx = x[2]; s.b1.vy = x[3];

s.b2.x = x[4]; s.b2.y = x[5];
s.b2.vx = x[6]; s.b2.vy = x[7];

s.b3.x = -(s.b1.x + s.b2.x);
s.b3.y = -(s.b1.y + s.b2.y);
s.b3.vx = -(s.b1.vx + s.b2.vx);
s.b3.vy = -(s.b1.vy + s.b2.vy);
```

Output:
```
 Running Genetic Algorithm...
 ----------------------------
 Generation =    0 | X1 =   0.22928 | X2 =  -0.66313 | X3 =   0.19716 | X4 =  -0.72244 | X5 =   0.29795 | X6 =   1.35009 | X7 =  -0.71725 | X8 =   0.27788 | F(x) =     -0.53168

 ...
 
 Generation = 5000 | X1 =  -0.09531 | X2 =   1.38268 | X3 =   0.74383 | X4 =   0.22881 | X5 =  -0.04538 | X6 =  -0.62516 | X7 =   0.44741 | X8 =   0.17897 | F(x) =     -0.01358

```

The result is a smooth hierarchical star system with a converged -0.01358 error in distance.

![](results/barycentric.gif)

### Rotational Symmetry
Most three-star systems aren't symmetric. If we want to search for a symmetric orbit, we have to put a constraint that mirrors the orbit at $t=0$. We do this by constraining the GA to only look for one body's position and velocity, and the other two followed its configuration at a different angle.

Here, I have the GA only configure one body, and have the other get the same configuration but rotated along the initial position axis. One body would get rotated 120° and the other 240°:

$$\theta_2 = \frac{2\pi}{3}, \quad \theta_3 = \frac{4\pi}{3}$$

$$x_2 = x_1 \cos(\theta_2) - y_1 \sin(\theta_2), \quad y_2 = x_1 \sin(\theta_2) + y_1 \cos(\theta_2)$$

$$v_{x2} = v_{x1} \cos(\theta_2) - v_{y1} \sin(\theta_2), \quad v_{y2} = v_{x1} \sin(\theta_2) + v_{y1} \cos(\theta_2)$$

$$x_3 = x_1 \cos(\theta_3) - y_1 \sin(\theta_3), \quad y_3 = x_1 \sin(\theta_3) + y_1 \cos(\theta_3)$$

$$vx_3 = v_{x1} \cos(\theta_3) - v_{y1} \sin(\theta_3), \quad v_{y3} = v_{x1} \sin(\theta_3) + v_{y1} \cos(\theta_3)$$

This would reduce the GA degree of freedom to only 4, which are only one body's initial state.

```c++
const double PI = std::acos(-1.0);
const double a120 = 2.0 * PI / 3.0;
const double a240 = 4.0 * PI / 3.0;

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
```

Output:
```
 Running Genetic Algorithm...
 ----------------------------
 Generation =    0 | X1 =  -1.76001 | X2 =   0.39658 | X3 =  -0.13622 | X4 =   0.33812 | F(x) =     -0.03077

 ...
 
 Generation = 5000 | X1 =  -0.68315 | X2 =  -1.55273 | X3 =  -0.38961 | X4 =   0.14220 | F(x) =     -0.00011

```

The result is a beautifully converged -0.00011 error in distance.   
![](results/rotational.gif)

The GA balances the centrifugal force against gravity. The surviving orbits are perfect, periodic rings. 

### Collinear Symmetry
If we took this symmetric constraint further, we could enforce a collinear configuration at $t=0$. In this configuration, the $y$ coordinates of all bodies are set in one line.

$$y_1 = y_2 = y_3 = 0$$

Then, we mirror one body with the other, with the same velocity:

$$x_1 = -x_2$$
$$\vec{v_1} = \vec{v_2}$$

Finally, the third body should be in the center, while balancing the momentum from the other bodies:

$$x_3 = 0$$
$$\vec{v_3} = -(\vec{v_1} + \vec{v_2})$$

By doing this, we reduce the degree of freedom to 3. That is the position($x$) and velocity($x,y$) for only one body.

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

Output:
```
 Running Genetic Algorithm...
 ----------------------------
 Generation =    0 | X1 =  -0.41715 | X2 =  -0.57896 | X3 =   0.89772 | F(x) =     -1.22002
 Generation =   10 | X1 =   0.37101 | X2 =   0.36820 | X3 =   0.83992 | F(x) =     -0.14973
 ...
 Generation = 4990 | X1 =   0.29258 | X2 =  -0.63359 | X3 =   0.98799 | F(x) =     -0.00636
 Generation = 5000 | X1 =   0.29258 | X2 =  -0.63359 | X3 =   0.98799 | F(x) =     -0.00636
```

And we get the famous figure-eight orbit for the three-body problem discovered by *Cristopher Moore* in 1993!
![](results/collinear.gif)

## Final Note
I had a lot of fun making this project. The program ran smoothly, considering it executed 5000 3-body simulations, and was done within a minute on a potato laptop, thanks to the parallelization. 
These are only some examples of how GA can be used to optimize constrained systems within the three-body problem. The current method to search for periodic families of orbits uses more advanced methods, such as topological tools with high-performance computing. However, the fact that giving constraints to random searches can lead to beautiful symmetries and replicate previously known orbits is astonishing. This exploration using an evolutionary algorithm shows how random numbers can lead us to finding order in chaos.  
