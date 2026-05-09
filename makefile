CXX = g++
CXXFLAGS = -O3 -std=c++11 -I./src
GALGOFLAGS = -fopenmp -I./lib/GALGO-2.0/src

default:
	$(CXX) $(CXXFLAGS) ${GALGOFLAGS} src/main.cpp -o evolve
	$(CXX) $(CXXFLAGS) src/simulate.cpp -o simulate

barycentric:
	$(CXX) $(CXXFLAGS) ${GALGOFLAGS} src/Barycentric/main.cpp -o evolve
	$(CXX) $(CXXFLAGS) src/Barycentric/simulate.cpp -o simulate

rotational:
	$(CXX) $(CXXFLAGS) ${GALGOFLAGS} src/Rotational/main.cpp -o evolve
	$(CXX) $(CXXFLAGS) src/Rotational/simulate.cpp -o simulate

collinear:
	$(CXX) $(CXXFLAGS) ${GALGOFLAGS} src/Collinear/main.cpp -o evolve
	$(CXX) $(CXXFLAGS) src/Collinear/simulate.cpp -o simulate

clean:
	rm -f evolve simulate