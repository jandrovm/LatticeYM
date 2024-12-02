#include <iostream>
#include <cmath>
#include <fstream>

#include "Lattice.hh"
#include "SimulationMT.hh"
#include "Simulation.hh"

using namespace std;

#define beta 2.7
#define Nt 4
#define Ns 4
#define dim 4

int main(void){
    int sizes[4][2] = {{4,8},{6,12},{8,16},{10,20}};
    string modes[3] = {"hot_start", "cold_start", "warm_start"};
    for(size_t i = 1; i < 3; i++){
        SimulationMT* sim = new SimulationMT(beta, dim, sizes[0][0], sizes[0][1], modes[i], "wilson_loop", 3, true);
        cout << "Simulation begin: \n";
        sim->RunMT();
        delete sim;
    }
    
    return 0;
}