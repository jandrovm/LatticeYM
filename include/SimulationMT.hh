#ifndef SIMULATIONMT_HH
#define SIMULATIONMT_HH

#include <iostream>
#include <random>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <vector>
#include <thread>

#include "SimulationMT.hh"
#include "Simulation.hh"
#include "SU3.hh"
#include "utils.hh"
#include "Link.hh"
#include "Lattice.hh"
#include "ThreadPool.hh"

class SimulationMT : public Simulation {
    public:
        SimulationMT();
        ~SimulationMT(){};
        SimulationMT(double coupling, size_t dimension, size_t spaceN, size_t timeN, const std::string mode, std::string obs, size_t threads, bool verb);
        size_t StepMCMT();
        void RunMT();
    private:
        size_t nThreads;

        size_t UpdateMCMT(ThreadPool&, size_t);
        void CallableUpdate(Lattice&, size_t&);
};

#endif