#ifndef SIMULATION_HH
#define SIMULATION_HH

#include "Link.hh"

class Simulation{
    public:
        Simulation();
        ~Simulation(); //  dim, Ns,  Nt
        Simulation(double, size_t, size_t, size_t, const std::string, std::string, bool);
        Simulation(double, size_t, size_t, size_t, const std::string, std::string, bool, bool);
        size_t StepMC();
        double GetAction() const;
        void MeasureObservable(std::string obsName);
        void Run();
        bool GetGo() const{ return go; }
    protected:
        size_t dim;
        size_t Nt;
        size_t Ns;
        size_t totalN;
        double beta;
        double spread;
        size_t time;
        bool go;
        bool verbose;
        Link* link;
        std::string actionFile;
        std::string configFile;
        std::string* observable_names;

        bool UpdateMC();
        void InitFiles(std::string* fileNames, size_t size) const;
        void WriteConfig() const;
        void WriteEDensity() const;
        Simulation& GetInstance(){ return *this;}
};

#endif