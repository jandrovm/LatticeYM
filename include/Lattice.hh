#ifndef LATTICE_HH
#define LATTICE_HH

#include "SU3.hh"
#include "utils.hh"
#include <cstring>

class Lattice{
    public:
        
        Lattice();
        ~Lattice();
        Lattice(size_t, size_t);
        Lattice(size_t, size_t, size_t);
        Lattice(const Lattice&);
        Lattice(int*, size_t, size_t, size_t);

        int* n;
        size_t GetNs() const { return Ns; }
        size_t GetNt() const { return Nt; }
        size_t GetDim() const { return dim; }    
        Lattice TransportMu(size_t , int);
        Lattice& MultiIndex(const size_t I); //Sort of constructor
        void PrintLattice();

        Lattice& operator=( const Lattice&);
        Lattice& operator+=( const Lattice&);
        Lattice& operator-=( const Lattice&);
        bool operator==( const Lattice&);

    private:
        size_t Ns;
        size_t Nt;
        size_t totalN;
        size_t dim;
        void AllocLattice();
};
Lattice operator+(const Lattice&, const Lattice&);
Lattice operator-(const Lattice&, const Lattice&);

size_t SingleIndex(const Lattice& m);
int Parity(const Lattice& m);
#endif