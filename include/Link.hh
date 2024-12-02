#ifndef LINK_HH
#define LINK_HH

#include "SU3.hh"
#include "utils.hh"
#include "Lattice.hh"
#include <cstring>
#include <complex>

class Link{
    public:
        
        Link();
        ~Link();
        Link(const std::string, const size_t, const size_t);
        Link(const std::string, const size_t, const size_t, const size_t, double, bool);
        Link(SU3**, const size_t, const size_t, const size_t, double);
        Link(SU3**, const size_t, const size_t, const size_t, double, bool);
        void ProposeChange(SU3*, size_t);
        void ProposeChangeSingleLink(SU3&, size_t, size_t);
        double ActionVariation(SU3* X, Lattice n, const double beta, SU3* staple) const;
        double ActionVariationSingleLink(SU3& X, const Lattice& n, size_t mu, const double beta, SU3& staple) const;
        void SetSpread(double);
        size_t GetNs() const { return Ns; }
        size_t GetNt() const { return Nt; }
        size_t GetDim() const { return dim; }
        size_t GetPad() const { return pad; }
        SU3& GetLink(size_t mu, size_t I) const { return U[mu][I]; }
        SU3 GetStaple(Lattice n, size_t mu) const;
        double GetEDensity() const;
        double MeasureWilsonLoop(size_t , size_t, size_t, bool) const;
        Link& operator=(const Link&);
        

    private:
        size_t Ns;
        size_t Nt;
        size_t totalN;
        size_t dim;
        bool MT;
        size_t pad;
        SU3** U;
        //For the small random SU3 matrices, around unity
        double* initSpread = NULL;
        
        const std::string mode;
        std::complex<double> WilsonLoop(size_t , size_t, size_t) const;
        void AllocLink();
};

#endif