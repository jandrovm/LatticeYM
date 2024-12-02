#include <cmath>
#include <iostream>
#include <cstring>

#include "Lattice.hh"
#include "utils.hh"

using namespace std;

void Lattice::AllocLattice(){
    n = new int[dim];
    return;
}
Lattice::Lattice(){
    Ns = Nt = 1u;
    dim = 1u;
    AllocLattice();
    totalN = 1u;
    n[0] = 0;
}
Lattice::~Lattice(){
    delete [] n;
}
Lattice::Lattice(size_t di, size_t nS, size_t nT){
    
    dim = di;
    Ns = nS;
    Nt = nT;
    totalN = Nt;
    for(size_t d = 0; d < dim-1; ++d){ totalN *= Ns; } 

    AllocLattice();
    for(size_t mu = 0; mu < dim; ++mu){
        n[mu] = 0;
    }
}
Lattice::Lattice(int* m, const size_t di, const size_t nS, const size_t nT){
    
    dim = di;
    Ns = nS;
    Nt = nT;
    totalN = Nt;
    for(size_t d = 0; d < dim-1; ++d){ totalN *= Ns; } 

    AllocLattice();
    for(size_t mu = 0; mu < dim; ++mu){
        while(m[mu] < 0){ mu == dim-1u ? m[mu] += Nt : m[mu] += Ns; }
        mu == dim-1u ? m[mu] %= Nt : m[mu] %= Ns;
        n[mu] = m[mu];
    }
}
Lattice::Lattice(size_t nS, size_t nT){
    
    dim = 4u;
    Ns = nS;
    Nt = nT;
    totalN = 1u;
    for(size_t d = 0; d < dim-1u; ++d){ totalN *= Ns; }
    totalN *= Nt; 

    AllocLattice();
    for(size_t mu = 0; mu < dim; ++mu){
        n[mu] = 0;
    }
}
Lattice::Lattice(const Lattice& m){
    
    Ns = m.Ns;
    Nt = m.Nt;
    dim = m.dim;
    AllocLattice();
    totalN = m.totalN;
    for(size_t mu = 0; mu < dim; ++mu){
        n[mu] = m.n[mu];
    }
}
Lattice Lattice::TransportMu(size_t mu, int sgn){
    int sgnTrue = sgn/abs(sgn); //Just in case
    Lattice temp = *this;
    if(sgnTrue == 1){
        temp.n[mu]++;
        if( mu == dim-1u ){
            temp.n[mu] = temp.n[mu]%Nt;
        }
        else{ temp.n[mu] = temp.n[mu]%Ns; }
    }
    else if(sgnTrue == -1){
        temp.n[mu]--;
        if( mu == dim-1u ){
            while( temp.n[mu] < 0 ){ temp.n[mu] += static_cast<int>(Nt); }
            temp.n[mu] = temp.n[mu]%static_cast<int>(Nt);
        }
        else{ 
            while( temp.n[mu] < 0 ){ temp.n[mu] += static_cast<int>(Ns); } 
            temp.n[mu] = temp.n[mu]%static_cast<int>(Ns); 
        }
    }
    return temp;
}
Lattice& Lattice::MultiIndex(const size_t I){ 
    Lattice* temp = this; 
    ijkn(temp->n, Ns, I, dim); 
    return *temp; 
}
void Lattice::PrintLattice(){
    cout << "(";
    for(size_t d = 0; d < dim-1u; ++d){
        cout << n[d] << ", ";
    }
    cout << n[dim-1u] << ")\n";
    return;
}
Lattice& Lattice::operator=(const Lattice& m){
    if( this == &m ){
        return *this;
    }
    dim = m.dim;
    Ns = m.Ns;
    Nt = m.Nt;
    totalN = m.totalN;
    for(size_t mu = 0; mu < dim; ++mu){
        n[mu] = m.n[mu];
    }
    return *this;
}
Lattice& Lattice::operator+=(const Lattice& m){
    for(size_t mu = 0; mu < dim-1u; ++mu){
        n[mu] += m.n[mu];
        while( n[mu] < 0 ){ n[mu] += static_cast<int>(Ns); }
        n[mu] = n[mu]%static_cast<int>(Ns);
    }
    n[dim-1u] += m.n[dim-1u];
    while( n[dim-1u] < 0 ){ n[dim-1] += static_cast<int>(Nt); }
    n[dim-1u] = n[dim-1u]%static_cast<int>(Nt);
    return *this;
}
Lattice& Lattice::operator-=(const Lattice& m){
    //Negative numbers give a different modulo than positive numbers
    //different ordering in Z_n
    for(size_t mu = 0; mu < dim-1u; ++mu){
        n[mu] -= m.n[mu];
        while(n[mu] < 0){ n[mu] += static_cast<int>(Ns); }
        n[mu] = n[mu]%static_cast<int>(Ns);

    }
    n[dim-1u] -= m.n[dim-1u];
    while(n[dim-1u] < 0){ n[dim-1u] += static_cast<int>(Nt); }
    n[dim-1u] = n[dim-1u]%static_cast<int>(Nt);
    return *this;
}
bool Lattice::operator==(const Lattice& m){
    bool equal = (dim == m.dim) && (Ns == m.Ns) && (Nt == m.Nt);
    for(size_t mu = 0; mu < dim; mu++){
        equal = equal && (n[mu] == m.n[mu]);
    }
    return equal;
}
Lattice operator+(const Lattice& m1, const Lattice& m2){
    Lattice temp(m1);
    return (temp += m2);
}
Lattice operator-(const Lattice& m1, const Lattice& m2){
    Lattice temp(m1);
    return (temp -= m2);
}

size_t SingleIndex(const Lattice& m){ 
    return I(m.n, m.GetNs(), m.GetDim()); 
}
int Parity(const Lattice& m){
    //In a checkerboard pattern returns 0 for positive parity ('even' sites) and 1 for negative ('odd' sites)
    int dim = m.GetDim();
    int even = 0;
    for(int mu = 0; mu < dim; ++mu){
        even = even ^ (m.n[mu]%2);
    }
    return even;
}