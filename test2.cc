#include "Link.hh"
#include "SU3.hh"

#include <iostream>
#include <cmath>
#include <cstring>
#include <random>

using namespace std;

#define Nt static_cast<size_t>(4u)
#define Ns static_cast<size_t>(4u)
#define dim static_cast<size_t>(4u)
#define beta 5.76

int main(void){
    double spread = 0.05;
    string mode = "hot_start";
    Link* link = new Link(mode, dim,  Ns, Nt, spread,false);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> ni(0,Ns-1);
    uniform_int_distribution<int> n4(0,Nt-1); 
    
    cout << "Lattice pointer declared succesfully!\n";
    cout << "\n" << "-----------------------------------------------------------" << "\n";
    cout << "\n" << "-----------------------------------------------------------" << "\n\n";
    //link->SetSpread(spread);
    cout << "Spread of initial SU(3) matrices around the identity correctly set!\n";
    cout << "\n" << "-----------------------------------------------------------" << "\n";
    cout << "\n" << "-----------------------------------------------------------" << "\n\n";
    Lattice latSite(dim,Ns,Nt);
    for(size_t mu = 0; mu < dim-1; ++mu){
        latSite.n[mu] = ni(gen);
    }
    latSite.n[dim-1] = n4(gen);
    string mode_output;
    if(mode == "cold_start"){ mode_output = "Cold Start: all matrices set to unity."; }
    else if(mode == "hot_start"){ mode_output = "Hot Start: matrices uniformly set in all SU(3) (QR factorisation)."; }
    else if(mode == "cold_start"){ mode_output = "Warm Start: matrices uniformly set around unity (spread of 0.005 around Id)."; }
    cout << "\n" << "-----------------------------------------------------------" << "\n";
    cout << "SU(3) link variables declared in "+mode_output+"\n";
    cout << Ns << "^" << dim-1 << "x" << Nt << "lattice with SU(3) gauge group.\n";
    cout << "\n" << "-----------------------------------------------------------" << "\n\n";
    SU3* X = new SU3[dim];

    cout << "Random lattice point: " << endl;
    latSite.PrintLattice();
    cout << "Links at that site in all " << dim << " directions." << endl;
    for(size_t mu = 0; mu < dim; ++mu){
        link->GetLink(mu,SingleIndex(latSite)).PrintSU3();
        cout << "With determinant: " << link->GetLink(mu,SingleIndex(latSite)).DetSU3() << "." << endl;
    }

    cout << "\n" << "-----------------------------------------------------------" << "\n";
    cout << "\n" << "-----------------------------------------------------------" << "\n\n";

    cout << "Small random SU(3) (close to I) link change matrix in all " << dim << " directions: " << endl;
    for(size_t mu = 0; mu < dim; ++mu){
        X[mu].SmallRandSU3(spread);
        X[mu].PrintSU3();
        cout << "With determinant: " << X[mu].DetSU3() << "." << endl;
    }

    cout << "\n" << "-----------------------------------------------------------" << "\n";
    cout << "\n" << "-----------------------------------------------------------" << "\n\n";
    link->ProposeChange(X,SingleIndex(latSite));
    cout << "Change in said link variables after a random transformation close to unity: \n";
    
    for(size_t mu = 0; mu < dim; ++mu){
        link->GetLink(mu,SingleIndex(latSite)).PrintSU3();
        cout << "With determinant: " << link->GetLink(mu,SingleIndex(latSite)).DetSU3() << "." << endl;
    }

    cout << "\n" << "-----------------------------------------------------------" << "\n";
    cout << "\n" << "-----------------------------------------------------------" << "\n\n";
    for(size_t n = 0; n < 10; ++n){
        SU3* staple = new SU3[dim];
        for(size_t mu = 0; mu < dim; ++mu){
            staple[mu] = link->GetStaple(latSite,mu);
            X[mu].SmallRandSU3(spread);
        }
        cout << "Action variation after U_mu(n)-->X_mu U_mu(n) with X close to unity: " << 
        link->ActionVariation(X,latSite,5.76,staple) << " .\n";
        delete [] staple;
    }
    
    delete [] X;
    delete link;
    return 0;
}