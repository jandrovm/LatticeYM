#include <cmath>
#include <iostream>
#include <cstring>

#include "Link.hh"
#include "Lattice.hh"
#include "SU3.hh"
#include "utils.hh"

using namespace std;

#define PAD_STD 2u

void Link::AllocLink(){
    U = new SU3*[dim];
    for(size_t mu = 0; mu < dim; ++mu){
        U[mu] = new SU3[totalN];
    }
}
Link::Link(){
    Ns = Nt = 1;
    dim = 1;
    initSpread = NULL;
    AllocLink();
    totalN = 1;
    MT = false;
    pad = 1;
    U[0][0] = SU3(1.);
    return;
}
Link::Link(SU3** V, const size_t di, const size_t nS, const size_t nT, double spread){
    Nt = nT; 
    Ns = nS; 
    dim = di;
    totalN = 1;
    initSpread = &spread;
    MT = false;
    pad = 1;
    for(size_t mu = 0; mu < dim; ++mu){mu < dim-1 ? totalN *= Ns : totalN *= Nt;}
    AllocLink();
    for(size_t mu = 0; mu < dim; ++mu){
        for(size_t I = 0; I < totalN; ++I){
            U[mu][I] = V[mu][I];
        }
    }
}
Link::Link(SU3** V, const size_t di, const size_t nS, const size_t nT, double spread, bool multiThread){
    Nt = nT; 
    Ns = nS; 
    dim = di;
    totalN = 1;
    initSpread = &spread;
    MT = multiThread;
    if(MT){
        pad = PAD_STD;
    }
    else pad = 1;
    for(size_t mu = 0; mu < dim; ++mu){mu < dim-1 ? totalN *= Ns : totalN *= Nt;}
    totalN *= pad;
    AllocLink();
    totalN /= pad;
    for(size_t mu = 0; mu < dim; ++mu){
        for(size_t I = 0; I < totalN; ++I){
            U[mu][pad*I] = V[mu][pad*I];
        }
    }
}
Link::~Link(){
    for(size_t mu = 0; mu < dim; ++mu){
        delete [] U[mu];
    }
    delete [] U;
}
Link::Link(const string mode, const size_t di, const size_t nS, const size_t nT, double spread, bool multiThread){
    //Hot start: all matrices initialized randomly via Haar measure
    //Cold start: all matrices set to unity
    //Warm start: all matrices random close to unity
    MT = multiThread;
    if(MT){
        pad = PAD_STD;
    }
    else pad = 1;

    bool correctMode = (mode=="hot_start")||(mode=="cold_start")||(mode=="warm_start");
    assert(correctMode);
    dim = di;
    Ns = nS;
    Nt = nT;
    cout << "Ns=" << Ns << ", Nt=" << Nt << ", dim=" << dim-1 << "+" << 1 << ".\n";
    totalN = Nt;
    for(size_t d = 0; d < dim-1; ++d){ totalN *= Ns; }
    cout << "\n"; 
    cout << "Total number of link variables U_\\mu: " << totalN << "\n";
    cout << "\n";
    initSpread = &spread;
    //Just to allocate space, we need to undo this or else we get double the time dimension
    totalN *= pad;
    AllocLink();
    totalN /= pad;

    cout << "Initializing link variables in "+mode+ ".\n";
    for(size_t II = 0; II < totalN; ++II){
        for(size_t mu = 0; mu < dim; ++mu){
            if( mode == "hot_start" ){
                U[mu][pad*II].RandSU3QR();
            }
            else if( mode == "cold_start" ){
                U[mu][pad*II] = SU3(1.);
            }
            else if( mode == "warm_start" ){
                U[mu][pad*II].SmallRandSU3(*initSpread);
            }
        }
    }
    return;
}
Link::Link(const string mode, const size_t nS, const size_t nT){
    //Hot start: all matrices initialized randomly via Haar measure
    //Cold start: all matrices set to unity
    //Warm start: all matrices random close to unity

    bool correctMode = (mode=="hot_start")||(mode=="cold_start")||(mode=="warm_start");
    assert(correctMode);
    dim = 4u;
    Ns = nS;
    Nt = nT;
    totalN = nS*nS*nS*nT;
    MT = false;
    pad = 1;
    if(!initSpread){ *initSpread = 0.01; } 
    AllocLink();
    
    for(size_t II = 0; II < totalN; ++II){
        for(size_t mu = 0; mu < dim; ++mu){
            if( mode == "hot_start" ){
                U[mu][II].RandSU3QR();
            }
            else if( mode == "cold_start" ){
                U[mu][II]=SU3(1.);
            }
            else if( mode == "warm_start" ){
                U[mu][II].SmallRandSU3(*initSpread);
            }
        }
    }
    return;
}

void Link::SetSpread(double s){
    initSpread = &s;
    return;
}
void Link::ProposeChange(SU3* X, size_t I){
    assert( I >= 0 && I < totalN );
    for(size_t mu = 0; mu < dim; ++mu){
        U[mu][pad*I] = X[mu]*U[mu][pad*I];
        //U[mu][I].AssertDeterminant();
    }
    return;
}
void Link::ProposeChangeSingleLink(SU3& X, size_t I, size_t mu){
    assert( I >= 0 && I < totalN );
    assert( mu >= 0 && mu < dim );
    U[mu][pad*I] = X*U[mu][pad*I];
        //U[mu][I].AssertDeterminant();
    return;
}
//-------------------------------------------------------------------------------------
//-----------------------Observable-related methods------------------------------------
//-------------------------------------------------------------------------------------
double Link::ActionVariation(SU3* X, Lattice n, const double beta, SU3* staple) const {
    SU3* temp = new SU3(0.);
    size_t I = SingleIndex(n);
    //assert( I >= 0);
    //assert( I < totalN );
    //assert(n.GetDim() == dim && n.GetNs() == Ns && n.GetNt() == Nt );
    
    for(size_t mu = 0; mu < dim; ++mu){
        *temp += (X[mu]-SU3(1.))*U[mu][pad*I]*staple[mu];
        //cout << "U_mu(n)_{11} = " << U[mu][pad*I](0,0) << ", Staple_{11} = " << staple[mu](0,0) << endl; 
    }
    double dS = -beta/3.*(temp->TraceSU3()).real();
    delete temp;
    //cout << "Action variation: " << dS << endl;
    return dS;
}

double Link::ActionVariationSingleLink(SU3& X, const Lattice& n, size_t mu, const double beta, SU3& staple) const {
    SU3* temp = new SU3(0.);
    size_t I = SingleIndex(n);
    //assert( I >= 0);
    //assert( I < totalN );
    //assert(n.GetDim() == dim && n.GetNs() == Ns && n.GetNt() == Nt );
    
    //These act on the link U_mu only
    *temp += (X-SU3(1.))*U[mu][I]*staple;
    double dS = -beta/3.*(temp->TraceSU3()).real();
    delete temp;

    return dS;
}
SU3 Link::GetStaple(Lattice n, size_t mu) const{
    SU3 staple(0.);
    Lattice tempLat(dim, Ns, Nt);
    size_t nPlusMu, nPlusNu, nMinusNu, nPlusMuMinusNu;
    size_t I = SingleIndex(n);
    for(size_t nu = mu+1u; nu < dim; ++nu){
        nPlusMu = SingleIndex(n.TransportMu(mu,+1));
        nPlusNu = SingleIndex(n.TransportMu(nu,+1));
        nMinusNu = SingleIndex(n.TransportMu(nu,-1));
        tempLat = n.TransportMu(mu,+1);
        nPlusMuMinusNu = SingleIndex(tempLat.TransportMu(nu,-1));
        assert(nPlusMu < totalN && nPlusNu < totalN && nPlusMuMinusNu < totalN && nMinusNu < totalN);
        assert(nPlusMu >= 0 && nPlusNu >= 0 && nPlusMuMinusNu >= 0 && nMinusNu >= 0 );
        staple += U[nu][pad*nPlusMu]*U[mu][pad*nPlusNu].Adj()*U[nu][pad*I].Adj() 
        + U[nu][pad*nPlusMuMinusNu].Adj()*U[mu][pad*nMinusNu].Adj()*U[nu][pad*nMinusNu];
    }
    return staple;
}

complex<double> Link::WilsonLoop(size_t II, size_t R, size_t T) const {
    //Needs changing to avoid illegal site access (negative numbers or outside lattice)
    assert(II >=0 && II < totalN );
    assert(Ns/4 != 0 && Nt/4 != 0 && Ns%2 == 0 && Nt%2 == 0);
    assert( R < Ns && T < Nt);
    int k0 = 1;
    int n0 = 1;
    int k1 = k0 + R;
    int n1 = n0 + T;
    size_t d = dim;
    SU3 W(1.);
    // U_\mu(n)U_\nu(n+\mu)U_\mu^\dag(n+\nu)U_\nu^\dag(n)
    int* i = new int[dim];
    ijkn(i,Ns,II,dim);
    i[dim-2] = k0;
    i[dim-1] = n0;
    Lattice m(i,dim,Ns,Nt), nextm(m); 
    //Careful, transport routine returns a copy, it doesnt change the Lattice variable
    for(int k = k0; k < k1; ++k){ W *= U[d-2][pad*SingleIndex(nextm)]; nextm = nextm.TransportMu(dim-2,+1); }
    for(int n = n0; n < n1; ++n){ W *= U[d-1][pad*SingleIndex(nextm)]; nextm = nextm.TransportMu(dim-1,+1); }
    for(int k = k1; k > k0; --k){ nextm = nextm.TransportMu(dim-2,-1); W *= U[d-2][pad*SingleIndex(nextm)].Inv(); }
    for(int n = n1; n > n0; --n){ nextm = nextm.TransportMu(dim-1,-1); W *= U[d-1][pad*SingleIndex(nextm)].Inv(); }

    delete [] i;
    return 1./3.*W.TraceSU3();
}
double Link::MeasureWilsonLoop(size_t II, size_t R, size_t T, bool go) const {
    if(go){
        return WilsonLoop(II,R,T).real();
    }
    else return 0;
}
double Link::GetEDensity() const{
    SU3 UMuNu1, UMuNu2, UMuNu3, UMuNu4;
    Lattice n(dim,Ns,Nt);
    size_t IPlusMu, IPlusNu, IMinusNu, IMinusMu, IPlusMuMinusNu, IPlusNuMinusMu, IMinusMuMinusNu;
    Lattice tempPlus(n),tempMinus(n);
    double sum = 0.;
    for(size_t I = 0; I < totalN; ++I){
        n.MultiIndex(I);
        SU3 CMuNu;
        for(size_t mu = 0; mu < dim; ++mu){
            tempPlus = n.TransportMu(mu,+1);
            tempMinus = n.TransportMu(mu,-1);
            IPlusMu = SingleIndex(tempPlus);
            IMinusMu = SingleIndex(tempMinus);
            for(size_t nu = mu+1; nu < dim; ++nu){
                IPlusNu = SingleIndex(n.TransportMu(nu,1));
                IMinusNu = SingleIndex(n.TransportMu(nu,-1));
                IPlusMuMinusNu = SingleIndex(tempPlus.TransportMu(nu,-1));
                IPlusNuMinusMu = SingleIndex(tempMinus.TransportMu(nu,+1));
                IMinusMuMinusNu = SingleIndex(tempMinus.TransportMu(nu,-1));
                //Four 'leafs' of the clover term
                UMuNu1 = U[mu][pad*I]*U[nu][pad*IPlusMu]*U[mu][pad*IPlusNu].Adj()*U[nu][pad*I].Adj();
                UMuNu2 = U[mu][pad*I]*U[nu][pad*IPlusMuMinusNu].Adj()*U[mu][pad*IMinusNu].Adj()*U[nu][pad*IMinusNu];
                UMuNu3 = U[mu][pad*IMinusMu].Adj()*U[nu][pad*IMinusMu]*U[mu][pad*IPlusNuMinusMu]*U[nu][pad*I].Adj();
                UMuNu4 = U[mu][pad*IMinusMu].Adj()*U[nu][pad*IMinusMuMinusNu].Adj()*U[mu][pad*IMinusMuMinusNu]*U[nu][pad*IMinusNu];
                CMuNu = UMuNu1 + UMuNu3 - UMuNu2 - UMuNu4;
                CMuNu = (CMuNu - CMuNu.Adj() - 1./3.*(CMuNu - CMuNu.Adj()).TraceSU3()*SU3(1.));
                sum += ((CMuNu*CMuNu).TraceSU3()).real();
            }
        }
    }
    assert(totalN != 0);
    return -sum/(64.*totalN);
}
//----------------------------------------------------------------
//--------Operator overloads (only operator=() in this case)------
//----------------------------------------------------------------
Link& Link::operator=(const Link& link){
    Nt = link.Nt;
    Ns = link.Ns;
    dim = link.dim;
    totalN = link.totalN;
    initSpread = link.initSpread;
    if(this == &link){
        return *this;
    }
    else for(size_t I = 0; I < totalN; I++){
        for(size_t mu = 0; mu < dim; mu++){
            U[mu][pad*I] = link.U[mu][pad*I]; 
        }
    }
    return *this; 
}
