#include <iostream>
#include <random>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <sstream>

#include "Simulation.hh"
#include "SU3.hh"
#include "utils.hh"
#include "Link.hh"
#include "Lattice.hh"

using namespace std;

//Until I find a better algorithm to determine the thermalization
#define nTerm 5000
#define nMC 500
#define nUp 10 //Local updates to reduce correlation
#define decorr 50 //Sweeps between configurations

Simulation::Simulation(){
    spread = 0.05;
    beta = 5.67;
    time = 0;
    dim = 4;
    Ns = 4;
    Nt = 4;
    totalN = Ns*Ns*Ns*Nt;
    verbose = true;
    cout << "Initializing lattice: " << 4 << "^" << 3 << "x" << 4 << " sites.\n";
    go = true;

    actionFile = "hot_start_action.dat";
    configFile = "config.dat";
    string* fileNames = new string[2];
    fileNames[0] = actionFile;
    fileNames[1] = configFile;
    InitFiles(fileNames,2);
    delete [] fileNames;

    link = new Link("hot_start",4,4,4,spread,false);
    observable_names = new string;
    *observable_names = "wilson_loop";
    cout << "Lattice initialized! \n";
}
Simulation::Simulation(double coupling, size_t dimension, size_t spaceN, size_t timeN, const string mode, string obs, bool verb){
    spread = 0.01;
    beta = coupling;
    time = 0;
    verbose = verb;
    Nt = timeN;
    Ns = spaceN;
    dim = dimension;
    totalN = 1;

    for(size_t mu = 0; mu < dim-1; ++mu){ totalN *= Ns; }
    totalN *= Nt;
    if(Nt%2==1){
        Nt += 1;
    }
    if(Ns%2==1){
        Ns += 1;
    }
    cout << "Initializing lattice: " << Ns << "^" << dim-1 << "x" << Nt << " sites.\n";
    link = new Link(mode,dim, Ns, Nt, spread,false);
    
    actionFile = mode+"_action.dat";
    configFile = "config.dat";
    string* fileNames = new string[2];
    fileNames[0] = actionFile;
    fileNames[1] = configFile;
    InitFiles(fileNames,2);
    delete [] fileNames;
    
    go = true;
    observable_names = new string;
    *observable_names = obs;

    cout << "Lattice initialized! \n";
}
//This next one is only for MT mode, where padding must be implemented to avoid false sharing
Simulation::Simulation(double coupling, size_t dimension, size_t spaceN, size_t timeN, const string mode, string obs, bool verb, bool padding){
    spread = 0.01;
    beta = coupling;
    time = 0;
    verbose = verb;
    Nt = timeN;
    Ns = spaceN;
    dim = dimension;
    totalN = 1;

    for(size_t mu = 0; mu < dim-1; ++mu){ totalN *= Ns; }
    totalN *= Nt;
    if(Nt%2==1){
        Nt += 1;
    }
    if(Ns%2==1){
        Ns += 1;
    }
    cout << "Initializing lattice: " << Ns << "^" << dim-1 << "x" << Nt << " sites.\n";
    link = new Link(mode,dim, Ns, Nt, spread,padding);
    
    actionFile = mode+"_action.dat";
    configFile = "config.dat";
    string* fileNames = new string[2];
    fileNames[0] = actionFile;
    fileNames[1] = configFile;
    InitFiles(fileNames,2);
    delete [] fileNames;
    
    go = true;
    observable_names = new string;
    *observable_names = obs;

    cout << "Lattice initialized! \n";
}
Simulation::~Simulation(){
    delete link;
    delete observable_names;
}

//-------------------------------------------------------------------
//---------------------------Update methods--------------------------
//-----------------------------(MonteCarlo)--------------------------
//-------------------------------------------------------------------

bool Simulation::UpdateMC(){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> p(0.,1.);
    uniform_int_distribution<size_t> I(0,totalN-1);
    size_t site = I(gen);
    
    Lattice* siteLat = new Lattice(dim,Ns,Nt);
    siteLat->MultiIndex(site);

    SU3* staple = new SU3[dim];
    SU3* X = new SU3[dim];
    for(size_t mu = 0; mu < dim; ++mu){
        staple[mu]=link->GetStaple(*siteLat,mu);
    }

    bool accepted = false;
    for(size_t n = 0; n < nUp; ++n){
        for(size_t mu = 0; mu < dim; ++mu){
            X[mu].SmallRandSU3(spread);
        }  
        double dS = link->ActionVariation(X,*siteLat,beta,staple);
        double prob = min(1.,exp(-dS));
        //cout << "exp(-dS)=" << exp(-dS) << "\n";
        if( prob >= p(gen) ){
            //We accept change
            link->ProposeChange(X,site);
            accepted = true;
        }
        else accepted = false;
    }
    delete [] X;
    delete [] staple;
    delete siteLat;
    //Else link is not updated
    return accepted;
}
size_t Simulation::StepMC(){
    if(!go){
        cout << "Simulation end.\n";
        return 0;
    }
    bool accept;
    for(size_t n = 0; n < totalN; ++n){
        accept = UpdateMC();
    }
    /*fstream file;
    file.open(actionFile, ios::in | ios::out | ios::app);
    double action = link->GetEDensity();
    if(!file){
        cout << "Error outputting action data!\n";
    }
    else file << time << setw(15) << action << "\n";
    if(verbose){cout << time << setw(15) << action << "\n";}
    file.close();*/

    if(time == nTerm + decorr*nMC){ go = false; }
    time++;
    if(accept){
        return 1;
    }
    else return 0;
}
void Simulation::Run(){
    string name;
    size_t acceptance = 0;
    if(*observable_names == "wilson_loop"){
        name = "Wilson loop";
    }
    else name = "Undefined observable";
    for(size_t n = 0; n < nTerm; ++n){
        StepMC();
        if(verbose){
            cout << "Thermalization: " << (float)100.*(n+1.)/(nTerm) << "%.\n";
        }
    }
    for(size_t n = 0; n < nMC; ++n){
        for(size_t m = 0; m < decorr; ++m){
            acceptance += StepMC();
        }
        WriteConfig(); 
        if(verbose){
            cout << "Generating configurations: " << (float)100.*(n+1)/(nMC) << "%.\n";
        }
    }
    go = false;
    cout << "Acceptance ratio: " << 100.*(float)acceptance/(decorr*nMC) << "%\n";
    MeasureObservable(*observable_names);

    return;
}
void Simulation::InitFiles(string* fileNames, size_t size) const{
    for(size_t i = 0; i < size; i++){
        fstream file;
        file.open(fileNames[i], ios::out | ios::trunc);
        if(!file){
            cout << "Error: Output file '"+fileNames[i]+"' not initialized!\n";
        }
        file.close();
    }
    return;
}

//---------------------------------------------------------------------
//---------------------Measurement and observable methods--------------
//---------------------------------------------------------------------

double Simulation::GetAction() const{
    SU3 UMuNu;

    Lattice n(dim,Ns,Nt);
    size_t IPlusMu, IPlusNu;
    double sum = 0.;
    for(size_t I = 0; I < totalN; ++I){
        n.MultiIndex(I);
        for(size_t mu = 0; mu < dim; ++mu){
            IPlusMu = SingleIndex(n.TransportMu(mu,1));
            for(size_t nu = mu+1; nu < dim; ++nu){
                IPlusNu = SingleIndex(n.TransportMu(nu,1));
                UMuNu = link->GetLink(mu,I)*link->GetLink(mu,IPlusMu)*(link->GetLink(mu,IPlusNu)).Adj()*(link->GetLink(mu,I)).Adj();
                sum += ((SU3(1.)-UMuNu).TraceSU3()).real();
            }
        }
    }
    return sum/3.;
}


//-----------------------------------------------------------------------------------
//----------------------------Data output/readout methods----------------------------
//-----------------------------------------------------------------------------------

void Simulation::WriteConfig() const{
    
    fstream file;
    file.open(configFile, ios::out | ios::app);
    //In each line we write U_0(n), U_1(n), U_2(n), U_3(n)
    //Each U_\mu(n) is written as (U_\mu(n))_{ij} in real and imag parts
    size_t pad = link->GetPad();
    for(size_t I = 0; I < totalN; ++I){
        //Each line containts 2x4x9 elements: 4 lorentz components, 9 matrix elements, 2 from complex numbers
        for(size_t mu = 0; mu < dim; ++mu){
            for(size_t i = 0; i < 3; ++i){
                for(size_t j = 0; j < 3; ++j){
                    file << link->GetLink(mu,pad*I)(i,j).real() << setw(15) << link->GetLink(mu,pad*I)(i,j).imag() << setw(15);
                }
            }
        }
        file << "\n";
    }
    file.close();
    return;
}
void Simulation::WriteEDensity() const{
    fstream file;
    file.open(actionFile, ios::in | ios::out | ios::app);
    double action = link->GetEDensity();
    if(!file){
        cout << "Error outputting action data!\n";
    }
    else file << time << setw(15) << action << "\n";
    if(verbose){cout << time << setw(15) << action << "\n";}
    file.close();
}
void Simulation::MeasureObservable(string obsName){
    //The idea is that meanObs (mean of observable) is computed in simulation time
    ifstream file;
    file.open(configFile, ios::in);
    if(!file){ cout << "Error measuring from config file.\n"; return;}
    
    SU3** U = new SU3*[dim];
    Lattice lat(dim,Ns,Nt);
    for(size_t mu = 0; mu < dim; ++mu){
        mu == dim-1 ? (lat.n[mu] = Nt/2) : (lat.n[mu] = Ns/2); 
    }
    for(size_t mu = 0; mu < dim; ++mu){
        U[mu] = new SU3[link->GetPad()*totalN];
    }
    size_t II = SingleIndex(lat);
    //We overdid one
    double* O = new double[nMC];
    size_t counter = 0;
    string meanName, varName;
    complex<double>** M = new complex<double>*[3];
    for(size_t i = 0; i < 3; i++){
        M[i] = new complex<double>[3]; 
    }
    //We need to read the file first
    for(size_t n = 0; n < nMC; n++){
        string s;
        for(size_t I = 0; I < totalN && getline(file,s); ++I){
            istringstream sin(s);
            for(size_t mu = 0; mu < dim; ++mu){
                for(size_t i = 0; i < 3; i++){
                    for(size_t j = 0; j < 3; j++){
                        double re, im;
                        //sin works like cin
                        sin >> re >> im;
                        M[i][j] = complex<double>(re,im);
                    }
                }
                U[mu][link->GetPad()*I] = SU3(M,false);
            }
            //cout << "Observable computation: " << (float)100.*(n*totalN + I)/(nMC*totalN-1.) << "%\n";
        }
        Link* L = new Link(U,dim,Ns,Nt,spread,true);
        if(*observable_names == "wilson_loop"){
            meanName = "<W>";
            varName = "Var[W]";
            O[n] = L->MeasureWilsonLoop(II,2,Nt-2,true);
        }
        delete L;
    }
    file.close();
    for(size_t i = 0; i < 3; i++){
        delete [] M[i];
    }
    delete [] M;
    double meanO = mean(O,nMC);
    cout << meanName+" = " << meanO << ", "+varName+" = " << sqrt(variance(O,meanO,nMC)) << ".\n";
    for(size_t mu = 0; mu < dim; ++mu){
        delete [] U[mu];
    }
    delete [] U;
    return;
}