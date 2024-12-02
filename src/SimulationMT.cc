#include "SimulationMT.hh"

using namespace std;

//Until I find a better algorithm to determine the thermalization
#define nTerm 1000
#define nMC 1000
#define nUp 10 //Local updates to reduce correlation
#define decorr 25 //Sweeps between configurations

SimulationMT::SimulationMT() :
    Simulation(){
    Simulation(5.76,4,4,4,"hot_mode","wilson_loop",true,true);
    nThreads = 2;
    cout << "Threads: " << nThreads << "\n";
}
SimulationMT::SimulationMT(double coupling, size_t dimension, size_t spaceN, size_t timeN, const string mode, string obs, size_t threads, bool verb) : 
Simulation(coupling,dimension,spaceN,timeN,mode,obs,verb,true){  
    nThreads = threads;
    cout << "Threads: " << nThreads << "\n";
}

//-------------------------------------------------------------------
//--------------------Multithread Update methods---------------------
//-----------------------------(MonteCarlo)--------------------------
//-------------------------------------------------------------------
void SimulationMT::CallableUpdate(Lattice& siteLat, size_t& acceptCount){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> p(0.,1.);
    SU3* staple = new SU3[dim];
    SU3* X = new SU3[dim];
    for(size_t mu = 0; mu < dim; mu++){
        *staple = this->link->GetStaple(siteLat,mu);
        
    }
    //cout << "Thread ID: " << this_thread::get_id() << endl;
    
    bool accepted = false;
    acceptCount = 0;
    for(size_t n = 0; n < nUp; ++n){    
        for(size_t mu = 0; mu < dim; ++mu){
            X[mu].SmallRandSU3(spread);
        }
        double dS = this->link->ActionVariation(X,siteLat,beta,staple);
        double prob = min(1.,exp(-dS));
        //cout << "exp(-dS)=" << exp(-dS) << "\n";
        size_t site = SingleIndex(siteLat);
        if( prob >= p(gen) ){
            //We accept change
            this->link->ProposeChange(X,site);
            acceptCount++;
        }
    }   
    delete [] X;
    delete [] staple;
    return;
}
size_t SimulationMT::UpdateMCMT(ThreadPool& workerThreads, size_t threadNum){
    size_t acceptCount = 0;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> I(0,static_cast<int>(totalN)-1);

    //Since the wilson action contains only nn interactions, we must update only sites
    //of equal parity, i.e, only half the lattice available for parallel updates 
    size_t site = I(gen); 
    
    Lattice* siteLat = new Lattice[threadNum];
    siteLat[0] = Lattice(dim,Ns,Nt);
    siteLat[0].MultiIndex(site);
    
    int* n = new int[dim];
    for(size_t i = 1; i < threadNum; ++i){
        //This goes shifting the other sites about the first one like:
        //site[1].n[0] = site[0].n[0]+2 and the rest of coordinates equal,
        //site[2].n[1] = site[0].n[1]+2 and the rest of coordinates equal,
        // ...
        //site[dim].n[dim-1] = site[0].n[dim-1]+2 and the rest of coordinates equal,
        //site[dim+1].n[0] = site[0].n[0]+4 and the rest of coordinates equal,
        //site[dim+2].n[1] = site[0].n[1]+4 and the rest of coordinates equal,
        // ...
        //site[2*dim].n[dim-1] = site[0].n[dim-1]+4 and the rest of coordinates equal,
        //site[2*dim+1].n[0] = site[0].n[0]+6 and the rest of coordinates equal, etc.
        for(size_t mu = 0; mu < dim-1; ++mu){
            if(mu == (i-1)%dim){
                n[mu] = siteLat[0].n[mu] + 2*((i-1)/dim+1);
                n[mu] %= Ns;
            }
            else n[mu] = siteLat[0].n[mu];
        }
        if( dim-1 == (i-1)%dim){
            n[dim-1] = siteLat[0].n[dim-1] + 2*((i-1)/dim+1);
            n[dim-1] %= Nt;
        }
        else n[dim-1] = siteLat[0].n[dim-1];

        siteLat[i] = Lattice(n,dim,Ns,Nt);
    }
    delete [] n;
    size_t* count = new size_t[threadNum];
    for(size_t i = 0; i < threadNum; ++i){
        count[i] = 0;
        workerThreads.Enqueue( bind(&SimulationMT::CallableUpdate, this, siteLat[i], count[i]) );
    }
    
    for(size_t i = 0; i < threadNum; ++i){
        acceptCount += count[i];
    }
    delete [] count;
    delete [] siteLat;
    //Else link is not updated
    assert(acceptCount <= threadNum*nUp);
    return acceptCount;
}

size_t SimulationMT::StepMCMT(){
    if(!go){
        cout << "SimulationMT end.\n";
        return 0;
    }
    size_t accept = 0;
    ThreadPool workers(nThreads);
    for(size_t n = 0; n < totalN/nThreads; ++n){
        accept += UpdateMCMT(workers,nThreads);
    }
    size_t remaining = totalN%nThreads;
    if( remaining != 0 ){
        accept += UpdateMCMT(workers,remaining);
    }
    assert(accept <= nUp*totalN);

    if(time == nTerm + decorr*nMC){ go = false; }
    time += 1;
    return accept;
}
void SimulationMT::RunMT(){
    
    string name;
    size_t acceptance = 0;
    if(*observable_names == "wilson_loop"){
        name = "Wilson loop";
    }
    else name = "Undefined observable";
    for(size_t n = 0; n < nTerm; ++n){
        StepMCMT();
        if(verbose){
            cout << "Thermalization: " << (float)100.*(n+1.)/(nTerm) << "%.\n";
        }
    }
    for(size_t n = 0; n < nMC; ++n){
        for(size_t m = 0; m < decorr; ++m){
            acceptance += StepMCMT();
        }
        WriteConfig(); 
        if(verbose){
            cout << "Generating configurations: " << (float)100.*(n+1)/(nMC) << "%.\n";
        }
        WriteEDensity();
    }
    go = false;
    cout << "Acceptance ratio: " << 100.*(float)acceptance/(nUp*decorr*nMC*totalN) << "%\n";
    MeasureObservable(*observable_names);

    return;
}
