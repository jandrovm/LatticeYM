#include <cmath>
#include <complex>
#include <random>
#include <iomanip>
#include <iostream>

#include "SU3.hh"
#include "utils.hh"

using namespace std;

//---------------------------------------------------------------------
//---------------------------Constructors------------------------------
//---------------------------------------------------------------------

void SU3::AllocSpace(){
    U = new complex<double>*[3];
    for(size_t i = 0; i < 3; ++i){
        U[i] = new complex<double>[3];
    }
    return;
}
SU3::SU3(){
    AllocSpace();
    U[0][0]=U[1][1]=U[2][2]=complex<double>(1.,0.);
    U[0][1]=U[1][2]=U[0][2]=complex<double>(0.,0.);
    U[1][0]=U[2][1]=U[2][0]=complex<double>(0.,0.);
    return;
}
SU3::SU3(const double a){
    AllocSpace();
    U[0][0]=U[1][1]=U[2][2]=complex<double>(a,0.);
    U[0][1]=U[1][2]=U[0][2]=complex<double>(0.,0.);
    U[1][0]=U[2][1]=U[2][0]=complex<double>(0.,0.);
    return;
}
SU3::SU3(complex<double>* u, complex<double>* v){
    AllocSpace();    
    //We define the third vector orthogonal so as to avoid problems
    complex<double>* w = new complex<double>[3];
    w[0]=u[1]*v[2]-u[2]*v[1];
    w[1]=u[2]*v[0]-u[0]*v[2];
    w[2]=u[0]*v[1]-u[1]*v[0];

    complex<double>** V = new complex<double>*[3];
    for(size_t i = 0; i < 3; ++i){
        V[i] = new complex<double>[3];
    }
    for(size_t j = 0; j < 3; ++j){
        V[0][j]=u[j];
        V[1][j]=v[j];
        V[2][j]=w[j];
    }
    gramSchmidt(V,3);
    complex<double> deter = det(V,3);
    assert( abs(deter) !=0 );
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            U[i][j]=V[i][j]/pow(deter,1/3.);
        }
    }
    for(size_t i = 0; i < 3; ++i){
        delete [] V[i];
    }
    delete [] V;
    delete [] w;
    return;
}
SU3::SU3(complex<double>* u, complex<double>* v, complex<double>*w){
    AllocSpace();
    complex<double>** V = new complex<double>*[3];
    for(size_t i = 0; i < 3; ++i){
        V[i] = new complex<double>[3];
    }
    for(size_t j = 0; j < 3; ++j){
        V[0][j]=u[j];
        V[1][j]=v[j];
        V[2][j]=w[j];
    }
    gramSchmidt(V,3);
    complex<double> deter = det(V,3);
    assert( abs(deter) !=0 );
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            U[i][j]=V[i][j]/pow(deter,1/3.);
        }
    }
    for(size_t i = 0; i < 3; ++i){
        delete [] V[i];
    }
    delete [] V;
    return;
}
SU3::SU3(complex<double>** M, bool detNorm){
    AllocSpace();
    // det controls wether we want to normalise by the determinant
    complex<double> deter = det(M,3);
    if( abs(deter) != 0 && detNorm ){
        for(size_t i = 0; i < 3; ++i){
            for(size_t j = 0; j < 3; ++j){
                U[i][j]=M[i][j]/pow(deter,1/3.);
            }
        }
    }
    else{
        for(size_t i = 0; i < 3; ++i){
            for(size_t j = 0; j < 3; ++j){
                U[i][j]=M[i][j];
            }
        }
    }
    return;
}
SU3::SU3(double** M, bool detNorm){
    AllocSpace();
    //In case we want to initialize with reals only
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            U[i][j]=complex<double>(M[i][j],0.);
        }
    }
    complex<double> deter = det(U,3);
    if( abs(deter) != 0 && detNorm ){
        for(size_t i = 0; i < 3; ++i){
            for(size_t j = 0; j < 3; ++j){
                U[i][j]/=pow(deter,1/3.);
            }
        }
    }
    else{
        for(size_t i = 0; i < 3; ++i){
            for(size_t j = 0; j < 3; ++j){
                U[i][j]=complex<double>(M[i][j],0.);
            }
        }
    }
    return;
}
SU3::SU3(const SU3& M){
    AllocSpace();
    complex<double> deter = det(M.U,3);
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            U[i][j]=M.U[i][j];
        }
    }
    return;
}
SU3::~SU3(){
    for(size_t i = 0; i < 3; ++i){
        delete [] U[i];
    }
    delete [] U;
}

//---------------------------------------------------------------------------
//-----------------------Member functions------------------------------------
//---------------------------------------------------------------------------
SU3 SU3::Adj() const{
    //This one is NEVER passed by reference, it messes values up
    SU3 temp;
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            temp.U[i][j]=conj(U[j][i]);
        }
    }
    return (temp);
}
SU3 SU3::Inv() const{
    //This one is NEVER passed by reference, it messes values up
    //The inverse Inv routine is meant to be more accurate than Adj, since 
    //errors of order O(d) of the latter are turned into O(d^2)
    SU3 temp = *this;
    
    return (2.*temp.Adj() - temp.Adj()*temp*temp.Adj());
}
complex<double> SU3::DetSU3() const{
    return det(U,3);
}
SU3& SU3::RandSU3QR(){
    //Implementation of the QR algorithm with a single-valuedness correction
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> u(0.,1.), v(0.,1.);
    complex<double>** M = new complex<double>*[3];
    complex<double>** R = new complex<double>*[3];
    complex<double>** temp = new complex<double>*[3];
    double r, phi;
    
    for(size_t i = 0; i < 3; ++i){
        M[i] = new complex<double>[3];
        R[i] = new complex<double>[3];
        temp[i] = new complex<double>[3];
        for(size_t j = 0; j < 3; ++j){
            r=sqrt(2.*log(1./(u(gen)+0.0000000001)));
            phi = 2*3.141592653589793*v(gen);
            M[i][j]=temp[i][j]=complex<double>(r*cos(phi),r*sin(phi));
        }
    }
    //gramSchmidt(M,3);
    houseHolder(M,3u);
    //Since M=QR and Q^-1=Q^+, R=Q^+M
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            //We are only interested in \Lambda the diagonal unit entries of R
            if(i==j){
                complex<double> sum(0.,0.);
                for(size_t k = 0; k < 3; ++k){
                    sum += conj(M[k][i])*temp[k][j];
                }
                assert(abs(sum) != 0);
                R[i][j]=sum/abs(sum);
            }
            else{ R[i][j] = complex<double>(0.,0.);}
        }
    }
    //Now Q'=Q\Lambda is distributed with Haar measure (and numerically stable)
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            //We are only interested in \Lambda the diagonal unit entries of R
            
            complex<double> sum(0.,0.);
            for(size_t k = 0; k < 3; ++k){
                sum += M[i][k]*R[k][j];
            }
            temp[i][j]=sum;
        }
    }
    //Already determinant one
    SU3 Q(temp, false);
    Q /= pow(Q.DetSU3(),1./3.);
    for(size_t i = 0; i < 3; ++i){
        delete [] M[i];
        delete [] R[i];
        delete [] temp[i];
    } 
    return (*this = Q);
}
SU3& SU3::RandSU3(const double spread){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> xi(0.,1.);
    complex<double>* u = new complex<double>[3];
    complex<double>* v = new complex<double>[3];
    //We generate two random complex vectors
    u[0]=complex<double>(spread*(2.*xi(gen)-1.),spread*(2.*xi(gen)-1.));
    v[0]=complex<double>(spread*(2.*xi(gen)-1.),spread*(2.*xi(gen)-1.));
    u[1]=complex<double>(spread*(2.*xi(gen)-1.),spread*(2.*xi(gen)-1.));
    v[1]=complex<double>(spread*(2.*xi(gen)-1.),spread*(2.*xi(gen)-1.));
    u[2]=complex<double>(spread*(2.*xi(gen)-1.),spread*(2.*xi(gen)-1.));
    v[2]=complex<double>(spread*(2.*xi(gen)-1.),spread*(2.*xi(gen)-1.));
    
    SU3 temp(u,v);
    return (*this = temp);
}
SU3& SU3::SmallRandSU3(const double spread){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> xi(0.,1.);
    uniform_int_distribution<int> sgn(0,1);
    double* r = new double[3];
    double* s = new double[3];
    double* t = new double[3];
    //This has to be a double, not an int!!!!! Made me fail once
    double r0, s0, t0;
    r0=t0=s0=sqrt(1.-spread*spread);
    r0 *= 2*sgn(gen)-1;
    s0 *= 2*sgn(gen)-1;
    t0 *= 2*sgn(gen)-1;

    r[0]=2.*xi(gen)-1;
    r[1]=2.*xi(gen)-1;
    r[2]=2.*xi(gen)-1;
    if(!normalise(r,3u)){
        cout << "Error initializing SU(3) close to I.\n";
    }
    
    s[0]=2.*xi(gen)-1;
    s[1]=2.*xi(gen)-1;
    s[2]=2.*xi(gen)-1;
    if(!normalise(s,3u)){
        cout << "Error initializing SU(3) close to I.\n";
    }
    
    t[0]=2.*xi(gen)-1;
    t[1]=2.*xi(gen)-1;
    t[2]=2.*xi(gen)-1;
    if(!normalise(t,3u)){
        cout << "Error initializing SU(3) close to I.\n";
    }
    complex<double>** R = new complex<double>*[3];
    complex<double>** S = new complex<double>*[3];
    complex<double>** T = new complex<double>*[3];
    for(size_t i = 0; i < 3u; ++i){
        R[i] = new complex<double>[3];
        S[i] = new complex<double>[3];
        T[i] = new complex<double>[3];
        r[i] *= spread;
        s[i] *= spread;
        t[i] *= spread;        
    }
    R[0][0]=complex<double>(r0,-r[2]); R[1][1]=complex<double>(r0,r[2]); R[2][2]=complex<double>(1.,0.);
    R[0][1]=complex<double>(-r[1],-r[0]); R[1][0]=complex<double>(r[1],-r[0]); 
    R[0][2]=R[2][0]=R[1][2]=R[2][1]=complex<double>(0.,0.);
    S[0][0]=complex<double>(s0,-s[2]); S[1][1]=complex<double>(1.,0.); S[2][2]=complex<double>(s0,s[2]);
    S[0][2]=complex<double>(-s[1],-s[0]); S[2][0]=complex<double>(s[1],-s[0]); 
    S[0][1]=S[1][0]=S[1][2]=S[2][1]=complex<double>(0.,0.);
    T[1][1]=complex<double>(t0,-t[2]); T[2][2]=complex<double>(t0,t[2]); T[0][0]=complex<double>(1.,0.);
    T[1][2]=complex<double>(-t[1],-t[0]); T[2][1]=complex<double>(t[1],-t[0]); 
    T[0][1]=T[2][0]=T[1][0]=T[0][2]=complex<double>(0.,0.);

    //These already are SU(2) by construction
    SU3 RU(R, false), SU(S, false), TU(T, false);

    for(size_t i = 0; i < 3; ++i){
        delete [] R[i];
        delete [] S[i];
        delete [] T[i];
    }
    delete [] R; delete [] r;
    delete [] S; delete [] s;
    delete [] T; delete [] t;
    if( xi(gen) >= 0.5){
        return (*this = RU*SU*TU);
    }
    else return (*this = TU.Adj()*SU.Adj()*RU.Adj());
}
complex<double> SU3::TraceSU3() const{
    return U[0][0] + U[1][1] + U[2][2]; 
}
double SU3::NormSU3() const{
    double sum = 0;
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            sum += norm(U[i][j]);
        }
    }
    return sum;
}
double SU3::Dist(const SU3 V) const{
    double sum = 0;
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            sum += abs( U[i][j] - V.U[i][j] )*abs( U[i][j] - V.U[i][j] );
        }
    }
    return sum;
}
double SU3::AbsSU3() const{
    SU3 temp = *this;
    return sqrt(temp.NormSU3());
}
void SU3::AssertDeterminant(){
    SU3 temp = *this;
    temp /= pow(temp.DetSU3(),1./3.);
    *this = temp; 
}
void SU3::PrintSU3() const{
    cout << setprecision(5);
    for(size_t i = 0; i < 3; ++i){
        cout << "( " << setw(3);
        for(size_t j = 0; j < 3; ++j){
            string sign = " + ";
            if(U[i][j].imag() == 0){
                cout << U[i][j].real() << setw(15);
            }
            else if(U[i][j].imag() < 0){
                sign = " - ";
                if(U[i][j].real() == 0){
                    cout << U[i][j].imag() << "i" << setw(17);
                }
                else{
                    cout << U[i][j].real() << sign << abs(U[i][j].imag()) << "i" << setw(17);    
                }
            }
            else if(U[i][j].imag() > 0){
                if(U[i][j].real() == 0){
                    cout << U[i][j].imag() << "i" << setw(17);
                }
                else{
                    cout << U[i][j].real() << sign << abs(U[i][j].imag()) << "i" << setw(17);    
                }
            }
        }
        cout << " )\n";
    }
    return;
}
SU3 SU3::ConstructIdentity(){
    SU3 temp;
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            if(i==j){
                temp.U[i][j]=complex<double>(1.,0.);
            }
            else temp.U[i][j]=complex<double>(0.,0.);
        }
    }
    return temp;
}
//---------------------------------------------------------------------------
//-------------------------Operator overloads--------------------------------
//---------------------------------------------------------------------------

SU3& SU3::operator=(const SU3& V){
    if(this == &V){
        return *this;
    }
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            U[i][j]=V.U[i][j];
        }
    }
    return *this;
}
SU3& SU3::operator+=(const SU3& V){
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            U[i][j] += V.U[i][j];
        }
    }
    return *this;
}
SU3& SU3::operator-=(const SU3& V){
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            U[i][j] -= V.U[i][j];
        }
    }
    return *this;
}
SU3& SU3::operator*=(const SU3& V){
    SU3 temp;

    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            temp.U[i][j] = complex<double>(0.,0.);
            for(size_t k = 0; k < 3; ++k){
                temp.U[i][j] += (U[i][k]*V.U[k][j]);
            }
        }
    }
    return (*this = temp);
}
SU3& SU3::operator*=(complex<double> a){
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            U[i][j] *= a;
        }
    }
    return *this;
}
SU3& SU3::operator*=(double a){
    return (*this *= complex<double>(a,0.));
}
SU3& SU3::operator/=(complex<double> a){
    assert( norm(a)!=0 );
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            U[i][j] /= a;
        }
    }
    return *this;
}
SU3& SU3::operator/=(double a){
    assert( a!=0 );
    return (*this /= complex<double>(a,0.));
}
SU3 SU3::operator^(int n){
    assert( n >= 0 );
    if( n == 0 ){
        //Constructs the identity
        SU3 M(1.);
        return M;
    }
    else if( n == 1 ){
        return *this;
    }
    else{
        SU3 temp = *this;
        for(size_t m = 1; m < n; m++){
           temp *= temp;
        }
        return temp;
    }
}
//---------------------------------------------------------------------
//---------Non-member binary operators (+,-,*,/)-----------------------
//---------------------------------------------------------------------

SU3 operator+(const SU3& U1, const SU3& U2){
    SU3 temp = U1;
    return (temp += U2);
}
SU3 operator-(const SU3& U1, const SU3& U2){
    SU3 temp = U1;
    return (temp -= U2);
}
SU3 operator*(const SU3& U1, const SU3& U2){
    SU3 temp = U1;
    return (temp *= U2);
}
SU3 operator*(const SU3& U1, complex<double> a){
    SU3 temp = U1;
    return (temp *= a);
}
SU3 operator*(complex<double> a, const SU3& U1){
    SU3 temp(U1);
    return (temp *= a);
}
SU3 operator*(const SU3& U1, double a){
    SU3 temp(U1);
    return (temp *= complex<double>(a,0.));
}
SU3 operator*(double a, const SU3& U1){
    SU3 temp(U1);
    return (temp *= complex<double>(a,0.));
}
SU3 operator/(const SU3& U1, complex<double> a){
    SU3 temp(U1);
    assert( norm(a)!=0 );
    return (temp /= a);
}
SU3 operator/(const SU3& U1, double a){
    SU3 temp(U1);
    assert( a != 0);
    return (temp /= complex<double>(a,0.));
}
//---------------------------------------------------------------------
//--------------Non-member functions (miscellaneous)-------------------
//---------------------------------------------------------------------

double Dist(SU3 U, const SU3 V){
    return U.Dist(V);
}