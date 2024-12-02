#include <cmath>
#include <complex>
#include <iostream>
#include <assert.h>

#include "utils.hh"

using namespace std;

bool gramSchmidt(complex<double>** v, size_t n){
    if(norm(det(v,n)) == 0){   
        cout << "Warning: zero determinant. Dependent system!\n";
        return false;
    }
    complex<double>** temp = new complex<double>*[n];
    for(size_t i = 0; i < n; i++){
        temp[i] = new complex<double>[n];
    }
    for(size_t i = 0; i < n; i++){
        for(size_t j = 0; j < n; j++){
            temp[i][j]=v[i][j];
            for(size_t k = 0; k < i; k++){
                temp[i][j]-=scalarProd(v[k],v[i],n)*v[k][j];
            }
        }
        normalise(temp[i],n);
        copyVector(temp[i],v[i],n);
    }
    for(size_t i = 0; i < n; i++){
        delete [] temp[i];
    }
    delete [] temp;
    return true;
}
bool houseHolder(complex<double>** v, size_t n){
    if(norm(det(v,n)) == 0){
        cout << "Warning: zero determinant. Rank < " << n << " matrix\n";
        return false;
    }
    //QR factorisation routine which returns the unitary matrix Q in the original pointer v
    //More stable numerically than Gram-Schmidt, but it does the same orthogonalization
    complex<double>** A = new complex<double>*[n];
    complex<double>** Q = new complex<double>*[n];
    complex<double>** tempQ = new complex<double>*[n];
    for(size_t i = 0; i < n; i++){
        A[i] = new complex<double>[n];
        Q[i] = new complex<double>[n];
        tempQ[i] = new complex<double>[n];
        for(size_t j = 0; j < n; j++){
            A[i][j] = v[i][j];
            if(j == i){
                Q[i][j] = complex<double>(1.,0.);
            }
            else Q[i][j] = complex<double>(0.,0.);
        }
    }
    for(size_t m = 0; m < n-1u; m++){
        complex<double>* u = new complex<double>[n-m];
        //Matrix A^(m) is going to be the m-th minor operation at A^(m-1)_00
        //This is, a matrix of one less row (0) and column (0) every iteration
        //In practice we just start at a higher index m for both i,j. Final
        //multiplications are carried out the usual way for the full nxn matrices.
        double sum = 0;
        for(size_t i = 0; i < n; i++){
            for(size_t j = 0; j < n; j++){
                if( j == i ){
                    tempQ[i][j] = complex<double>(1.,0.);
                }
                else tempQ[i][j] = complex<double>(0.,0.);
            }
        }
        for(size_t i = m; i < n; i++){
            sum += norm(A[i][0]);
        }
        //We construct n-m dim vector u = x - alpha e_1 where x is the first column of A starting at pivot m
        u[0] = A[0][0] + complex<double>(cos(arg(A[0][0])),sin(arg(A[0][0])))*sqrt(sum);
        for(size_t i = m+1; i < n; i++){
            u[i-m] = A[i-m][0]; 
        }
        normalise(u,n-m);
        //We fill in matrix Q_m
        for(size_t i = m; i < n; i++){
            for(size_t j = m; j < n; j++){
                tempQ[i][j] -= 2.*u[i]*conj(u[j]);
            }
        }
        multiplyMatrix(tempQ,A,A,n);
        //Now Q_mA is upper diagonal in the first m columns
        // Q = Q1 ... Qm
        multiplyMatrix(Q,tempQ,Q,n);
        delete [] u;
    }
    for(size_t i = 0; i < n; i++){
        for(size_t j = 0; j < n; j++){
            v[i][j] = Q[i][j];
        }
        delete [] A[i];
        delete [] Q[i];
        delete [] tempQ[i];
    }

    delete [] A;
    delete [] Q;
    delete [] tempQ;
    return true;
}
void multiplyMatrix(complex<double>** m1, complex<double>** m2, complex<double>** mOut, size_t n){
    complex<double>** temp = new complex<double>*[n];
    for(size_t i = 0; i < n; ++i){
        temp[i] = new complex<double>[n];
        for(size_t j = 0; j < n; ++j){
            temp[i][j] = complex<double>(0.,0.);
            for(size_t k = 0; k < n; ++k){
                temp[i][j] += m1[i][k]*m2[k][j]; 
            }
            //assert(abs(mOut[i][j]) != 0);
        }
    }
    for(size_t i = 0; i < n; ++i){
        for(size_t j = 0; j < n; ++j){
            mOut[i][j] = temp[i][j];
        }
        delete [] temp[i];
    }
    delete [] temp;
    return;
}
void copyVector(complex<double>* in, complex<double>* out, size_t n){
    for(size_t i = 0; i < n; i++){
        out[i]=in[i];
    }
    return;
}
void printVector(complex<double>* v, size_t n){
    cout << "(";
    for(size_t i = 0; i < n-1; i++){
        if(v[i].imag() < 0){
            cout << v[i].real() << "-" << abs(v[i].imag()) << "i, ";
        }
        else if(v[i].imag() > 0){
            cout << v[i].real() << "+" << v[i].imag() << "i, ";
        }
        else{
            cout << v[i].real() << ", ";
        }
    }
    if(v[n-1].imag() < 0){
        cout << v[n-1].real() << "-" << abs(v[n-1].imag()) << "i )";
    }
    else if(v[n-1].imag() > 0){
        cout << v[n-1].real() << "+" << v[n-1].imag() << "i )";
    }
    else{
        cout << v[n-1].real() << " )";
    } 
}
bool normalise(complex<double>* v, size_t n){
    double sum = 0.;
    for(size_t i = 0; i < n; ++i){
        sum += norm(v[i]);
    }
    if(sum == 0){
        return false;
    }
    else for(size_t i = 0; i < n; ++i){
        v[i]/=sqrt(sum);
    }
    return true;
}
bool normalise(double* v, size_t n){
    double sum = 0.;
    for(size_t i = 0; i < n; ++i){
        sum += v[i]*v[i];
    }
    if(sum == 0){
        return false;
    }
    else for(size_t i = 0; i < n; ++i){
        v[i]/=sqrt(sum);
    }
    return true;
}
complex<double> scalarProd(complex<double>* u, complex<double>* v, size_t n){
    //Sesquilinear, so order is clearly important
    complex<double> sum = complex<double>(0.,0.);
    for(size_t i = 0; i < n; i++){
        sum += conj(u[i])*v[i];
    }
    return sum;
}
void minor(complex<double>** M, complex<double>** m, size_t i, size_t j, size_t ncols){
    assert(j >= 0 && j < ncols);
    assert(i >= 0 && i < ncols);
    for(size_t I = 0; I < ncols; I++){
        for(size_t J = 0; J < ncols; J++){
            if(I < i){
                if(J < j){
                    m[I][J]=M[I][J];
                }
                else if(J > j){
                    m[I][J-1]=M[I][J];
                }
            }
            else if(I > i){
                if(J < j){
                    m[I-1][J]=M[I][J];
                }
                else if(J > j){
                    m[I-1][J-1]=M[I][J];
                }
            }
        }
    }
    return;
}
complex<double> det(complex<double>** m, size_t ncols){
    //We take a minor decomposition starting at the given order
    size_t order = ncols;
    if (ncols==2u){
        return m[0][0]*m[1][1]-m[1][0]*m[0][1];
    }
    else if(ncols == 1u){
        return m[0][0];
    }
    /*else if(ncols == 3u){
        return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) 
            - m[0][1]*(m[1][0]*m[2][2]-m[2][0]*m[1][2]) 
            + m[0][2]*(m[1][0]*m[2][1]-m[2][0]*m[1][1]); 
    }*/
    else{
        complex<double> sum; 
        sum = complex<double>(0.,0.);
        int factor = 1;
        complex<double>** temp = new complex<double>*[ncols-1];
        for(size_t k = 0; k < ncols-1; k++){
            temp[k] = new complex<double>[ncols-1];
        }
        for(size_t j = 0; j < ncols; j++){  
            //We compute the determinant of the (0,j) th minor (first row)  
            minor(m,temp,0,j,ncols);
            sum += ((double)factor)*m[0][j]*det(temp, ncols-1);
            factor *= -1;
        }
        for(size_t k = 0; k < ncols-1; k++){
            delete [] temp[k];
        }
        delete [] temp;
        return sum;
    } 
}
size_t I(int* i, size_t Ns, size_t dim){
    size_t temp = 0;
    size_t factor = 1;
    for(size_t d = 0; d < dim; d++){
        temp += factor*i[d];
        factor *= Ns;
    }
    return temp;
}
void ijkn(int* i, size_t Ns, size_t I, size_t dim){
    size_t factor = 1;
    
    for(size_t d = 0; d < dim-1; d++){ factor *= Ns; }
    for(size_t d = dim-1; d > 0; --d){
        
        i[d] = static_cast<int>(I/factor);
        I = I%factor;
        factor /= Ns;
    }
    i[0] = static_cast<int>(I);
    return;
}
double min(const double a1, const double a2){
    if(a1 <= a2){
        return a1;
    }
    else return a2;
}
int max(const int a1, const int a2){
    if(a1 >= a2){
        return a1;
    }
    else return a2;
}
double mean(double* v, const size_t N){
    double sum = 0;
    for(size_t n = 0; n < N; n++){
        sum += v[n];
    }
    return sum/N;
}
double variance(double* v, double meanV, const size_t N){
    double sum = 0;
    for(size_t n = 0; n < N; n++){
        sum += (v[n]-meanV)*(v[n]-meanV);
    }
    return sum/N;
}