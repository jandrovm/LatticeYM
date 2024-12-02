#include <iostream>
#include <cmath>
#include <assert.h>
#include "SU3.hh"
#include "utils.hh"

#define I complex<double>(0.,1.)
using namespace std;

int main(void){
    bool test_orto = true;
    bool test_unitary = false;
    bool test_SU3_operations = false;
    if(test_orto){
        complex<double>* u = new complex<double>[3]; 
        complex<double>* v = new complex<double>[3]; 
        complex<double>* w = new complex<double>[3];

        u[0] = 1.5*I+0.3; u[1] = complex<double>(1.,0.24); u[2] = complex<double>(-0.1314,3.5395);
        v[0] = complex<double>(4.,0.); v[1] = 3.2*I+0.12; v[2] = complex<double>(-0.2,-0.57);
        w[0] = complex<double>(0.235,-12.245); w[1] = complex<double>(1.,0.); w[2] = complex<double>(-2.,-1124124.44);

        cout << "Initial vectors: \n";
        cout << "u = ";
        cout << u[0].real() << "+" << u[0].imag() << "i" << "\n";
        printVector(u,3);
        cout << ", ";

        cout << "v = ";
        printVector(v,3);
        cout << ", ";

        cout << "w = ";
        printVector(w,3);
        cout << "\n";

        complex<double>** M = new complex<double>*[3];
        for(int i = 0; i < 3; i++){
            M[i] = new complex<double>[3];
        }
        copyVector(u,M[0],3);
        copyVector(v,M[1],3);
        copyVector(w,M[2],3);

        gramSchmidt(M,3);

        cout << "Orthonormalized (Gram-Schmidt) vectors: \n";
        cout << "u' = ";
        printVector(M[0],3);
        cout << ", ";

        cout << "v' = ";
        printVector(M[1],3);
        cout << ", ";

        cout << "w' = ";
        printVector(M[2],3);
        cout << "\n";

        copyVector(u,M[0],3);
        copyVector(v,M[1],3);
        copyVector(w,M[2],3);

        houseHolder(M,3);

        cout << "Orthonormalized (Householder) vectors: \n";
        cout << "u' = ";
        printVector(M[0],3);
        cout << ", ";

        cout << "v' = ";
        printVector(M[1],3);
        cout << ", ";

        cout << "w' = ";
        printVector(M[2],3);
        cout << "\n";

        cout << "Scalar products: \n";
        cout << "< u', v' > = " << scalarProd(M[0],M[1],3) << "\n";
        cout << "< u', w' > = " << scalarProd(M[0],M[2],3) << "\n";
        cout << "< v', w' > = " << scalarProd(M[1],M[2],3) << "\n";

        //Delete all pointers
        delete [] u;
        delete [] v;
        delete [] w;
        for(int i = 0; i < 3; i++){
            delete [] M[i];
        }
        delete [] M;
    }
    if(test_unitary){

        SU3 Id;
        Id.ConstructIdentity();
        SU3 R1, R2, R3;
        SU3 Id1, Id2, Id3;
        
        R1.RandSU3(1.);
        R2.RandSU3(1.);
        R3.RandSU3(1.);

        Id1 = R1*R1.Adj();
        Id2 = R2*R2.Adj();
        Id3 = R3*R3.Adj();
        
        double d = Id1.Dist(Id);
        complex<double> determ = R1.DetSU3();
        cout << "||R_1R_1^+-I_3||^2 = " << d << ", det(R_1) = " << determ << "\n";
        d = Id3.Dist(Id);
        determ = R2.DetSU3();
        cout << "||R_2R_2^+-I_3||^2 = " << d << ", det(R_2) = " << determ << "\n";
        d = Id3.Dist(Id);
        determ = R3.DetSU3();
        cout << "||R_3R_3^+-I_3||^2 = " << d << ", det(R_3) = " << determ << "\n";

        cout << "R_1 = \n";
        R1.PrintSU3();
        cout << ",\n";
        cout << "R_2 = \n";
        R2.PrintSU3();
        cout << ",\n";
        cout << "R_3 = \n";
        R3.PrintSU3();
        cout << ".\n\n";

        cout << "R_1R_1^+ = \n";
        Id1.PrintSU3();
        cout << ",\n";
        cout << "R_2R_2^+ = \n";
        Id2.PrintSU3();
        cout << ",\n";
        cout << "R_3R_3^+ = \n";
        Id3.PrintSU3();
        cout << ".\n\n";
    
    }
    if(test_SU3_operations){
        complex<double>** M = new complex<double>*[3];
        for(int i = 0; i < 3; i++){
            M[i] = new complex<double>[3];
        }
        M[0][0]=M[1][1]=complex<double>(0.,0.);
        M[2][2]=complex<double>(1.,1.);
        M[1][0]=I; M[0][1]=-I;
        M[2][0]=-M[0][2]=complex<double>(0.,0.);
        M[2][1]=-M[1][2]=complex<double>(0.,0.);
        
        SU3 A(M,true);
        //(A^1).printSU3();
        cout << "\nMatrix A:\n";
        A = A;
        cout << "A^0 = \n"; (A^0).PrintSU3(); cout << ",\n";
        cout << "A^1 = \n"; (A^1).PrintSU3(); cout << ",\n";
        cout << "A^2 = \n"; (A^2).PrintSU3(); cout << ",\n";
        cout << "A^3 = \n"; (A^3).PrintSU3(); cout << ",\n";
        for(int i = 0; i < 3; i++){
            delete [] M[i];
        }
        delete [] M;
        cout << "Small random SU(3): \n";
        SU3 X, Y, Z;
        X.SmallRandSU3(0.1);
        Y.SmallRandSU3(0.2);
        Z.SmallRandSU3(0.01);
        cout << "eps = 0.1. X = \n";
        X.PrintSU3();
        cout << ",\n";
        cout << "eps = 0.2. Y = \n";
        Y.PrintSU3();
        cout << ",\n";
        cout << "eps = 0.01. Z = \n";
        Z.PrintSU3();
        cout << ",\n";
        cout << "X + Y + Z = \n";
        (X+Z+Y).PrintSU3();
        cout << ",\n";
        cout << "X - Y - Z = \n";
        (X-Y-Z).PrintSU3();
        cout << ",\n";

        cout << "A^3(X - Y - Z)+2X = \n";
        ((A^3)*(X-Y-Z)+2*X).PrintSU3();
        cout << ",\n";

        cout << "Uniform random SU(3): \n";
        SU3 X1, Y1, Z1;
        X1.RandSU3QR();
        Y1.RandSU3QR();
        Z1.RandSU3QR();
        cout << "A = \n";
        X1.PrintSU3();
        cout << ",\n";
        cout << "B = \n";
        Y1.PrintSU3();
        cout << ",\n";
        cout << "C = \n";
        Z1.PrintSU3();
        cout << ",\n";
    }    
    return 0;
}
