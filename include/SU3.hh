#ifndef SU3_HH
#define SU3_HH

#include <complex>
#include <assert.h>

class SU3{
    public:
        SU3();
        ~SU3();
        SU3(const double a);
        SU3(std::complex<double>* u, std::complex<double>* v);
        SU3(std::complex<double>* u, std::complex<double>* v, std::complex<double>* w);
        SU3(std::complex<double>** , bool);
        SU3(double** , bool);
        SU3(const SU3&);
        
        SU3 Adj() const;
        SU3 Inv() const;
        SU3& RandSU3(const double spread);
        SU3& RandSU3QR();
        SU3& SmallRandSU3(const double spread);
        SU3 ConstructIdentity();
        double Dist(const SU3) const;
        std::complex<double> TraceSU3() const;
        //Just like complex numbers these are mag. squared, and abs. value, respectively.
        double NormSU3() const;
        double AbsSU3() const;
        std::complex<double> DetSU3() const;
        void PrintSU3() const;
        void AssertDeterminant();

        SU3& operator=(const SU3&);
        inline std::complex<double>& operator()(int x,int y) { assert(x < 3 && y < 3 && x >= 0 && y >= 0); return U[x][y]; }

        SU3& operator+=(const SU3&);
        SU3& operator-=(const SU3&);
        SU3& operator*=(const SU3&);
        SU3& operator*=(std::complex<double>);
        SU3& operator*=(double);
        SU3& operator/=(std::complex<double>);
        SU3& operator/=(double);
        SU3 operator^(int);

    private:
        //We may need to toggle between group and algebra operations 
        std::complex<double>** U;
        void AllocSpace();
};

SU3 operator+(const SU3&, const SU3&);
SU3 operator-(const SU3&, const SU3&);
SU3 operator*(const SU3&, const SU3&);
SU3 operator*(const SU3&, std::complex<double>);
SU3 operator*(std::complex<double>, const SU3&);
SU3 operator*(const SU3&, double);
SU3 operator*(double, const SU3&);
SU3 operator/(const SU3&, double);
SU3 operator/(const SU3&, std::complex<double>);

double Dist(const SU3&, const SU3&);
#endif
