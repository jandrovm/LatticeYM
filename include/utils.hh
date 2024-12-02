#ifndef UTILS_HH
#define UTILS_HH

bool gramSchmidt(std::complex<double>** v, size_t n);
bool houseHolder(std::complex<double>** v, size_t n);
void copyVector(std::complex<double>* in, std::complex<double>* out, size_t n);
void multiplyMatrix(std::complex<double>** m1, std::complex<double>** m2, std::complex<double>** mOut, size_t n);
bool normalise(std::complex<double>* v, size_t n);
bool normalise(double* v, size_t n);
double min(const double a1, const double a2);
double mean(double*, const size_t);
double variance(double*, double meanV, const size_t);
int max(const int a1, const int a2);
std::complex<double> scalarProd(std::complex<double>* u, std::complex<double>* v, size_t n);
std::complex<double> det(std::complex<double>** m, size_t ncols);
void minor(std::complex<double>** M, std::complex<double>** m, size_t i, size_t j, size_t ncols);
void printVector(std::complex<double>* v, size_t n);
size_t I(int* i, size_t Ns, size_t dim);
void ijkn(int* i, size_t Ns, size_t I, size_t dim);

#endif