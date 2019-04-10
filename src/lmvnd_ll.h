#ifndef LMVND_LL    // To make sure you don't declare the function more than once by including the header multiple times.
#define LMVND_LL

double lmvnd(std::vector<double> Sigma, std::vector<double> X, std::vector<double> mean, int imax );

#endif

#ifndef LGD    // To make sure you don't declare the function more than once by including the header multiple times.
#define LGD

double lgammad(double x, double a, double b);

#endif