#ifndef PARTICLES_H
#define PARTICLES_H
#include <vector>
#include <complex>

using Complex = std::complex<double>;

struct Particles{
    const int N; //passed from global scope
    std::vector<double> w,x,y;
    Particles(int N):x(N),y(N),w(N),N(N){}
    };

#endif