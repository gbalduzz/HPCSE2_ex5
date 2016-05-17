#ifndef PARTICLES_H
#define PARTICLES_H
#include <vector>
#include <complex>

struct Particles{
  int N; //passed from global scope
  double  *x,*y,*w;

  Particles(){N=0;}
  Particles(int N0);
  ~Particles();
  Particles subEnsamble(int i0,int l);
  void resize(int N);
};

Particles::Particles(int N0):N(N0)
{
  x = new double[N];
  y = new double[N];
  w = new double[N];
}

Particles::~Particles()
{
  if(N) {
    delete[] x;
    delete[] y;
    delete[] w;
  }
}

Particles Particles::subEnsamble(int i0, int l) {
  assert(i0+l <= N);
  Particles p2;
  p2.x = x+i0;
  p2.y = y+i0;
  p2.w = w+i0;
  p2.N = l;
  return p2;
}

void Particles::resize(int N0) {
  if(N == N0) return;
  if(N){
    delete[] x;
    delete[] y;
    delete[] w;
  }
  N=N0;
  x = new double[N];
  y = new double[N];
  w = new double[N];
}

#endif
