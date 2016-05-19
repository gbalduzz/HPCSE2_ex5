#ifndef PARTICLES_H
#define PARTICLES_H
#include <vector>
#include <complex>
#include <assert.h>

struct Particles{
  double  *x,*y,*w;
  int N;
  bool is_a_copy;

  Particles(){is_a_copy=true;}
  Particles(int N0);
  inline Particles(double* ,double* ,double* ,const int );
  ~Particles();
  inline Particles subEnsamble(int i0,int l) const;
  void resize(int N);
};

Particles::Particles(int N0):N(N0)
{
  is_a_copy=false;
  x = new double[N];
  y = new double[N];
  w = new double[N];
}

Particles::~Particles()
{
  if(not is_a_copy) {
    delete[] x;
    delete[] y;
    delete[] w;
  }
}

Particles::Particles(double* x0,double* y0,double* w0,const int N0):
x(x0),y(y0),w(w0),N(N0),is_a_copy(true)
{}

Particles Particles::subEnsamble(int i0, int l) const{
  assert(i0+l <= N);
  return Particles(x+i0,y+i0,w+i0,N);;
}

void Particles::resize(int N0) {
  N=N0;
  if(not is_a_copy) {
    delete[] x;
    delete[] y;
    delete[] w;
  }
  is_a_copy=false;
  x = new double[N];
  y = new double[N];
  w = new double[N];
}

#endif
