#ifndef PARTICLES_H
#define PARTICLES_H
#include <vector>
#include <complex>

struct Particles{
  double  *x,*y,*w;
  int N;
  bool is_a_copy;

  Particles(){is_a_copy=true;}
  Particles(int N0);
  Particles(double* ,double* ,double* ,const int );
  ~Particles();
  Particles subEnsamble(int i0,int l) const;
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
  Particles p2;
  p2.x = x+i0;
  p2.y = y+i0;
  p2.w = w+i0;
  p2.N = l;
  return p2;
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
