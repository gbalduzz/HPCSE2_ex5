#pragma once
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
  friend void swap(Particles& a,Particles& b);
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
  return Particles(x+i0,y+i0,w+i0,l);
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

void swap(Particles& a,Particles& b){
  assert(not (a.is_a_copy and b.is_a_copy));
  double *xtmp(a.x),*ytmp(a.y),*wtmp(a.w);
  int ntmp=a.N;
  a.x=b.x; a.y=b.y; a.w=b.w;
  a.N=b.N;
  b.x=xtmp; b.y=ytmp; b.w=wtmp;
  b.N=ntmp;
}

#include<iostream>
using std::cout; using std::endl;
void Print(const Particles& p,int n=0){
  n= n? std::min(n,p.N) : p.N;
  cout<<"x\ty\tw\n";
  for(int i=0;i<n;i++) cout<<p.x[i]<<"\t"<<p.y[i]<<"\t"<<p.w[i]<<endl;
}
