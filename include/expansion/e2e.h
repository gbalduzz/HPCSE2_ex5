#pragma once
#include <unordered_map>
#include <complex>
#include <assert.h>
int computeBinomial(int a,int b);

int binomial(int a, int b){
  const static int MAX_B=50;
  assert(b<=MAX_B);
  static std::unordered_map<int,int> known;
  int* result = &known[MAX_B*a+b];
  if(not *result) *result= computeBinomial(a,b);
  return *result;
}

int computeBinomial(int a,int b){
  if(b==0 or a==b) return 1;
  else return binomial(a-1,b)+binomial(a-1,b-1);
}


template<int ord>
void e2e(const double* c_re,const double* c_im, double z_re,double z_im,double* new_re,double* new_im){
  const std::complex<double> z(z_re,z_im);
  new_re[0]+=c_re[0];
  for(int l=1;l<ord;l++) {
    new_re[l]+=c_re[0]*(std::pow(z,l)).real()/l;
    for(int k=1;k<l;k++){
      const std::complex<double> a(c_re[k],c_im[k]);
      new_re[l]+=(a*std::pow(z,l-k)).real()*binomial(l-1,k-1);
      new_im[l]+=(a*std::pow(z,l-k)).imag()*binomial(l-1,k-1);
    }
  }
}
