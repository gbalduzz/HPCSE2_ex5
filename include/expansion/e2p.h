#pragma once

#include <omp.h>
#pragma omp declare simd
template<int K>
double e2p(const double z_re,const double z_im,const double* c_re,const double* c_im){
  //c_re[0]=Q
  double result= c_re[0]*0.5*std::log(z_re*z_re+z_im*z_im);
  
  double zk_re=1,zk_im=0;
#pragma unroll
  for(int k=1;k<K+1;k++) {
    const double temp=z_re*zk_re-z_im*zk_im;
    zk_im=z_re*zk_im+z_im*zk_re;
    zk_re=temp;
    result += (c_re[k]*zk_re+c_im[k]*zk_im)/(zk_re*zk_re+zk_im*zk_im);
  }
  
  return result;
}


//buffered version
template<int ord>
double buffered_e2p(const double* rzs,const double* izs,
                         const double** rxps,const double** ixps,const int bufcount){
  double result=0;
#pragma omp simd reduction(+:result)
  for(int i=0;i<bufcount;i++) result += e2p<ord>(rzs[i],izs[i],rxps[i],ixps[i]);
  return result;
}
