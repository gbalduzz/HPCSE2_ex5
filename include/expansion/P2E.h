//
// Created by giovanni on 25.04.16.
//
#pragma once
#include "particles.h"
#include <complex>
#include <cmath>

template<int k,int K>
struct sum_kth_coeff{
        inline static void execute (double* c_re,double* c_im,
                                    const double x,const double y,const double w,const double z_re,const double z_im){
          c_re[k]-=w*z_re;
          c_im[k]-=w*z_im;            //compute z=(x+i y)^k
          const double new_re=z_re*x-z_im*y;
          const double new_im=z_re*y+z_im*x;
          sum_kth_coeff<k+1,K>::execute(c_re,c_im,x,y,w,new_re,new_im);
        }
};

template<int k>
struct sum_kth_coeff<k,k>{
    inline static void execute (double* c_re,double* c_im,
                                const double x,const double y,const double w,const double z_re,const double z_im){
        c_re[k]-=w*(z_re*x-z_im*y);
        c_im[k]-=w*(z_re*y+z_im*x);
    }//end recursion
};

template <int k>
void P2E(const Particles& p,double xcom,double ycom,double* c_re,double* c_im){ //c_re and c_im must have k+1 reserved spaces
       for(int i=0;i<p.N;i++){
         const double x=p.x[i]-xcom;
         const double y=p.y[i]-ycom;
         const double w=p.w[i];
         c_re[0]+=w;
         sum_kth_coeff<1,k>::execute(c_re,c_im,x,y,w,x,y);
    }
    //divide by coeff order
#pragma omp simd
    for(int i=2;i<k+1;i++) {c_re[i]/=i;c_im[i]/=i;}

}
