#pragma once
#include <omp.h>
#include <cmath>
#include "helper_methods.h"

double p2p( const Particles& p,
            double xtarget, double ytarget)
{
  double res(0);
#pragma omp simd reduction(+:res)
//#pragma ivdep
  for(int i=0;i<p.N;i++){
   res+=p.w[i]*std::log(squareDistance(p.x[i],xtarget,p.y[i],ytarget));
 }
  return res/2.;
}
