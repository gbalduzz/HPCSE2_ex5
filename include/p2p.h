#pragma once
#include <cmath>

inline double squareDistance(double x1,double x2,double y1,double y2){
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}

double p2p( const Particles& p,
            double xtarget, double ytarget)
{
  double res(0);
 for(int i=0;i<p.N;i++){
   res+=p.w[i]*std::log(squareDistance(p.x[i],xtarget,p.y[i],ytarget));
 }
  return res/2.;
}
