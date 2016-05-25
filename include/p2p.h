#pragma once
#include "p2p.ispc.h"

inline double p2p( const Particles& p,
            double xtarget, double ytarget)
{
  return ispc::p2p(p.x,p.y,p.w,xtarget,ytarget,p.N);
}
