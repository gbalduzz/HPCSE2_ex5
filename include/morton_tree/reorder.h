//
// Created by giovanni on 08.03.16.
//
#ifndef HPCSE2_EX2_REORDER_H
#define HPCSE2_EX2_REORDER_H

void reorder(const int N, const uint* const keys, const float * const x,const float * const y,const float* const m,const uint* const label,float *xsorted,float *ysorted,float* msorted)
{
  #pragma omp parallel for schedule(static)
    for(int i=0;i<N;i++){
        xsorted[i]=x[keys[i]];
        ysorted[i]=y[keys[i]];
        msorted[i]=m[keys[i]];
    }
}

#endif //HPCSE2_EX2_REORDER_H
