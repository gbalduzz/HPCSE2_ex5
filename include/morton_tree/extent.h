//
// Created by giovanni on 07.03.16.
//

#ifndef EX2_EXTENT_H
#define EX2_EXTENT_H
#include<omp.h>
#include<vector>
using pairf=std::pair<float,float>;

pairf serial_MinMax(const float * const  x,int N);

pairf MinMax(const float * const  x,int N)
//parallel implementation
{
const int n_threads=omp_get_max_threads();
  std::vector<pairf> block_minmax(n_threads);
  const int block_size=N/n_threads;
  const int last_block=block_size+N%n_threads;
#pragma omp parallel
  {
      int t=omp_get_thread_num();
      block_minmax[t]=serial_MinMax(x+t*block_size, (t==n_threads-1) ? block_size : last_block);
  }
  double min=block_minmax[0].first;
  double max=block_minmax[0].second;
  for(int t=0;t<n_threads;t++){
    if(block_minmax[t].first<min) min=block_minmax[t].first;
    if(block_minmax[t].second>max) max=block_minmax[t].second;
  }
  return pairf(min,max);
}

void extent(const int N, const float * const x,const float * const y,float& xmin,float& ymin,float& ext)
{
pairf x_mm=MinMax(x,N);
pairf y_mm=MinMax(y,N);
xmin=x_mm.first;
ymin=y_mm.first;
float  dx(x_mm.second-x_mm.first);
float  dy(y_mm.second-y_mm.first);
ext= dx>dy ? dx : dy;
}

pairf serial_MinMax(const float * const  x,int N)
{
    double min(x[0]),max(x[0]);
    for(int i=1;i<N;i++) {
        if (x[i] > max) max = x[i];
        else if (x[i] < min) min = x[i];
    }
    return pairf(min,max);
}

#endif //EX2_EXTENT_H
