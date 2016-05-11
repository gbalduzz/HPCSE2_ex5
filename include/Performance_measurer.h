#pragma once
#include"profiler.h"
#include <cmath>

template<typename F,typename... Args>
std::pair<double,double> Performance(const int n_meas,F& f,Args&... args)
// return average and std error of time measurements of function F with arguments args
{
    double mean,mean2,meas;
    for(int i=0;i<n_meas;i++){
        Profiler pr("");
        f(args...);
        meas=pr.stop(false);
        mean+=meas;
        mean2+=meas*meas;
    }
    mean/=n_meas; mean2/=n_meas;
    return std::pair<double,double>(mean,(mean2-mean*mean)/std::sqrt(n_meas));
}