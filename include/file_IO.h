#pragma once
#include <complex>
#include <string>
#include <cmath>
#include <fstream>
#include <assert.h>
#include "particles.h"

void read_from_file(Particles& p, const std::string& filename){
    std::ifstream inp(filename.c_str(),std::ios::binary);
    assert(inp);
    const int size=N*sizeof(double);
    inp.read((char*)p.x.data(),size);
    inp.read((char*)p.y.data(),size);
    inp.read((char*)p.w.data(),size);
    inp.close();
}

#include <iostream>
using std::cout; using std::endl;
bool check_from_file(double res,int k, const std::string& filename){
    std::ifstream inp(filename.c_str());
    assert(inp);
    std::string line;
    for(int i=0;i<k;i++) std::getline(inp,line);
    double ref;
    inp>>ref>>ref;
    inp.close();
    bool check =  std::abs((ref-res)/ref) < 1e-6;
    if(not check)
    {
        cout.precision(15);
        cout<<"Computation error"<<endl;
        cout<<"Stream:\t"<<res<<endl;
        cout<<"File value:\t"<<ref<<endl;
    }
    return check;
}