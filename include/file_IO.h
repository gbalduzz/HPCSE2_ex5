#pragma once
#include <string>
#include <cmath>
#include <fstream>
#include <assert.h>
#include "particles.h"
#include <iostream>
using std::cout; using std::endl;

void read_from_file(const std::string& filename, Particles& particles, Particles& targets){
  std::ifstream inp(filename.c_str(),std::ios::binary);
  assert(inp);
  int N;
  inp.read((char *) &N,sizeof(int));
  particles.resize(N);
  int size=N*sizeof(double);
  inp.read((char*)particles.x,size);
  inp.read((char*)particles.y,size);
  inp.read((char*)particles.w,size);

  inp.read((char *) &N,sizeof(int));
  targets.resize(N);
  size = N*sizeof(double);
  inp.read((char*)targets.x,size);
  inp.read((char*)targets.y,size);
  inp.read((char*)targets.w,size);

  inp.close();
  }

void writeToFile(const Particles& p, const std::string& filename){
  std::ofstream out(filename.c_str());
  if(not out)throw("output file not available");
  for(int i=0;i<p.N;i++) out<<p.x[i]<<"\t"<<p.y[i]<<"\t"<<p.w[i]<<endl;
}

/*using std::cout; using std::endl;
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
  *  return check;
}*/


