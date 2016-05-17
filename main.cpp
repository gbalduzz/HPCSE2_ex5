const static int k=20;
#include <iostream>
#include <string>
#include "include/Performance_measurer.h"
#include "include/file_IO.h"
#include "include/particles.h"
#include "include/profiler.h"
#include "include/p2p.h"
using std::cout; using std::endl;
using std::vector;


int main(int argc, char** argv) {
  Particles particles,targets;
  read_from_file( "/home/giovanni/uni/HPCSE/ex5/diegoBinaryN400",//PROJECT_SOURCE_DIR
                  particles,targets);

  cout<<"N# of particles: "<<particles.N<<endl;
  cout<<"N# of targets: "<<targets.N<<endl;

  cout<<"\nParticles:\n";
  for(int i=0;i<10;i++){
    cout<<particles.x[i]<<"\t"<<particles.y[i]<<"\t"<<particles.w[i]<<endl;
  }
  cout<<"\nTargets:\n";
  for(int i=0;i<10;i++){
    cout<<targets.x[i]<<"\t"<<targets.y[i]<<"\t"<<targets.w[i]<<endl;
  }
  cout<<endl;


  //compute some expansion directly with p2p
  targets.N=10;
  for(int i=0;i<targets.N;i++){
    cout<<"*"; cout.flush();
    targets.w[i]=p2p(particles,targets.x[i],targets.y[i]);
  }
  cout<<endl;
  writeToFile(targets,"p2p.dat");
}

