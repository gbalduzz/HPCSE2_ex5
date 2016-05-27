#include <iostream>
#include <string>
#include "include/file_IO.h"
#include "include/particles.h"
#include "include/timing.h"
#include "morton_tree/build_tree.h"
#include "expansion/node_to_expansion.hpp"
using std::cout; using std::endl;
using std::vector;
void generateRandomData(Particles& ,int seed);
bool checkDifference(Particles&, Particles&);
bool checkResult(const Particles& targets,const int Np);

int main(int argc, char** argv) {
  constexpr int exp_order = 8;
  int Np=1e5;
  int Nt=1e5;
  if(argc==3){
    Np=atoi(argv[1]);
    Nt=atoi(argv[2]);
  }
  const int maxnodes= 4*Np/exp_order;
  Particles particles(Np),targets(Nt);
  generateRandomData(particles,0);
  generateRandomData(targets,42);
  assert(checkDifference(particles,targets));
  
  cout<<"N# of particles: "<<particles.N<<endl;
  cout<<"N# of targets: "<<targets.N<<endl;
  cout<<"parallelizing over "<<omp_get_max_threads()<<" threads"<<endl<<endl;
  
  reset_and_start_timer();
  //compute target locations with multipole expansion
  Tree<exp_order> tree(particles, maxnodes, 2*exp_order);
  ReorderIP(targets);
  potential<exp_order>(2., tree, targets);
  double dt = get_elapsed_mcycles();

  cout<<"Elapsed milolion cycles for Expansion: "<<dt<<endl;

  writeTime(omp_get_num_threads(),dt,"timing_e2p_e4_devel.out");
  Print(targets,5);

  /*
  //compute target locations with direct evaluations
  Profiler pr2("Direct evaluation");
  #pragma omp parallel for default(shared) schedule(static)
  for(int i=0;i<targets.N;i++) targets.w[i]=p2p(particles,targets.x[i],targets.y[i]);
  pr2.stop();
  writeToFile(targets,"direct.out");
  Print(targets,5);
  */ 
}

#include <random>
void generateRandomData(Particles& p,int seed){
  std::mt19937_64 mt(seed);
  std::uniform_real_distribution<double> ran(0,1);
  for(int i=0;i<p.N;i++){
    p.x[i]=ran(mt);
    p.y[i]=ran(mt);
    p.w[i]=2*ran(mt)-1;
  }
}

bool checkDifference(Particles& a, Particles& b){
  for(int i=0;i<a.N;i++)
    for(int j=0;j<b.N;j++)
      if((a.x[i]==b.x[i]) && (a.y[i]==b.y[i])) return false;
  return true;
}

bool checkResult(const Particles& targets,const int Np){
  if(Np != 1e5){cout<<"test not run \n"; return 0;}
  if(std::abs(targets.w[0]-(-7.92112))> 1e-4) throw("Test failed");
  if(std::abs(targets.w[1]-(-8.19657))> 1e-4) throw("Test failed");
  return 1;
}
