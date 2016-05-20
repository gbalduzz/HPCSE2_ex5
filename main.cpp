#define DETAILED_PROFILING
#include <iostream>
#include <string>
#include "include/Performance_measurer.h"
#include "include/file_IO.h"
#include "include/particles.h"
#include "include/profiler.h"
#include "morton_tree/build_tree.h"
#include "expansion/node_to_expansion.hpp"
using std::cout; using std::endl;
using std::vector;
void generateRandomData(Particles& ,int seed);
bool checkDifference(Particles&, Particles&);

int main(int argc, char** argv) {
  constexpr int exp_order = 8;
  int Np=1e4;
  int Nt=1e4;
  if(argc==3){
    Np=atoi(argv[1]);
    Nt=atoi(argv[2]);
  }
  const int maxnodes= 4*Np/exp_order;
  Particles particles(Np),targets(Nt);
  generateRandomData(particles,0);
  generateRandomData(targets,42);
  assert(checkDifference(particles,targets));
  //targets.x[0]=1;targets.y[0]=1;

  cout<<"N# of particles: "<<particles.N<<endl;
  cout<<"N# of targets: "<<targets.N<<endl;
  cout<<"parallelizing over "<<omp_get_max_threads()<<" threads"<<endl<<endl;

  //compute target locations with multipole expansion
  {
    Profiler pr("Expansion");
#ifdef DETAILED_PROFILING
    Profiler p2("Tree creation");
#endif
    Tree tree(particles, maxnodes, exp_order);
    tree.computeMassAndExpansions<exp_order>();
#ifdef DETAILED_PROFILING
    p2.stop();
#endif
    ReorderIP(targets);
#ifdef DETAILED_PROFILING
    Profiler p3("potentail eval");
#endif
    potential<exp_order>(2., tree, targets);
  }
  writeToFile(targets,"expansion.out");
  Print(targets,5);

  //compute target locations with direct evaluations
  Profiler pr2("Direct evaluation");
#pragma omp parallel for default(shared) schedule(static)
  for(int i=0;i<targets.N;i++) targets.w[i]=p2p(particles,targets.x[i],targets.y[i]);
  pr2.stop();
  writeToFile(targets,"direct.out");
  Print(targets,5);

  //tree.PrintInfo(10);
}


#include <random>
void generateRandomData(Particles& p,int seed){
  std::mt19937_64 mt(seed);
  std::uniform_real_distribution<double> ran(0,1);
  for(int i=0;i<p.N;i++){
    p.x[i]=ran(mt);
    p.y[i]=ran(mt);
    p.w[i]=ran(mt);//2*ran(mt)-1;
  }
}

bool checkDifference(Particles& a, Particles& b){
  for(int i=0;i<a.N;i++)
    for(int j=0;j<b.N;j++)
      if((a.x[i]==b.x[i]) && (a.y[i]==b.y[i])) return false;
  return true;
}
