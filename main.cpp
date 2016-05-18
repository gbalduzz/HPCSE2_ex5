const static int k=20;
#include <iostream>
#include <string>
#include "include/Performance_measurer.h"
#include "include/file_IO.h"
#include "include/particles.h"
#include "include/profiler.h"
#include "morton_tree/buildTree.h"
#include "expansion/node_to_expansion.hpp"
using std::cout; using std::endl;
using std::vector;
void generateRandomData(Particles& ,int seed);
bool checkDifference(Particles&, Particles&);

int main(int argc, char** argv) {
  constexpr int exp_order = 8;
  int Np=10000;
  int Nt=1000;
  if(argc==3){
    Np=atoi(argv[1]);
    Nt=atoi(argv[2]);
  }
  Particles particles(Np),targets(Nt);
  generateRandomData(particles,0);
  generateRandomData(targets,42);
  assert(checkDifference(particles,targets));

  cout<<"N# of particles: "<<particles.N<<endl;
  cout<<"N# of targets: "<<targets.N<<endl;

  //compute target locations with multipole expansion
  Profiler pr("Expansion");
  Particles p_ordered(particles.N);
  Tree<exp_order> tree;
  buildTree<exp_order>(particles,exp_order,p_ordered,tree);
  potential(2.,p_ordered,tree,targets);
  pr.stop();
  writeToFile(targets,"expansion.out");

  //compute target locations with direct evaluations
  Profiler pr2("Direct evaluation");
#pragma omp parallel for default(shared) schedule(static)
  for(int i=0;i<targets.N;i++) targets.w[i]=p2p(p_ordered,targets.x[i],targets.y[i]);
  pr2.stop();
  writeToFile(targets,"direct.out");

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
