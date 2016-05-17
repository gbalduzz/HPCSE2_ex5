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


int main(int argc, char** argv) {
  constexpr int exp_order = 8;

  Particles particles,targets;
  read_from_file( "/home/giovanni/uni/HPCSE/ex5/diegoBinaryN400",//PROJECT_SOURCE_DIR
                  particles,targets);

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

