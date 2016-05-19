#include <iostream>
#include <string>
#include "include/particles.h"
#include "morton_tree/build_tree.h"
using std::cout; using std::endl;
using std::vector;


void check_near(double a,double b);

int main() {
  constexpr int exp_order = 2;
  constexpr int Np=6;
  const int maxnodes= 2*Np;

  double x[Np]={0,0.25,0.75,0.25,0.75,1};
  double y[Np]={0,0.25,0.25,0.75,0.75,1};
  double w[Np]={0,1,1,1,1,0};
  Particles particles(x,y,w,Np);

  cout<<"N# of particles: "<<particles.N<<endl;

  //compute target locations with multipole expansion
  Tree tree(particles,maxnodes,exp_order);
  tree.computeMassAndExpansions<exp_order>();

  tree.PrintInfo(6);

  check_near(tree[0].r2,2*0.5*0.5);
  for(int i=1;i<5;i++) check_near(tree[i].r2,2*0.25*0.25);
}

void check_near(double a,double b){
  if(fabs(a-b)>1e-5) throw("check failed");
}


