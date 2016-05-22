#include <iostream>
#include <morton_tree/build_tree.h>
using std::cout; using std::endl;
void PrintExp(const double* re,const double* im,const int N);
void generateRandomData(Particles& ,int seed);

int main(){
  constexpr int exp_order = 4;
  int Np=10;
  int Nt=1;
  const int maxnodes= 4*Np/exp_order;
  Particles particles(Np),targets(Nt);
  generateRandomData(particles,0);
  generateRandomData(targets,42);
  //targets.x[0]=1;targets.y[0]=1;

  cout<<"N# of particles: "<<particles.N<<endl;
  cout<<"N# of targets: "<<targets.N<<endl;
  cout<<"parallelizing over "<<omp_get_max_threads()<<" threads"<<endl<<endl;

  Tree<exp_order> tree(particles, maxnodes, 2*exp_order);

  //  ReorderIP(targets);
   // Profiler p3("potentail eval");
   // potential<exp_order>(2., tree, targets);
  tree.PrintInfo(tree.size());
  cout<<"Node 0:\n";
  PrintExp(tree.getReExpansion(0),tree.getImExpansion(0),exp_order);
  for(int i= 0;i<4;i++){
    const int child_id = tree[0].child_id + i;
    cout<<"child "<<i+1<<endl;
     PrintExp(tree.getReExpansion(child_id),tree.getImExpansion(child_id),exp_order);
  }

  constexpr int k=2;
  double c_re[k+1]={1,0.5,0.25};
  double c_im[k+1]={0,0.5,0.25};
  double new_re[k+1]={0,0,0};
  double new_im[k+1]={0,0,0};
  e2e<k>(c_re,c_im,0,0,new_re,new_im);
  cout<<"Recomputed expansion\n";
  PrintExp(new_re,new_im,k);

  return 0;
}

void PrintExp(const double* re,const double* im,const int N){
  cout<<"Re: ";
  for(int i=0;i<N+1;i++) cout<<re[i]<<"\t";
  cout<<endl<<"Im: ";
  for(int i=0;i<N+1;i++) cout<<im[i]<<"\t";
  cout<<endl;
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
