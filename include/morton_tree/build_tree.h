#pragma once
#include <vector>
#include <assert.h>
#include <cmath>
#include "morton_tree/node.h"
#include "particles.h"
#include "profiler.h"
#include "helper_methods.h"
#include "expansion/P2E.h"
#include "profiler.h"
using std::vector;

class Tree
{
  double ext, xmin, ymin;
  int maxnodes;
  int   currnnodes=1;
  const int K; //max leaf capacity
  vector<Node> nodes;
  vector<double> re_expansions;
  vector<double> im_expansions;
  Particles p;
  vector<uint> label;

public:

  Tree(Particles& unordered_particles,const int maxn, const int leaf_cap);
  template<int order>
  void computeMassAndExpansions();
  inline const Node& operator[](const int i)const{return nodes[i];}
  inline const Particles& getParticles()const{return p;}
  inline const double* getReExpansion(const int i)const{assert(i<currnnodes); return &re_expansions[i*(exp_order+1)];}
  inline const double* getImExpansion(const int i)const{assert(i<currnnodes); return &im_expansions[i*(exp_order+1)];}
  void PrintInfo(int nprint);
  int size()const{return currnnodes;}

  int exp_order=0;
private:
  void build();
  void build_leaf(const int nodeid, const int s, const int e);
  void build_tree(const int nodeid);
  void labelAndReorder(Particles&);
};


Tree::Tree(Particles& up,const int maxnd, const int leaf_cap):
maxnodes(maxnd),K(leaf_cap),nodes(maxnodes),p(up.N),label(up.N)
{ 
  labelAndReorder(up);
  nodes[0].setup(0, p.N, 0, 0);
#pragma omp parallel
#pragma omp single nowait
    build_tree(0);
}


void Tree::labelAndReorder(Particles &p_unord) {
#ifdef DETAILED_PROFILING
  Profiler pr("label and reorder");
#endif
  vector<int> keys(p.N);
  extent(p.N,p_unord.x,p_unord.y,xmin,ymin,ext);
  morton(p.N,p_unord.x,p_unord.y,xmin,ymin,ext,label.data());
  sort(p.N,label.data(),keys.data());
  reorder(p.N, keys.data(), p_unord.x, p_unord.y, p_unord.w, p.x, p.y, p.w);
}


void Tree::build_leaf(const int nodeid, const int s, const int e)
{
  Node * node = nodes.data() + nodeid;

  double wx, wy;
  leaf_setup(p.x + s, p.y + s, p.w + s, e - s, node->mass, wx, wy);
  const double w = node->mass;
  node->xcom = w ? wx / w : 0;
  node->ycom = w ? wy / w : 0;
}


void Tree::build_tree(const int nodeid)
{
  Node * const node = nodes.data() + nodeid;

  const int s = node->part_start;
  const int e = node->part_end;
  const int l = node->level;
  const int mId = node->morton_id;

  const bool leaf = e - s <= K || l + 1 > LMAX;

  if (leaf)
  {
    build_leaf(nodeid, s, e);
  }
  else
  {
    int childbase;
#pragma omp atomic capture
    {
      childbase = currnnodes; currnnodes += 4;
    }
    assert(nodeid < childbase);
    assert(childbase + 4 <= maxnodes);

    node->child_id = childbase;

    for(int c = 0; c < 4; ++c)
    {
      const int shift = 2 * (LMAX - l - 1);

      const int key1 = mId | (c << shift);
      const int key2 = key1 + (1 << shift) - 1;

      const size_t indexmin = c == 0 ? s : lower_bound_vec(s, e, key1, label);
      const size_t indexsup = c == 3 ? e : upper_bound_vec(s, e, key2, label);

      const int chId = childbase + c;
      nodes[chId].setup(indexmin, indexsup, l + 1, key1);

#pragma omp task firstprivate(chId) if (indexsup - indexmin > 1e3 && c < 3)
      {
        build_tree(chId);
      }
    }
  }
}

template<int exp_order>
void Tree::computeMassAndExpansions() {
  this->exp_order = exp_order;
  re_expansions.resize(currnnodes*(exp_order+1));
  im_expansions.resize(currnnodes*(exp_order+1));
  for(int i=currnnodes-1;i>-1;i--){
    Node* const node = &nodes[i];
    const int child_id =node->child_id;
    if(child_id){ //combine childrens coms
      assert(node->mass==0);
      double mass(0),xcom(0),ycom(0);
      for(int j=0;j<4;j++) {
        const double mass_term = nodes[child_id+j].mass;
        mass += mass_term;
        xcom += mass_term * nodes[child_id+j].xcom;
        ycom += mass_term * nodes[child_id+j].ycom;
      }
      node->setCom(mass,xcom/mass,ycom/mass);
    }
    else if(node->mass == 0) continue;
    const double h = ext / (1 << node->level);
    const uint mId = node->morton_id;
    const double step = ext/(1 << LMAX);
    const double x0 = xmin + step * decodeId(mId);
    const double y0 = ymin + step * decodeId(mId >> 1);
    //compute radius
    node->r2 = computeRadius(x0,y0,h,node->xcom,node->ycom);
    //compute expansion
    const int offset = (exp_order+1)*i;
    const int s = node->part_start;
    const int e = node->part_end;
    P2E<exp_order>(p.subEnsamble(s,e-s),
                   node->xcom,node->ycom,&re_expansions[offset],&im_expansions[offset]);
  }
}


#include <iostream>
using std::cout; using std::endl;

void Tree::PrintInfo(int nprint){
  nprint = std::min(nprint,currnnodes);
  cout<<"Number of nodes: "<<currnnodes<<endl;
  for(int i=0; i<nprint;i++){
    const int s =nodes[i].part_start;
    const int e =nodes[i].part_end;
    cout<<"i="<<i<<" level: "<<nodes[i].level<<"\tfirst child: "<<nodes[i].child_id<<
    "\t r2: "<<nodes[i].r2<<
    "\t mass: "<<nodes[i].mass<<
    "\t start,end "<<s<<" , "<<e<<
    "\tcom: "<<nodes[i].xcom<<" , "<<nodes[i].ycom<<endl;
  }
}

void ReorderIP(Particles &p_unord) {
  const int N=p_unord.N;
  Particles p(N);
  vector<int> keys(N);
  vector<uint> label(N);
  double ext,xmin,ymin;
  extent(N,p_unord.x,p_unord.y,xmin,ymin,ext);
  morton(N,p_unord.x,p_unord.y,xmin,ymin,ext,label.data());
  sort(N,label.data(),keys.data());
  reorder(N, keys.data(), p_unord.x, p_unord.y, p_unord.w, p.x, p.y, p.w);
  swap(p_unord,p);
}
