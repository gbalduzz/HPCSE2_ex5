#pragma once
#include <vector>
#include <assert.h>
#include <cmath>
#include "morton_tree/node.h"
#include "particles.h"
#include "profiler.h"
#include "helper_methods.h"
#include "expansion/P2E.h"
#include "expansion/e2e.h"
#include "profiler.h"
using std::vector;

template<int order>
class Tree
{
  double ext, xmin, ymin;
  int maxnodes;
  int   currnnodes=1;
  const int K; //max leaf capacity
  vector<Node> nodes;
  Particles p;
  vector<uint> label;
  vector<double> re_expansions;
  vector<double> im_expansions;

public:

  Tree(Particles& unordered_particles,const int maxn, const int leaf_cap);
  inline const Node& operator[](const int i)const{return nodes[i];}
  inline const Particles& getParticles()const{return p;}
  inline const double* getReExpansion(const int i)const{assert(i<currnnodes); return &re_expansions[i*(order+1)];}
  inline const double* getImExpansion(const int i)const{assert(i<currnnodes); return &im_expansions[i*(order+1)];}
  inline  double* getReExpansion(const int i){assert(i<currnnodes); return &re_expansions[i*(order+1)];}
  inline  double* getImExpansion(const int i){assert(i<currnnodes); return &im_expansions[i*(order+1)];}
  void PrintInfo(int nprint);
  int size()const{return currnnodes;}

private:
  void build();
  void build_leaf(const int nodeid, const int s, const int e);
  void build_tree(const int nodeid);
  void labelAndReorder(Particles&);
  void computeMassAndExpansion();
};

template<int order>
Tree<order>::Tree(Particles& up,const int maxnd, const int leaf_cap):
maxnodes(maxnd),K(leaf_cap),nodes(maxnodes),p(up.N),label(up.N),
re_expansions(maxnodes*(order+1)),im_expansions(maxnodes*(order+1))
{ 
  labelAndReorder(up);
  nodes[0].setup(0, p.N, 0, 0);
#pragma omp parallel
#pragma omp single nowait
    build_tree(0);
//compute parent nodes properties
  //TODO parallelize over threads
  computeMassAndExpansion();
}

template<int order>
void Tree<order>::labelAndReorder(Particles &p_unord) {
  vector<int> keys(p.N);
  extent(p.N,p_unord.x,p_unord.y,xmin,ymin,ext);
  morton(p.N,p_unord.x,p_unord.y,xmin,ymin,ext,label.data());
  sort(p.N,label.data(),keys.data());
  reorder(p.N, keys.data(), p_unord.x, p_unord.y, p_unord.w, p.x, p.y, p.w);
}

template<int order>
void Tree<order>::build_leaf(const int nodeid, const int s, const int e)
{
  Node * node = nodes.data() + nodeid;

  double wx, wy;
  leaf_setup(p.x + s, p.y + s, p.w + s, e - s, node->mass, wx, wy);
  const double w = node->mass;
  if(w) {
    node->xcom = wx / w;
    node->ycom = wy / w;
  }
}

template<int order>
void Tree<order>::build_tree(const int nodeid)
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
    P2E<order>(p.subEnsamble(s, e-s),
               node->xcom, node->ycom, getReExpansion(nodeid), getImExpansion(nodeid));
  }
  else
  {
    int childbase;
#pragma omp atomic capture
    {
      childbase = currnnodes; currnnodes += 4;
    }
    assert(nodeid < childbase);
    if(childbase + 4 > maxnodes) throw(std::logic_error("Nodes exceed maximum"));

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

template<int order>
void Tree<order>::computeMassAndExpansion()
{
    for (int i = currnnodes - 1; i > -1; i--) {
      Node *const node = &nodes[i];
      const int first_child = node->child_id;
      if (first_child) {
        assert(node->mass == 0);
        //combine children coms
        double mass(0), xcom(0), ycom(0);
        for (int j = 0; j < 4; j++) {
          const int child_id = first_child+j;
          const double mass_term = nodes[child_id].mass;
          mass += mass_term;
          xcom += mass_term * nodes[child_id].xcom;
          ycom += mass_term * nodes[child_id].ycom;
        }
        node->setCom(mass, xcom / mass, ycom / mass);
        //combine children expansions
        for(int j=0;j<4;j++) {
          const int child_id = first_child+j;
          const double z0_re = nodes[child_id].xcom - node->xcom;
          const double z0_im = nodes[child_id].ycom - node->ycom;
          e2e<order>(getReExpansion(child_id), getImExpansion(child_id), z0_re, z0_im,
                     getReExpansion(i), getImExpansion(i));
        }
      }
      else if (node->mass == 0) continue;
      const double h = ext / (1 << node->level);
      const uint mId = node->morton_id;
      const double step = ext / (1 << LMAX);
      const double x0 = xmin + step * decodeId(mId);
      const double y0 = ymin + step * decodeId(mId >> 1);
      //compute radius
      node->r2 = computeRadius(x0, y0, h, node->xcom, node->ycom);
    }
}


#include <iostream>
using std::cout; using std::endl;
template<int order>
void Tree<order>::PrintInfo(int nprint){
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
