#pragma once
#include <vector>
#include <assert.h>
#include <cmath>
#include "node.h"
#include "particles.h"
#include "morton_tree/tree_prepare.h"
#include "profiler.h"
#define LMAX 15
using std::vector;
using uint = unsigned int;
inline int decodeId(int x);

class Tree
{
  vector<Node> nodes;
  int   currnnodes=1, maxnodes;
  Particles p;
  vector<int> label;
  double ext, xmin, ymin;
  const int K; //max leaf capacity

public:

  Tree(Particles& unordered_particles,const int maxn, const int leaf_cap);
  void PrintInfo(int nprint);

private:
  void build();
  void build_leaf(const int nodeid, const int s, const int e, const double x0, const double y0, const double h);
  void build_tree(const int nodeid);
  void labelAndReorder(Particles&);
};

Tree::Tree(Particles& up,const int maxnd, const int leaf_cap):
maxnodes(maxnd),nodes(maxnodes),p(up.N),label(up.N),K(leaf_cap)
{ 
  labelAndReorder(up);
  nodes[0].setup(0, p.N, 0, 0);
#pragma omp parallel
#pragma omp single nowait
    build_tree(0);
}

void Tree::labelAndReorder(Particles &p_unord) {
  vector<int> keys(p.N);
  extent(p.N,p_unord.x,p_unord.y,xmin,ymin,ext);
  morton(p.N,p_unord.x,p_unord.y,xmin,ymin,ext,label.data());
  sort(p.N,label.data(),keys.data());
  reorder(p.N, keys.data(), p_unord.x, p_unord.y, p_unord.w, p.x, p.y, p.w);
}

void Tree::build_leaf(const int nodeid, const int s, const int e, const double x0, const double y0, const double h)
{
  Node * node = nodes.data() + nodeid;

  double w, wx, wy;
  leaf_setup(p.x + s, p.y + s, p.w + s, e - s, node->mass, w, wx, wy);

  node->xcom = w ? wx / w : (x0 + 0.5 * h);
  node->ycom = w ? wy / w : (y0 + 0.5 * h);
}

void Tree::build_tree(const int nodeid)
{
  Node * const node = nodes.data() + nodeid;

  const int s = node->part_start;
  const int e = node->part_end;
  const int l = node->level;
  const int mId = node->morton_id;

  const double h = ext / (1 << l);

  const double x0 = xmin + h * decodeId(mId);
  const double y0 = ymin + h * decodeId(mId >> 1);

  const bool leaf = e - s <= K || l + 1 > LMAX;

  if (leaf)
  {
    build_leaf(nodeid, s, e, x0, y0, h);
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

      const size_t indexmin = c == 0 ? s : lower_bound_vec(s, e, key1, label.data());
      const size_t indexsup = c == 3 ? e : upper_bound_vec(s, e, key2, label.data());

      const int chId = childbase + c;
      nodes[chId].setup(indexmin, indexsup, l + 1, key1);

#pragma omp task firstprivate(chId) if (indexsup - indexmin > 5e3 && c < 3)
      {
        build_tree(chId);
      }
    }
  }
}

inline int decodeId(int x)
{
  x &= 0x55555555;
  x = (x ^ (x >>  1)) & 0x33333333;
  x = (x ^ (x >>  2)) & 0x0f0f0f0f;
  x = (x ^ (x >>  4)) & 0x00ff00ff;
  x = (x ^ (x >>  8)) & 0x0000ffff;
  return x;
}

#include <iostream>
using std::cout; using std::endl;
void Tree::PrintInfo(int nprint){
  nprint = std::min(nprint,currnnodes);
  cout<<"Number of nodes: "<<currnnodes<<endl;
  for(int i=currnnodes-nprint; i<currnnodes;i++){
    const int s =nodes[i].part_start;
    const int e =nodes[i].part_end;
    cout<<"i="<<i<<" level: "<<nodes[i].level<<"\tfirst child: "<<nodes[i].child_id<<
    //"\t N_points: "<<e-s<<
    "\t mass: "<<nodes[i].mass<<
    "\t start,end "<<s<<" , "<<e<<
    "\tcom: "<<nodes[i].xcom<<" , "<<nodes[i].ycom<<endl;
  }
}
