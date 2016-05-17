#pragma once
#include "p2p.h"
#include "morton_tree/node.h"
#include "e2p.h"

template<int k>
double evalPoint2Node(double x,double y,const double one_over_theta2,const Particles& prt, const Tree<k>& tree,const int id) {

  if (squareDistance(x, tree[id].xcom, y, tree[id].ycom) > one_over_theta2 * tree[id].r2)
    return e2p<k>(x - tree[id].xcom, y - tree[id].ycom, tree[id].re_expansion, tree[id].im_expansion);

  if (tree[id].child_id < 0) //is a leaf
    return p2p(prt.subEnsamble(tree[id].part_start,tree[id].occupancy()), x, y);

  double stream = 0;
  for (int i = 0; i < 4; i++) stream += evalPoint2Node(x, y, one_over_theta2 ,prt,tree, tree[id].child_id + i);
  return stream;
}

template<int exp_order>
void potential(const double theta,const Particles& particles, const Tree<exp_order>& tree, Particles& targets)
//compute the potential at targets and store the result in  targets.w array
{
  const double oot2=1/(theta*theta);
  for(int i=0;i<targets.N;i++)
  targets.w[i] = evalPoint2Node(targets.x[i],targets.y[i],oot2,particles,tree,0);
}
