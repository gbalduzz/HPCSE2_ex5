#pragma once
#include "p2p.h"
#include "morton_tree/node.h"
#include "e2p.h"

//TODO MAYBE: spawn children only if needed
template<int k>
double evalPoint2Node(double x,double y,const double one_over_theta2,const Particles& prt, const Tree& tree,const int id) {

  if (squareDistance(x, tree[id].xcom, y, tree[id].ycom) > one_over_theta2 * tree[id].r2)
    return e2p<k>(x - tree[id].xcom, y - tree[id].ycom, tree.getReExpansion(id), tree.getImExpansion(id));

  if (tree[id].child_id == 0) { //is a leaf
    const int s = tree[id].part_start;
    const int e = tree[id].part_end;
    return p2p(prt.subEnsamble(s,e-s), x, y);
  }

  double stream = 0;
  for (int i = 0; i < 4; i++) stream += evalPoint2Node<k>(x, y, one_over_theta2 ,prt,tree, tree[id].child_id + i);
  return stream;
}

template<int exp_order>
void potential(const double theta,const Tree& tree, Particles& targets)
//compute the potential at targets and store the result in  targets.w array
{
  assert(tree.exp_order == exp_order);
  const double oot2=1/(theta*theta);
  for(int i=0;i<targets.N;i++)
  targets.w[i] = evalPoint2Node<exp_order>(targets.x[i],targets.y[i],oot2,tree.getParticles(),tree,0);
}
