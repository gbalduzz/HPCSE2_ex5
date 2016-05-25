#pragma once
#include "p2p.h"
#include "morton_tree/node.h"
#include "e2p.ispc.h"
constexpr int  BUFSIZE=16;

template<int ord>
void evaluate(double& result,const double xt,const double yt,
              const double theta, const Tree<ord>& tree) {
  const double one_over_theta2=1./(theta*theta);
  int stack[LMAX * 3];
  int bufcount = 0;
  result = 0;
  double rzs[BUFSIZE], izs[BUFSIZE]/*,ws[BUFSIZE]*/;
  const double *rxps[BUFSIZE], *ixps[BUFSIZE];

  int stackentry = 0, maxentry = 0;
  stack[0] = 0;
  result = 0;
  while (stackentry > -1) {
    const int nodeid = stack[stackentry--];
    const Node *const node = &tree[nodeid];
    //double tmp[2];
    const double r2 = squareDistance(xt, node->xcom, yt, node->ycom);
    if (node->r2 < one_over_theta2 * r2) {
      rzs[bufcount] = xt - node->xcom;
      izs[bufcount] = yt - node->ycom;
      rxps[bufcount] = tree.getReExpansion(nodeid);
      ixps[bufcount] = tree.getImExpansion(nodeid);
      bufcount++;
      if (bufcount == BUFSIZE) {
        bufcount = 0;
        result += ispc::e2p(rzs, izs, rxps, ixps, BUFSIZE,ord);
      }
    }
    else {
      if (not node->child_id){
        const int s= node->part_start;
        const int e= node->part_end;
        result+=p2p(tree.getParticles().subEnsamble(s,e-s),xt,yt);
      }
      else{
        for(int i=0;i<4;i++) stack[++stackentry] = node->child_id+i;
        maxentry= std::max(maxentry,stackentry);
      }
    }
  }
  if(bufcount) result+= ispc::e2p(rzs, izs, rxps, ixps, bufcount, ord);
}

template<int exp_order>
void potential(const double theta,const Tree<exp_order>& tree, Particles& targets)
//compute the potential at targets and store the result in  targets.w array
{
#pragma omp parallel for schedule(static,1)
  for(int i=0;i<targets.N;i++)
    evaluate<exp_order>(targets.w[i],targets.x[i],targets.y[i],theta,tree);
}
