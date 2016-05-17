//
// Created by giovanni on 23.03.16.
//
#pragma once
#include <vector>
#include <assert.h>
#include <cmath>
#include "expansion/P2E.h"
#include "node.h"
#include "extent.h"
#include "sort.h"
#include "morton.h"
#include "reorder.h"
#include "profiler.h"
#include "p2p.h"
using std::vector;
template<int order>
void create_children(const int parent_id, Tree<order>& tree,
                     const Particles& p,const uint* label,const int N);
template<int k>
void computeRootComExpansion(Tree<k>& tree, const Particles& p);

template<int exp_order>
void buildTree(const Particles& p, const int k, Particles& p_ordered,Tree<exp_order>& tree)
{
  //perform the exrcise 2 and find the morton labels of the points
  const int N=p.N;
  uint *keys = new uint[N];
  uint *label = new uint[N];
  double xmin, ymin, ext;
  {
#ifdef DEBUG
    Profiler p("old morton labelling");
#endif
    extent(N, p.x, p.y, xmin, ymin, ext);
    morton(N, p.x, p.y, xmin, ymin, ext, label);
    sort(N, label, keys);
    reorder(N, keys, p.x, p.y, p.w, label, p_ordered.x, p_ordered.y, p_ordered.w);
  }
  delete[] keys;

  tree.resize(1);
  tree.reserve(N*2);
  //create root node
  tree[0].level = 0;
  tree[0].morton_id = 0;
  tree[0].part_start = 0;
  tree[0].part_end = N - 1;

  const int max_level = sizeof(int) * 4; //number of bits over 2
  //if number of particles > k  : split
  for (int i = 0; i < tree.size(); i++) {
    if (tree[i].occupancy() > k) {
      if (tree[i].level == max_level) { //no more space for branching
        std::cout << "Warning: No more space for branching" << std::endl;
        tree[i].child_id = -5;
        continue;
      }
      tree[i].child_id = tree.size();
      create_children(i, tree, p_ordered, label, N);
    }
  }
  computeRootComExpansion(tree,p_ordered);

//compute radius;
  tree[0].r2=0;
  const double x_corner[]={xmin,xmin+ext,xmin,xmin+ext};
  const double y_corner[]={ymin,ymin,ymin+ext,ymin+ext};
  for(int i=0;i<4;i++){
    const double dist=squareDistance(tree[0].xcom,x_corner[i],tree[0].ycom,y_corner[i]);
    if (tree[0].r2<dist) tree[0].r2=dist;
  }

  delete[] label;
}

inline int get_new_id(uint parent_id,int level,int i)
{
  static const uint n_bits=sizeof(int)*8;
  assert(2*level <= n_bits);
  return parent_id  | i<< (n_bits-2*level);
}

uint create_mask(int level)
//create a mask that ignores the bits in the morton index corresponding to a higher level
// the mask is 2*level 1s followed by 0s
{
  static const uint n_bits=sizeof(int)*8;
  assert(2*level <= n_bits);
  return 2*level==n_bits ? -1 : //otherwise returns all 0 instead of all 1
         ((1 << 2*level)-1) << (n_bits-2*level);
}

template<int order>
void create_children(const int parent_id, Tree<order>& tree,
                     const Particles& p,const uint* label,const int N)
{
  int current_idx=tree[parent_id].part_start;
  uint mask=create_mask(tree[parent_id].level+1);
  for(int i=0;i<4;i++) {
    tree.push_back(Node<order>());
    int child_id=tree.size()-1;
    tree[child_id].level=tree[parent_id].level+1;
    tree[child_id].morton_id=get_new_id(tree[parent_id].morton_id,tree[child_id].level,i);
    //find points inside node
    if((label[current_idx] & mask) == tree[child_id].morton_id){//branch is not empty
      tree[child_id].part_start=current_idx;
      while((label[current_idx] & mask) == tree[child_id].morton_id
            &&  current_idx<N) {
        tree[child_id].xcom+=p.w[current_idx]*p.x[current_idx];
        tree[child_id].ycom+=p.w[current_idx]*p.y[current_idx];
        tree[child_id].mass+=p.w[current_idx];
        current_idx++;
      }
      tree[child_id].part_end=current_idx-1;
      tree[child_id].xcom/=tree[child_id].mass;
      tree[child_id].ycom/=tree[child_id].mass;
      //radius. Naive implementation
      tree[child_id].r2=0;
      for(int i=tree[child_id].part_start;i<=tree[child_id].part_end;i++){
        const double dist=squareDistance(p.x[i],tree[child_id].xcom,p.y[i],tree[child_id].ycom);
        if(tree[child_id].r2<dist) tree[child_id].r2=dist;
      }

      //compute expansion
      const int offset = tree[child_id].part_start;
      const int length = tree[child_id].part_end-offset;
      P2E<order>(p.subEnsamble(offset,length),
                 tree[child_id].xcom, tree[child_id].ycom,
                 tree[child_id].re_expansion,tree[child_id].im_expansion);
    }
    else{//is empty
      tree[child_id].part_start=-1;
      tree[child_id].part_end=-2;
    }
  }
}

template<int k>
void computeRootComExpansion(Tree<k>& tree, const Particles& p){
  for(int i=1;i<5;i++){
    tree[0].mass+=tree[i].mass;
    tree[0].xcom+=tree[i].mass*tree[i].xcom;
    tree[0].ycom+=tree[i].mass*tree[i].ycom;
  }
  tree[0].xcom/=tree[0].mass;
  tree[0].ycom/=tree[0].mass;
  //compute expansion
  P2E<k>(p,
             tree[0].xcom, tree[0].ycom,
             tree[0].re_expansion,tree[0].im_expansion);
}
