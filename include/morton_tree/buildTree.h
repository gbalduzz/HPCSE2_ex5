//
// Created by giovanni on 23.03.16.
//

#ifndef EX3_BUILD_H
#define EX3_BUILD_H
#include <vector>
#include <assert.h>
#include <cmath>
#include "node.h"
#include "extent.h"
#include "sort.h"
#include "morton.h"
#include "reorder.h"
#include "find_last.h"

//max depth at which tasks are spawned
#define PARALLEL_MAX_DEPTH 2

using std::vector;
using Tree = vector<Node>;
void create_children_recursively(const int parent_id,vector<Node>&tree,const double* x,const double* y,const double* m,const uint* label,const int N,const int k);
void compute_com(vector<Node>& tree,double *x,double* y,double* m);
inline int check_size(const vector<Node>& tree);

void buildTree(const Particles& p, const int k, Particles p_ordered,Tree& tree)
{
  //perform the exrcise 2 and find the morton labels of the points
  const int N=p.N;
  uint *keys = new uint[N];
  uint *label = new uint[N];
  double xmin, ymin, ext;
  {
    Profiler p("old morton labelling");
    extent(N, p.x, p.y, xmin, ymin, ext);
    morton(N, p.x, p.y, xmin, ymin, ext, label);
    sort(N, label, keys);
    reorder(N, keys, p.x, p.y, p.w, label, p_ordered.x, p_ordered.y, p_ordered.w);
  }
  delete[] keys;

  tree.clear();
  tree.reserve(N*2);
  {
    Profiler p("Tree creation");
    //create root tree[id]
    tree[0].level = 0;
    tree[0].morton_id = 0;
    tree[0].part_start = 0;
    tree[0].part_end = N - 1;

#pragma omp parallel
    {
#pragma omp single
      {
        Profiler p2("children recursion");
        {create_children_recursively(0, tree, p_ordered.x, p_ordered.y, p_ordered.w, label, N, k);}
#pragma omp taskwait
        p2.stop();
        Profiler p3("mass recursion");
        {compute_com(tree, xsorted, ysorted, mass_sorted);}
#pragma omp taskwait
        p3.stop();
      }
    }//end omp parallel

  }//end Profiler
  delete[] label;
  return tree;
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
  static const uint n_bits=sizeof(uint)*8;
  assert(2*level <= n_bits);
  return 2*level==n_bits ? -1 : //otherwise returns all 0 instead of all 1
         ((1 << 2*level)-1) << (n_bits-2*level);
}

void create_children_recursively(const int parent_id,vector<Node>&tree,const double* x,const double* y,const double* m,const uint* label,const int N,const int k) {
  static const int max_level = sizeof(int) * 4; //number of bits over 2
  static const int max_thread_depth=std::ceil(std::log2(omp_get_max_threads())/2.);
  if (tree[parent_id].occupancy() <= k) return;
  if (tree[parent_id].level == max_level) { //no more space for branching
    std::cout << "Warning: No more space for branching" << std::endl;
    tree[parent_id].child_id = -5;
    return;
  }
  uint mask = create_mask(tree[parent_id].level + 1);
  int first_child_id,start,end;

  //create enough space in the tree in a thread safe manner
#pragma omp critical
  {
    first_child_id = tree.size();
    for (int i = 0; i < 4; i++) tree.push_back(Node());
  }

  tree[parent_id].child_id = first_child_id;
  bool done=false;//flag to exit when a child occupies the whole parent
  const int index_order[] = {0, 3, 1, 2};//traverse children in this order
  for (int i : index_order) {
    int child_id = first_child_id + i;
    tree[child_id].level = tree[parent_id].level + 1;
    tree[child_id].morton_id = get_new_id(tree[parent_id].morton_id, tree[child_id].level, i);
    //find points inside tree[child_id]
    if(not done) switch (i) {

        case 0:
          tree[child_id].part_start = tree[parent_id].part_start;
          tree[child_id].part_end = find_last(tree[parent_id].part_start, tree[parent_id].part_end, label, mask,
                                              tree[child_id].morton_id);
          if(tree[child_id].part_end == tree[parent_id].part_end) done=true;
          break;

        case 3:  //could use result of 0, but is just an additional iteration
          tree[child_id].part_start = find_first(tree[parent_id].part_start, tree[parent_id].part_end, label,
                                                 mask, tree[child_id].morton_id);
          tree[child_id].part_end = tree[parent_id].part_end;
          if(tree[child_id].part_start == tree[parent_id].part_start) done=true;
          break;
        case 1:
          start= tree[first_child_id].part_end==-1 ?  tree[parent_id].part_start : tree[first_child_id].part_end+1;
          tree[child_id].part_start = start;
          tree[child_id].part_end = find_last(start, tree[parent_id].part_end, label, mask,
                                              tree[child_id].morton_id);
          break;
        default: //last one is determined by the others
          end=  tree[first_child_id+3].part_start==-1 ? tree[parent_id].part_end : tree[first_child_id+3].part_start-1;
          start=tree[first_child_id+1].part_end==-1   ?  tree[first_child_id].part_end==-1 ? tree[parent_id].part_start
                                                                                           : tree[first_child_id].part_end+1
                                                      : tree[first_child_id+1].part_end+1;
          if(end>=start) {//not empty
            tree[child_id].part_start = start;
            tree[child_id].part_end = end;
          }
          break;
      }
  }


  for (int i = 0; i < 4; i++)//iterate on children
#pragma omp task shared(tree) if(tree[parent_id].level<PARALLEL_MAX_DEPTH)
  { create_children_recursively(first_child_id + i, tree, x, y, m, label, N, k); }
}

void solve_dependencies(vector<Node>& tree,int id,double*x,double*y,double*m);

void compute_com(vector<Node>& tree,double *x,double* y,double* m){
  assert(tree.size()>=5); //at least one branching
  for(int i=1;i<5;i++)
#pragma omp task default(shared)
  {solve_dependencies(tree,i,x,y,m);}
#pragma omp taskwait

  for(int i=1;i<5;i++){
    tree[0].mass+=tree[i].mass;
    tree[0].xcom+=tree[i].mass*tree[i].xcom;
    tree[0].ycom+=tree[i].mass*tree[i].ycom;
  }
  tree[0].xcom/=tree[0].mass;
  tree[0].ycom/=tree[0].mass;
}

void solve_dependencies(vector<Node>& tree,int id,double*x,double*y,double*m){
  if(tree[id].child_id<0){
    if(! tree[id].occupancy()) return;
    for(int i=tree[id].part_start;i<=tree[id].part_end;i++){
      tree[id].mass+=m[i];
      tree[id].xcom+=m[i]*x[i];
      tree[id].ycom+=m[i]*y[i];
    }
  }
  else{
    for(int i=0;i<4;i++)
#pragma omp task shared(tree) if(tree[id].level<PARALLEL_MAX_DEPTH)
            solve_dependencies(tree,tree[id].child_id,x,y,m);
    int child_id;
    for(int i=0;i<4;i++){
      child_id=tree[id].child_id+i;
      tree[id].mass+=tree[child_id].mass;
      tree[id].xcom+=tree[child_id].mass*tree[child_id].xcom;
      tree[id].ycom+=tree[child_id].mass*tree[child_id].ycom;
    }
  }
  tree[id].xcom/=tree[id].mass;
  tree[id].ycom/=tree[id].mass;
}


#endif //EX3_BUILD_H
