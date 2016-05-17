#pragma once
#include <assert.h>
using uint= unsigned int;

template<int k>
struct Node{
  double xcom,ycom,mass;
  double r2; //squared radius from com
  int level;
  int child_id;
  int part_start,part_end;
  uint morton_id;
  double re_expansion[k+1];
  double im_expansion[k+1];

  Node():xcom(0),ycom(0),mass(0),child_id(-1),part_start(0),part_end(-5){}
  int occupancy()const {return std::max(part_end-part_start+1,0);}
};

#include <vector>
template<int k>
using Tree = std::vector<Node<k>>;
