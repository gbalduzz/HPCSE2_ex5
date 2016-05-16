#pragma once
#include <assert.h>
using uint= unsigned int;

struct Node{
  int level;
  uint morton_id;
  int child_id;
  int part_start,part_end;
  float mass,xcom,ycom;
  Node():mass(0),xcom(0),ycom(0),child_id(-1),part_start(0),part_end(-6){}

  int occupancy()const {return std::max(part_end-part_start+1,0);}
  //Node child(int i);
};


/*Node& Node::child(int i) {
  assert(i<4 && i>=0);
  return child_id
}*/
