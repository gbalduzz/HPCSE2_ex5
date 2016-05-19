#pragma once
using uint= unsigned int;

struct Node
{
  int part_start, part_end;
  int child_id;
  int level;
  uint morton_id;
  double r2; //square radius from com
  double mass=0; //sum of abs(w)
  double xcom, ycom;

  void setup(int part_start, int part_end, int level, int morton_id)
  {
    this->part_start = part_start;
    this->part_end = part_end;
    this->child_id = 0;
    this->level = level;
    this->morton_id = morton_id;
  }

  inline void setCom(double m,double x ,double y){
    mass=m;
    xcom=x;
    ycom=y;
  }
};
