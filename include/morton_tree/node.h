#pragma once
using uint= unsigned int;

struct Node
{
  int part_start, part_end;
  int child_id;
  double mass, xcom, ycom;
  int level, morton_id;

  void setup(int part_start, int part_end, int level, int morton_id)
  {
    this->part_start = part_start;
    this->part_end = part_end;
    this->child_id = 0;
    this->level = level;
    this->morton_id = morton_id;
  }
};
