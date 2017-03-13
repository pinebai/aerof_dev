#pragma once


struct V6NodeData {
  int tet;
  int face;
  double r;
  double t;
  V6NodeData() { tet = -1; face = -1; r = 0.0; t = 0.0; }
};


