#ifndef PROLONGATION
#define PROLONGATION
#include <iostream>
#include <octave/oct.h>
#include "FiniteDifferenceStencil.h"
class Prolongation {
 Matrix P[9];
 octave_idx_type mc,nc;
 octave_idx_type mf,nf;
 public:
  Prolongation(octave_idx_type mc, octave_idx_type nc, octave_idx_type mf, octave_idx_type nf);
  void reset(FiniteDifferenceStencil L);
  Matrix refine(Matrix uc);
  void refine(Matrix uc, Matrix* uf);
  Matrix coarsen(Matrix uf);
  void coarsen(Matrix uf, Matrix* uc);
  void Galerkin(FiniteDifferenceStencil Lf,FiniteDifferenceStencil Lc);
  FiniteDifferenceStencil Galerkin(FiniteDifferenceStencil Lf);
  octave_idx_type getCoarseM() {return mc;}
  octave_idx_type getCoarseN() {return nc;}
  Matrix operator[] (int k) {return P[k];}
  friend ostream& operator<< (ostream& o, Prolongation& P);
};
#endif
