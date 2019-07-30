#ifndef FINITEDIFFERENCESTEMCIL
#define FINITEDIFFERENCESTEMCIL
#include <iostream>
#include <octave/oct.h>
using namespace std;

/* implemented relaxation schemes: Jacobi; Red Black Gauss Sidel; Zebra along the first dimension; Zebra along the second dimension; */
namespace Relaxation{
  enum iteration_type {J , RB, Z1, Z2};
}

namespace TriDiagonalSolver {
  static void solve(size_t N, double d[], double l[], double u[], double rhs[]);
}

/* an object representing a finite difference 3x3 stencil including boundary conditions */
class FiniteDifferenceStencil{
 private:
/* the stencil is held internally as a Matrix array */
  Matrix L[9];
/* apply the stencil to a function u at location (i,j) */
  double applyAt(Matrix u, octave_idx_type i, octave_idx_type j);
 public:
/* names of the directions each matrix represents */
  string dirString[9]; 

/* allocate internal data structure */
  FiniteDifferenceStencil(octave_idx_type m,octave_idx_type n);

/* fill the internal data structure */
  void setStencil(Matrix *S);
  void setStencil(FiniteDifferenceStencil S);

/* apply the stencil to a function u */
  void apply(Matrix u,Matrix* Lu);
  Matrix apply(Matrix u);
  Matrix operator*(Matrix u) { return apply(u);}

/* calculate the difference between two stencils */
  FiniteDifferenceStencil operator-(FiniteDifferenceStencil A);

/* perform relaxation for a number of iterations */
  void relax(Matrix *u0, Matrix f, int iterations, Relaxation::iteration_type flag);
  Matrix relax(Matrix u0, Matrix f, int iterations, Relaxation::iteration_type flag);
  Matrix relaxZebra1(Matrix u0, Matrix f, int iterations);
  void relaxZebra1(Matrix* u0, Matrix f, int iterations);
  Matrix relaxZebra2(Matrix u0, Matrix f, int iterations);
  void relaxZebra2(Matrix* u0, Matrix f, int iterations);
  Matrix relaxRedBlack(Matrix u0, Matrix f, int iterations);
  void relaxRedBlack(Matrix* u0, Matrix f, int iterations);
  Matrix relaxJacobi(Matrix u0, Matrix f, int iterations);
  void relaxJacobi(Matrix* u0, Matrix f, int iterations);

/* shortcut to access the internal data structure */
  Matrix operator[] (int k) {return L[k];}

/* print object */
  friend ostream& operator<< (ostream& o, FiniteDifferenceStencil& L);
};
#endif
