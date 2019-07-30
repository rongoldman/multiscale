#ifndef MULTIGRIDSETTING
#define MULTIGRIDSETTING
#include "FiniteDifferenceStencil.h"
#define PRECONDITION_ITERATIONS 1
#define REASONABLE_ITERATION_COUNT 8
#define MULTIGRID_PRECOND
class MultiGridSetting{
 typedef Relaxation::iteration_type rtype;
 static const rtype defaultRelaxation = Relaxation::RB;
 public:
  int maxIterations;
  double tolerance;
  int frontIterations;
  rtype frontRelaxation;
  int backIterations;
  rtype backRelaxation;
  int innerIterations;
  bool absNorm;
  bool verbose;

  MultiGridSetting(int fi, rtype fr, int bi, rtype br ,int ii, int maxi, double tol);
  MultiGridSetting(int fi, rtype fr, int bi, rtype br, int maxi, double tol);
  MultiGridSetting(int fi, int bi, int ii, int maxi, double tol);
  MultiGridSetting(int fi, int bi, int maxi, double tol);
  MultiGridSetting(int maxi, double tol);
  MultiGridSetting(int maxi);
  MultiGridSetting();
  MultiGridSetting(octave_scalar_map mgs);
/* print object */
  friend ostream& operator<< (ostream& o, MultiGridSetting& mgs);
};
static Relaxation::iteration_type parseRelaxation(std::string);
#endif
