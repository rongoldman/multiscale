#ifndef MULTIGRIDLEVEL
#define MULTIGRIDLEVEL
#define EPS std::numeric_limits<double>::epsilon()
#include <octave/oct.h>
#include "FiniteDifferenceStencil.h"
#include "Prolongation.h"
#include "MultiGridSetting.h"
class MultiGridLevel{

// members :

/* dimesions of the grid level */
int n,m;

/* coarsening level (i.e. level=0 is the finest grid, level = 1 is twice as coarse etc'...) */
int level;

/* status of relaxVCycle or relaxAsPrecondition methods. 
   If positive gives the number of iterations performed.
   If status = -1 then the number of iterations was to small to achieve the specified convergence.
   If status = -2 the scheme diverges or the specified tolerence is too small w.r.t. roundoff.*/
public: int status; 

/* prolongation operator to translate function and stencils between this level and the coarser level (NULL if this is the coarsest level) */
public: Prolongation *P;

/* Finite difference elliptic operator prolonged to this level */
public: FiniteDifferenceStencil *L;

/* cached functions for inverting L with initial guess u and rhs f at this level */
public: Matrix *u,*f;

/* the coarser level (NULL if this is the coarsest level) */
public: MultiGridLevel *coarse;

// methods :

/* constructor allocating room for this level coarser ones */
MultiGridLevel(int m, int n, int lev);

/* destructor for this level and coarser ones */
~MultiGridLevel();

/* copy the contents of L to the FiniteDifferenceStencil of this level and then prolong it to coarser levels */
public: void setOperator(Matrix* L);
public: void setOperator(FiniteDifferenceStencil L);

/* prolong the FiniteDifferenceStencil of this level to coarser levels */
private: void calcLevels();

/* perform V-cycle multigrid relaxation. If initial guess u0 and right hand side f are not specified the cached values are used */
public: Matrix relaxVCycle(Matrix u0, Matrix f, MultiGridSetting opt);
public: Matrix relaxVCycle(MultiGridSetting opt);

/* Use the existing MultiGridLevel as precondition to invert a slightly different elliptic operator A. 
   This solves iteratively the problem Lu=f-(A-L)u via V-cycle relaxation of L (if |(A-L)| is small enough ). */
public: Matrix relaxAsPrecondition(Matrix u0, Matrix f, FiniteDifferenceStencil A, MultiGridSetting opt);

/* set the initial guess of the V-cycle relaxation */
public: void setIC(double val) {u->fill(val);}
public: void setIC(Matrix u0);

/* set the right hand side of the V-cycle relaxation */
public: void setRHS(Matrix rhs);

/* calculate |Lu-f| */
public: double remainder(MultiGridSetting opt);
};
double maxNorm(Matrix a);
#endif
