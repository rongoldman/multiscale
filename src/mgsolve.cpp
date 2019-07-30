#include <time.h>
#include <iostream>
#include <octave/oct.h>
#include "MultiGridLevel.h"
#include "MultiGridSetting.h"
static bool is_valid_mgs(octave_scalar_map mgs);
static bool is_bad_args(const octave_value_list& args);
DEFUN_DLD(mgsolve,args,,
 "-*- texinfo -*-\n\
@deftypefn{Loadable Function} {@var{u} =} mgsolve(@var{L},@var{f},@var{u0})\n\
Perform multi-grid relaxation to solve the symmetric elliptic equation\n\
@group\n\
@var{L}@var{u} = @var{f}.\n\
@end group\n\
Where @var{u0} and @var{f} are MxN matrices representing the initial guess and non-homogenious RHS.\n\
@var{L} is a (M,N,9) representation of finite difference 3x3 stencil\n\
@infotex\n\
\n\
@example\n\
@group\n\
Lu(i,j) = L(i,j,9)*u(i-1,j+1)+L(i,j,2)*u(i,j+1)+L(i,j,3)*u(i+1,j+1)\n\
         + L(i,j,8)*u(i-1,j  )+L(i,j,1)*u(i,j  )+L(i,j,4)*u(i+1,j  )\n\
         + L(i,j,7)*u(i-1,j-1)+L(i,j,6)*u(i,j-1)+L(i,j,5)*u(i+1,j-1).\n\
@end group\n\
@end example\n\
\n\
@end infotex\n\
Note that u is zero-padded for convinience but boundary conditions should be provided by @var{L} and @var{f} \n\
@end deftypefn") 
{
 if (! is_bad_args(args)) {
  NDArray a(args(0).array_value());
  octave_idx_type m=a.dim1();
  octave_idx_type n=a.dim2();
  MultiGridLevel mg(m,n,0);
  MultiGridSetting* mgs;
  if (args.length()==4) {
   mgs = new MultiGridSetting(args(3).scalar_map_value());
  } else {
   mgs = new MultiGridSetting();
  }
  if (mgs->verbose) cout<< *mgs<<endl;
  Matrix* L = new Matrix[9]; 
  for (int k=0;k<9;k++) {
   L[k]=Matrix(m,n);
   for (octave_idx_type i=0;i<m;i++) 
    for (octave_idx_type j=0;j<n;j++)
     L[k].xelem(i,j) = a.xelem(i,j,k);
  }
  Matrix f(m+2,n+2,0.0);
  Matrix u0(m+2,n+2,0.0);
  for (octave_idx_type i=1;i<=m;i++) {
   for (octave_idx_type j=1;j<=n;j++) {
    f.xelem(i,j)=args(1).matrix_value().xelem(i-1,j-1);
    u0.xelem(i,j)=args(2).matrix_value().xelem(i-1,j-1);
   }
  }
  mg.setOperator(L);
  Matrix u;
  u = mg.relaxVCycle(u0,f,*mgs);
  double r=mg.remainder(*mgs);
  octave_value_list retval(2);
  retval(0)=octave_value(u.extract(1,1,m,n));
  retval(1)=octave_value(mg.status);
  retval(2)=octave_value(r);
  delete[] L;
  return retval; 
 }
 return octave_value(args(0)); 
}
static bool is_bad_args(const octave_value_list& args) {
 if (args.length() < 3 ) {
  error("Three arguments required.");
  return true;
 }
 if (args.length() > 4 ) {
  error("Too many arguments.");
  return true;
 }
 if (! args(0).is_real_matrix() ) {
  error("First argument not a matrix.");
  return true;
 } 
 if (args(0).array_value().ndims() !=3) {
  error("First argument not a 3D matrix.");
  return true;
 }  
 if (args(0).array_value().dim3()!=9) {
  error("First argument not a MxNx9 matrix.");
  return true;
 }
 if (! args(1).is_real_matrix() ) {
  error("Second argument not a matrix.");
  return true;
 } 
 if (args(1).array_value().ndims() !=2) {
  error("Second argument not a 2D matrix.");
  return true;
 }  
 if (args(1).array_value().dim1()!=args(0).array_value().dim1()
   && args(1).array_value().dim2()!=args(0).array_value().dim2()) {
  error("Second argument dimensions do not match the First argument dimensions.");
  return true;
 }
 if (! args(2).is_real_matrix() ) {
  error("Third argument not a matrix.");
  return true;
 } 
 if (args(2).array_value().ndims() !=2) {
  error("Third argument not a 2D matrix.");
  return true;
 }  
 if (args(2).array_value().dim1()!=args(0).array_value().dim1()
   && args(2).array_value().dim2()!=args(0).array_value().dim2()) {
  error("Third argument dimensions do not match the First argument dimensions.");
  return true;
 }
 if (args.length()==4) {
  if (! args(3).is_map()) {
   error("Fourth argument is no a structure.");
   return true;
  }
  if (! is_valid_mgs(args(3).scalar_map_value())) {
   error("Invalid multigrid setting structure.");
   return true;
  }
 }
 return false;
}
static bool is_valid_mgs(octave_scalar_map mgs){
 if (! mgs.contains("Tol")) {
  error("Missing parameter Tol");
  return false;
 }
 octave_value par = mgs.getfield("Tol");
 if(!(par.is_empty()||par.is_numeric_type())){
  error("Tol is not numeric type");
  return false;
 }
 if (! mgs.contains("FrontIterations")) {
  error("Missing parameter FrontIterations");
  return false;
 }
 par = mgs.getfield("FrontIterations");
 if(!(par.is_empty()||par.is_numeric_type())){
  error("FrontIterations must be an integer");
  return false;
 }
 if (! mgs.contains("BackIterations")) {
  error("Missing parameter BackIterations");
  return false;
 }
 par = mgs.getfield("BackIterations");
 if(!(par.is_empty()||par.is_numeric_type())){
  error("BackIterations must be an integer");
  return false;
 }
 if (! mgs.contains("InnerIterations")) {
  error("Missing parameter InnerIterations");
  return false;
 }
 par = mgs.getfield("InnerIterations");
 if(!(par.is_empty()||par.is_numeric_type())){
  error("InnerIterations must be an integer");
  return false;
 }
 if (! mgs.contains("MaxIterations")) {
  error("Missing parameter MaxIterations");
  return false;
 }
 par = mgs.getfield("MaxIterations");
 if(!(par.is_empty()||par.is_numeric_type())){
  error("MaxIterations must be an integer");
  return false;
 }
 if (! mgs.contains("FrontRelaxation")) {
  error("Missing parameter FrontRelaxation");
  return false;
 }
 par = mgs.getfield("FrontRelaxation");
 if(!(par.is_empty()||par.is_string())){
  error("FrontRelaxation must be a string");
  return false;
 }
 if (! mgs.contains("BackRelaxation")) {
  error("Missing parameter BackRelaxation");
  return false;
 }
 par = mgs.getfield("BackRelaxation");
 if(!(par.is_empty()||par.is_string())){
  error("BackRelaxation must be a string");
  return false;
 }
 if (! mgs.contains("NormControl")) {
  error("Missing parameter NormControl");
  return false;
 }
 par = mgs.getfield("NormControl");
 if(!(par.is_empty()||par.isinf().bool_value()||par.int_value()==2)){
  error("NormControl must be 2 or inf");
  return false;
 }
 if (! mgs.contains("Verbose")) {
  error("Missing parameter Verbose");
  return false;
 }
 par = mgs.getfield("Verbose");
 if(!(par.is_empty()||par.is_bool_type())){
  error("Verbose must be boolean");
  return false;
 }
 return true;
}
