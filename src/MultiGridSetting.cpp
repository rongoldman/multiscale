#include "MultiGridSetting.h"
MultiGridSetting::MultiGridSetting(int fi, rtype fr, int bi, rtype br, int ii, int maxi, double tol) { 
 maxIterations=maxi;
 tolerance=tol;
 frontIterations=fi;
 frontRelaxation=fr;
 backIterations=bi;
 backRelaxation=br;
 innerIterations=ii;
 absNorm=false;
 verbose=false;
}
MultiGridSetting::MultiGridSetting(int fi, rtype fr, int bi, rtype br, int maxi, double tol) { 
 maxIterations=maxi;
 tolerance=tol;
 frontIterations=fi;
 frontRelaxation=fr;
 backIterations=bi;
 backRelaxation=br;
 innerIterations=20;
 absNorm=false;
 verbose=false;
}
MultiGridSetting::MultiGridSetting(int fi, int bi, int ii, int maxi, double tol) {
 maxIterations=maxi;
 tolerance=tol;
 frontIterations=fi;
 frontRelaxation=defaultRelaxation;
 backIterations=bi;
 backRelaxation=defaultRelaxation;
 innerIterations=ii;
 absNorm=false;
 verbose=false;
}
MultiGridSetting::MultiGridSetting(int fi, int bi, int maxi, double tol) {
 maxIterations=maxi;
 tolerance=tol;
 frontIterations=fi;
 frontRelaxation=defaultRelaxation;
 backIterations=bi;
 backRelaxation=defaultRelaxation;
 innerIterations=20;
 absNorm=false;
 verbose=false;
}
MultiGridSetting::MultiGridSetting(int maxi, double tol) {
 maxIterations=maxi;
 tolerance=tol;
 frontIterations=1;
 frontRelaxation=defaultRelaxation;
 backIterations=1;
 backRelaxation=defaultRelaxation;
 innerIterations=20;
 absNorm=false;
 verbose=false;
}
MultiGridSetting::MultiGridSetting(int maxi) {
 maxIterations=maxi;
 tolerance=0.0;
 frontIterations=1;
 frontRelaxation=defaultRelaxation;
 backIterations=1;
 backRelaxation=defaultRelaxation;
 innerIterations=20;
 absNorm=false;
 verbose=false;
}
MultiGridSetting::MultiGridSetting(){
 maxIterations=200;
 tolerance=0.0;
 frontIterations=1;
 frontRelaxation=defaultRelaxation;
 backIterations=1;
 backRelaxation=defaultRelaxation;
 innerIterations=20;
 absNorm=false;
 verbose=false;
}
MultiGridSetting::MultiGridSetting(octave_scalar_map mgs){
 Array<int> val = mgs.getfield("MaxIterations").int_vector_value();
 maxIterations=(!val.isempty())?val(0):200;
 val = mgs.getfield("InnerIterations").int_vector_value();
 innerIterations=(!val.isempty())?val(0):20;
 val = mgs.getfield("FrontIterations").int_vector_value();
 frontIterations=(!val.isempty())?val(0):1;
 val = mgs.getfield("BackIterations").int_vector_value();
 backIterations=(!val.isempty())?val(0):1;
 Array<double> rval=mgs.getfield("Tol").vector_value();
 tolerance=(!rval.isempty())?rval(0):0.0;
 octave_value oval = mgs.getfield("NormControl");
 absNorm=(!oval.isempty())?oval.isinf().bool_value():false;
 oval = mgs.getfield("FrontRelaxation");
 frontRelaxation = (!oval.isempty())?parseRelaxation(oval.string_value()):defaultRelaxation;
 oval = mgs.getfield("BackRelaxation");
 backRelaxation = (!oval.isempty())?parseRelaxation(oval.string_value()):defaultRelaxation;
 oval = mgs.getfield("Verbose");
 verbose=(!oval.isempty())?oval.bool_value():false;
}
ostream& operator<< (ostream& o, MultiGridSetting& mgs) {
 o<<"Set to run V-cycle ("; 
 std::string fstr;
 std::string bstr;
 std::string norm = mgs.absNorm?"oo":"2";
 switch (mgs.frontRelaxation) {
 case Relaxation::J :
  fstr="J ";
  break;
 case Relaxation::RB :
  fstr="RB";
  break;
 case Relaxation::Z1 :
  fstr="Z1";
  break;
 case Relaxation::Z2 :
  fstr="Z2";
 }
 o<<fstr<<"^"<<mgs.frontIterations; 
 o<<" P V P^{-1} "; 
 switch (mgs.backRelaxation) {
 case Relaxation::J :
  bstr=" J";
  break;
 case Relaxation::RB :
  bstr="RB";
  break;
 case Relaxation::Z1 :
  bstr="Z1";
  break;
 case Relaxation::Z2 :
  bstr="Z2";
 }
 o<<bstr<<"^"<<mgs.backIterations;
 o<< "). Stop after "<<mgs.maxIterations <<" times";
 if (mgs.tolerance>0.0) {
  o<<"or if ||Lu-f||_"<<norm<<" < ||f||_"<<norm<<" "<<mgs.tolerance<<" .";
 }
 return o;
}
static Relaxation::iteration_type parseRelaxation(std::string rs) {
 if (rs=="J") {
  return Relaxation::J;
 }
 if (rs=="RB") {
  return Relaxation::RB;
 }
 if (rs=="Z1") {
  return Relaxation::Z1;
 }
 if (rs=="Z2") {
  return Relaxation::Z2;
 }
 error("Unrecognized relaxation scheme.");
 return Relaxation::RB;
}
