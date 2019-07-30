#include <octave/oct.h>
#include "MultiGridLevel.h"
MultiGridLevel::MultiGridLevel(int m,int n,int lev) {
  this->m=m;
  this->n=n;
  this->level=lev;
  L=new FiniteDifferenceStencil(m,n);
  u=new Matrix(m+2,n+2,0.0);
  f=new Matrix(m+2,n+2,0.0);
  int mc=m/2;
  int nc=n/2;
  if ((mc>1) && (nc>1)){ // recursive allocation
    P=new Prolongation(mc,nc,m,n);
    coarse=new MultiGridLevel(mc,nc,lev+1);
  } else {
    P=NULL;
    coarse=NULL;
  }
}
MultiGridLevel::~MultiGridLevel() {
 if (coarse!=NULL) delete coarse;
 if (P!=NULL) delete P;
 if (L!=NULL) delete L;
 if (u!=NULL) delete u;
 if (f!=NULL) delete f;
}
void MultiGridLevel::setOperator(Matrix* Op) {
 L->setStencil(Op);
 calcLevels(); 
}
void MultiGridLevel::setOperator(FiniteDifferenceStencil Op) {
 L->setStencil(Op);
 calcLevels(); 
}
void MultiGridLevel::setIC(Matrix u0) {
 if( u==NULL || u==&u0 ) {
  u=&u0;
  return;
 }
 if(u0.dim1() != u->dim1() && u0.dim2() != u->dim2()) {
  //gripe_nonconformant ("operator =", u0.dim1()-2, u0.dim2()-2, u->dim1()-2, u->dim2()-2);
  octave::err_nonconformant ("operator =", u0.dim1()-2, u0.dim2()-2, u->dim1()-2, u->dim2()-2);
 }
 delete u;
 u=new Matrix(u0);
}
void MultiGridLevel::setRHS(Matrix rhs) {
 if( f==NULL || f==&rhs ) {
  f=&rhs;
  return;
 }
 if(rhs.dim1() != f->dim1() && rhs.dim2() != f->dim2()) {
  //gripe_nonconformant ("operator =", rhs.dim1()-2, rhs.dim2()-2, f->dim1()-2, f->dim2()-2);
  octave::err_nonconformant ("operator =", rhs.dim1()-2, rhs.dim2()-2, f->dim1()-2, f->dim2()-2);
 }
 delete f;
 f=new Matrix(rhs);
}
void MultiGridLevel::calcLevels(){
 if (P!=NULL) {
  P->reset(*L);
  P->Galerkin(*L,*(coarse->L));
  coarse->calcLevels();
 } 
}
Matrix MultiGridLevel::relaxVCycle(Matrix initial, Matrix rhs, MultiGridSetting opt) { 
 setIC(initial);
 setRHS(rhs);
 return relaxVCycle(opt);
}
Matrix MultiGridLevel::relaxVCycle(MultiGridSetting opt) {
 status=0;
 if (P==NULL || coarse==NULL) {
  L->relax(u,*f,opt.innerIterations/2,opt.frontRelaxation);
  L->relax(u,*f,opt.innerIterations/2,opt.backRelaxation);
 } else {
  Matrix Lu(m+2,n+2,0.0);
  double fnorm;
  if (level==0 && opt.tolerance>0.0) {
   if (!opt.absNorm) {
    fnorm=sqrt(((f->sumsq()).sum())(0));
   } else {
    fnorm=maxNorm(*f);
   }
  }
  int itmax=(level==0)?opt.maxIterations:1;
  for (int it=0;it<itmax;it++) {
   if (opt.tolerance>0.0) {
    if (opt.verbose) {
     if (!opt.absNorm) {
      if (level==0) cout<<"( "<<it<<") |f-Lu|_2="<<remainder(opt)<<" ? "<<opt.tolerance*fnorm<<endl;
     } else {
      if (level==0) cout<<"( "<<it<<") |f-Lu|_oo="<<remainder(opt)<<" ? "<<opt.tolerance*fnorm<<endl;
     }
    }
    if (level==0 && (remainder(opt)<=max(opt.tolerance*fnorm,10*EPS) ) ) break;
   } else {
    if (opt.verbose) {
     if (!opt.absNorm) {
      cout<<"|f-Lu|_2="<<remainder(opt)<<endl;
     } else {
      cout<<"|f-Lu|_oo="<<remainder(opt)<<endl;
     }
    }
   }
   coarse->u->fill(0.0);
   L->relax(u,*f,opt.frontIterations,opt.frontRelaxation);
   L->apply(*u,&Lu);
   P->coarsen(*f-Lu,coarse->f);
   *u=*u+P->refine(coarse->relaxVCycle(opt));
   L->relax(u,*f,opt.backIterations,opt.backRelaxation);
   status=it+1;
  }
  if (level==0 && opt.tolerance>0.0 && (remainder(opt)>max(opt.tolerance*fnorm,10*EPS) ) ) status=-1;;
 }
 return *u;
}
Matrix MultiGridLevel::relaxAsPrecondition(Matrix initial, Matrix rhs, FiniteDifferenceStencil op, MultiGridSetting opt) {
 status=0;
 if (level != 0) cout<<"TODO: add exception here for level ne 0"<<endl;
 double f0norm;
 double rem;
 double oldrem;
 setIC(initial);
 FiniteDifferenceStencil dL=op - *L;
 // V-cycle options: one iteration, don't check tolerance
 MultiGridSetting vopt=opt;
 vopt.maxIterations=PRECONDITION_ITERATIONS;
 vopt.tolerance=0.0;
 if (opt.tolerance>0.0) {
  if (!opt.absNorm) {
   f0norm=sqrt(((rhs.sumsq()).sum())(0));
  } else {
   f0norm=maxNorm(rhs);
  }
 }
 for (int it=0;it<=opt.maxIterations;it++) {
  setRHS(rhs-dL*(*u));
  if (opt.tolerance>0.0) {
   oldrem=rem;
   if (!opt.absNorm) {
    rem=sqrt( ( ((rhs-op*(*u)).sumsq()).sum() )(0) );
//    cout<<"|f-Lu|_2="<<rem<<endl;
   }else{
    rem=maxNorm(rhs-op*(*u));
//    cout<<"|f-Lu|_oo="<<rem<<endl;
   }
   if (rem<=max(opt.tolerance*f0norm,10*EPS) ) break;
   if (it==opt.maxIterations) {
    cout<<"WARNING: Tolerence level not reached."<<endl;
    status=-1;
    break;
   }
   if (rem>oldrem && it!=0) {
    cout<<"WARNING: Iteration is not converging anymore."<<endl
        <<"         Either roundoff is greater than tolerence level or preconditioning failed."<<endl;
    status=-2;
    break;
   }
  }
  relaxVCycle(vopt);
  status=it+1;
 }
 return *u;
}
double MultiGridLevel::remainder(MultiGridSetting opt) {
 Matrix r(m+2,n+2,0.0);
 L->apply(*u,&r);
 if (!opt.absNorm) {
  double nr=0.0;
  for (octave_idx_type i=1;i<=m;i++){
   double lsum=0.0;
   for (octave_idx_type j=1;j<=n;j++)
    lsum+=pow(f->xelem(i,j)-r.xelem(i,j),2);
   nr+=lsum;
  }
  return sqrt(nr);
 } else {
  double nr=0.0;
  for (octave_idx_type i=1;i<=m;i++){
   for (octave_idx_type j=1;j<=n;j++){
     nr=max(nr,abs(f->xelem(i,j)-r.xelem(i,j)));
   }
  }
  return nr;
 }
}
double maxNorm( Matrix a) {
 double nrm=0;
 octave_idx_type m = a.dim1();
 octave_idx_type n = a.dim2();
 for (octave_idx_type i=0;i<m;i++){
  for (octave_idx_type j=0;j<n;j++){
    nrm=max(nrm,abs(a.xelem(i,j)));
  }
 }
 return nrm; 
}
