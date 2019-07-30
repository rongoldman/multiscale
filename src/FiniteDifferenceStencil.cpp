#include <octave/oct.h>
#include "FiniteDifferenceStencil.h"
FiniteDifferenceStencil::FiniteDifferenceStencil(octave_idx_type m,octave_idx_type n) {
 for (int k=0;k<9;k++) L[k]=Matrix(m+2,n+2,(k==0?1.0:0.0));
 dirString[0] = "Center";
 dirString[1] = "Up";
 dirString[2] = "Up-Left";
 dirString[3] = "Left";
 dirString[4] = "Down-Left";
 dirString[5] = "Down";
 dirString[6] = "Down-Right";
 dirString[7] = "Right";
 dirString[8] = "Up-Right";

}
void FiniteDifferenceStencil::setStencil(Matrix* S) {
 for (int k=0;k<9;k++) {
  L[k].fill(k==0?1.0:0.0);
  octave_idx_type m = L[k].dim1()-2;
  octave_idx_type n = L[k].dim2()-2;
  if (S[k].dim1()==1 && S[k].dim2()==1) {
   L[k].fill(S[k](1,1),1,1,m,n);
  } else if(S[k].dim1()==m && S[k].dim2()==n) {
    L[k].insert(S[k],1,1);
  }else {
   //gripe_nonconformant ("setStencil", m, n, S[k].dim1(), S[k].dim2());
   octave::err_nonconformant ("setStencil", m, n, S[k].dim1(), S[k].dim2());
  }
 }
 /*{ octave_idx_type m = L[0].dim1()-2;
  octave_idx_type n = L[0].dim2()-2;
  for (octave_idx_type i=1;i<=m;i++) {
   for (octave_idx_type j=1;j<=n;j++){
    if (abs(L[0].xelem(i,j)<1e-9)) cout<< "L0("<<i <<","<< j<<")="<< L[0].xelem(i,j)<<endl;
   }
  }}*/
}
void FiniteDifferenceStencil::setStencil(FiniteDifferenceStencil S) {
 for (int k=0;k<9;k++) {
  L[k].fill(k==0?1.0:0.0);
  octave_idx_type m = L[k].dim1()-2;
  octave_idx_type n = L[k].dim2()-2;
  if(S[k].dim1()-2 ==m && S[k].dim2()-2 ==n) {
    L[k].insert(S[k],0,0);
  } else {
   //gripe_nonconformant ("setStencil", m, n, S[k].dim1()-2, S[k].dim2()-2);
   octave::err_nonconformant ("setStencil", m, n, S[k].dim1()-2, S[k].dim2()-2);
  }
 }
}
FiniteDifferenceStencil FiniteDifferenceStencil::operator-(FiniteDifferenceStencil A) {
 octave_idx_type n1=L[0].dim1()-2;
 octave_idx_type n2=L[0].dim2()-2;
 octave_idx_type a_n1=A[0].dim1()-2;
 octave_idx_type a_n2=A[0].dim2()-2;
 //if (n1!=a_n1 || n2 != a_n2) gripe_nonconformant ("operator -", n1, n2, a_n1, a_n2);
 if (n1!=a_n1 || n2 != a_n2) octave::err_nonconformant ("operator -", n1, n2, a_n1, a_n2);
 Matrix d[9];
 FiniteDifferenceStencil diff(n1,n2); 
 for (int k=0;k<9;k++) {
  d[k]=(L[k]-A[k]).extract(1,1,n1,n2);
 }
 diff.setStencil(d);
 return diff;
} 
void FiniteDifferenceStencil::apply(Matrix u,Matrix *Lu){
  octave_idx_type m=L[0].dim1()-2;
  octave_idx_type n=L[0].dim2()-2;
  for (octave_idx_type i=1;i<=m;i++) {
    for (octave_idx_type j=1;j<=n;j++){
      Lu->xelem(i,j)=applyAt(u,i,j);
    }
  }
}
Matrix FiniteDifferenceStencil::apply(Matrix u){
  octave_idx_type m=L[0].dim1()-2;
  octave_idx_type n=L[0].dim2()-2;
  Matrix Lu(m+2,n+2,0.0);
  apply(u,&Lu);
  return Lu;  
}
double FiniteDifferenceStencil::applyAt(Matrix u, octave_idx_type i, octave_idx_type j) {
 return L[0].xelem(i,j)*u.xelem(i  ,j  )
       +L[1].xelem(i,j)*u.xelem(i-1,j  )
       +L[2].xelem(i,j)*u.xelem(i-1,j+1)
       +L[3].xelem(i,j)*u.xelem(i  ,j+1)
       +L[4].xelem(i,j)*u.xelem(i+1,j+1)
       +L[5].xelem(i,j)*u.xelem(i+1,j  )
       +L[6].xelem(i,j)*u.xelem(i+1,j-1)
       +L[7].xelem(i,j)*u.xelem(i  ,j-1)
       +L[8].xelem(i,j)*u.xelem(i-1,j-1);
}
void FiniteDifferenceStencil::relax(Matrix *u0, Matrix f, int iterations, Relaxation::iteration_type flag) {
 switch (flag) {
  case Relaxation::Z1:
   relaxZebra1(u0, f, iterations);
   break;
  case Relaxation::Z2:
   relaxZebra2(u0, f, iterations);
   break;
  case Relaxation::RB:
   relaxRedBlack(u0, f, iterations);
   break;
  case Relaxation::J:
   relaxJacobi(u0, f, iterations);
   break;
 }
}
Matrix FiniteDifferenceStencil::relax(Matrix u0, Matrix f, int iterations, Relaxation::iteration_type flag) {
 switch (flag) {
  case Relaxation::Z1:
   return relaxZebra1(u0, f, iterations);
  case Relaxation::Z2:
   return relaxZebra2(u0, f, iterations);
  case Relaxation::RB:
   return relaxRedBlack(u0, f, iterations);
  case Relaxation::J:
   return relaxJacobi(u0, f, iterations);
 }
 cout<<"WTF "<<flag<<" not in "<<Relaxation::J<<","<<Relaxation::RB<<endl;
}
void FiniteDifferenceStencil::relaxRedBlack(Matrix *u, Matrix f,int it){
 octave_idx_type i0[4]={1,2,1,2};
 octave_idx_type j0[4]={1,2,2,1};
 octave_idx_type m=L[0].dim1()-2;
 octave_idx_type n=L[0].dim2()-2;
 Matrix Lu(m+2,n+2,0.0);
 for (int l=0;l<it;l++) {
  for (int k=0;k<2;k++) {
   for (octave_idx_type i=i0[k];i<=m;i+=2) {
    for (octave_idx_type j=j0[k];j<=n;j+=2) {
     Lu.xelem(i,j)=applyAt(*u,i,j);
    }
   }
  }
  for (int k=0;k<2;k++) {
   for (octave_idx_type i=i0[k];i<=m;i+=2) {
    for (octave_idx_type j=j0[k];j<=n;j+=2) {
      u->xelem(i,j)+=(f.xelem(i,j)-Lu.xelem(i,j))/L[0].xelem(i,j);
    }
   }
  }
  for (int k=2;k<4;k++) {
   for (octave_idx_type i=i0[k];i<=m;i+=2) {
    for (octave_idx_type j=j0[k];j<=n;j+=2) {
     Lu.xelem(i,j)=applyAt(*u,i,j);
    }
   }
  }
  for (int k=2;k<4;k++) {
   for (octave_idx_type i=i0[k];i<=m;i+=2) {
    for (octave_idx_type j=j0[k];j<=n;j+=2) {
      u->xelem(i,j)+=(f.xelem(i,j)-Lu.xelem(i,j))/L[0].xelem(i,j);
    }
   }
  }
 }
}
Matrix FiniteDifferenceStencil::relaxRedBlack(Matrix u, Matrix f,int it){
 Matrix* v = new Matrix(u);
 relaxRedBlack(v,f,it);
 return *v;
}
void FiniteDifferenceStencil::relaxJacobi(Matrix *u, Matrix f,int it){
 octave_idx_type m=L[0].dim1()-2;
 octave_idx_type n=L[0].dim2()-2;
 Matrix Lu(m+2,n+2,0.0);
 for (int l=0;l<it;l++) {
   apply(*u,&Lu);
   *u = *u + quotient(f-Lu , L[0]);
 }
}
Matrix FiniteDifferenceStencil::relaxJacobi(Matrix u, Matrix f,int it){
 Matrix v=Matrix(u);
 for (int l=0;l<it;l++) {
   v = v + quotient(f-apply(v) , L[0]);
 }
 return v;
}
Matrix FiniteDifferenceStencil::relaxZebra1(Matrix u, Matrix f,int it){
 Matrix* v = new Matrix(u);
 relaxZebra1(v,f,it);
 return *v;
}
Matrix FiniteDifferenceStencil::relaxZebra2(Matrix u, Matrix f,int it){
 Matrix* v = new Matrix(u);
 relaxZebra2(v,f,it);
 return *v;
}
void FiniteDifferenceStencil::relaxZebra1(Matrix *u, Matrix f,int it){
 octave_idx_type m=L[0].dim1()-2;
 octave_idx_type n=L[0].dim2()-2;
 for (octave_idx_type stripes = 1;stripes<3;stripes++) {
    for (octave_idx_type i = stripes;i<=m;i+=2){
        // Compute the residual for the entire line i
      double rhs[n];
      double mdiag[n];
      double udiag[n];
      double ldiag[n];
//      SparseMatrix trid(n,n);
      for (octave_idx_type j = 1; j<=n; j++){
        rhs[j-1] = f.xelem(i,j) - applyAt(*u,i,j);
        ldiag[j-1]=(j>1) ? L[7].xelem(i,j) : 0.0;
        mdiag[j-1]=L[0].xelem(i,j);
        udiag[j-1]=(j<n) ? L[3].xelem(i,j) : 0.0;
      }
      TriDiagonalSolver::solve(n,mdiag,ldiag,udiag,rhs);
      for (octave_idx_type j=1;j<=n;j++) {
        u->xelem(i,j)+=rhs[j-1];
      }
    }
 }
}
void FiniteDifferenceStencil::relaxZebra2(Matrix *u, Matrix f,int it){
 octave_idx_type m=L[0].dim1()-2;
 octave_idx_type n=L[0].dim2()-2;
 for (octave_idx_type stripes = 1;stripes<3;stripes++) {
    for (octave_idx_type j = stripes;j<=n;j+=2){
        // Compute the residual for the entire line i
      double rhs[m];
      double mdiag[m];
      double udiag[m];
      double ldiag[m];
//      ColumnVector rhs(m);
//      SparseMatrix trid(m,m);
      for (octave_idx_type i = 1; i<=m; i++){
//        rhs.xelem(i-1) = f.xelem(i,j) - applyAt(*u,i,j);
        rhs[i-1] = f.xelem(i,j) - applyAt(*u,i,j);
        ldiag[i-1]=(i>1) ? L[1].xelem(i,j) : 0.0;
        mdiag[i-1]=L[0].xelem(i,j);
        udiag[i-1]=(i<m) ? L[5].xelem(i,j) : 0.0;
//        if (i>1) {
//         trid.assign(i-1,i-2,Sparse<double>(1,1,L[1].xelem(i,j)));
//        }
//        trid.assign(i-1,i-1,Sparse<double>(1,1,L[0].xelem(i,j)));
//        if (i<m) {
//         trid.assign(i-1,i,Sparse<double>(1,1,L[5].xelem(i,j)));
//        }
      }
//      ColumnVector du=trid.solve(rhs);
      TriDiagonalSolver::solve(m,mdiag,ldiag,udiag,rhs);
      for (octave_idx_type i=1;i<=m;i++) {
//        u->xelem(i,j)+=du.xelem(i-1);
        u->xelem(i,j)+=rhs[i-1];
      }
    }
 }
}
ostream& operator<< (ostream& o, FiniteDifferenceStencil& L) {
 for (int k=0;k<9;k++) o<<L.dirString[k]<<endl<<L[k];
 return o;
}
void TriDiagonalSolver::solve(size_t N, double d[], double l[], double u[], double x[]) {
//    cout<< "begin"<<endl<<"l d u rhs"<<endl;
//    for (size_t i=0;i<N;i++) cout<<l[i]<<" "<<d[i]<<" "<<u[i] <<" "<<x[i]<<endl;
//    cout << endl;
    u[0] = u[0] / d[0];
    x[0] = x[0] / d[0];
   
    /* loop from 1 to N - 1 inclusive */
    for (size_t i = 1; i < N; i++) {
        double det = (d[i] - l[i] * u[i - 1]);
        u[i] = u[i] / det;
        x[i] = (x[i] - l[i] * x[i - 1]) / det;
    }
   
    /* loop from N - 2 to 0 inclusive, safely testing loop end condition */
    for (size_t i = N - 1; i-- > 0; )
        x[i] = x[i] - u[i] * x[i + 1];

//    for (size_t i=0;i<N;i++) cout<<x[i]<<" ";
//    cout<<endl;
}
