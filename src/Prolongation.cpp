#include <octave/oct.h>
#include "Prolongation.h"
#include "FiniteDifferenceStencil.h"

Prolongation::Prolongation(octave_idx_type mc, octave_idx_type nc, octave_idx_type mf, octave_idx_type nf) {
 this->mc=mc;
 this->nc=nc;
 this->mf=mf;
 this->nf=nf;
 for (int i=0;i<9;i++) P[i]=Matrix(mc+2,nc+2);
}
void Prolongation::reset(FiniteDifferenceStencil L){
 for (int k=0;k<9;k++) 
  P[k].fill(0.0);
 for (octave_idx_type i=1;i<=mc;i++) {
  for (octave_idx_type j=1;j<=nc;j++) {
   octave_idx_type i2=2*i;
   octave_idx_type j2=2*j;
   P[0].xelem(i,j)=1.0;
   P[1].xelem(i,j)=-(L[4].xelem(i2-1,j2)+L[5].xelem(i2-1,j2)+L[6].xelem(i2-1,j2))
             /(L[3].xelem(i2-1,j2)+L[0].xelem(i2-1,j2)+L[7].xelem(i2-1,j2));
   P[3].xelem(i,j)=-(L[6].xelem(i2,j2+1)+L[7].xelem(i2,j2+1)+L[8].xelem(i2,j2+1))
             /(L[5].xelem(i2,j2+1)+L[0].xelem(i2,j2+1)+L[1].xelem(i2,j2+1));
   P[5].xelem(i,j)=-(L[8].xelem(i2+1,j2)+L[2].xelem(i2+1,j2)+L[1].xelem(i2+1,j2))
             /(L[3].xelem(i2+1,j2)+L[0].xelem(i2+1,j2)+L[7].xelem(i2+1,j2));
   P[7].xelem(i,j)=-(L[2].xelem(i2,j2-1)+L[3].xelem(i2,j2-1)+L[4].xelem(i2,j2-1))
             /(L[5].xelem(i2,j2-1)+L[0].xelem(i2,j2-1)+L[1].xelem(i2,j2-1));
   P[2].xelem(i,j)=-(L[6].xelem(i2-1,j2+1)+P[1].xelem(i,j)*L[7].xelem(i2-1,j2+1)+P[3].xelem(i,j)*L[5].xelem(i2-1,j2+1))
             /(L[0].xelem(i2-1,j2+1));
   P[4].xelem(i,j)=-(L[8].xelem(i2+1,j2+1)+P[3].xelem(i,j)*L[1].xelem(i2+1,j2+1)+P[5].xelem(i,j)*L[7].xelem(i2+1,j2+1))
             /(L[0].xelem(i2+1,j2+1));
   P[6].xelem(i,j)=-(L[2].xelem(i2+1,j2-1)+P[5].xelem(i,j)*L[3].xelem(i2+1,j2-1)+P[7].xelem(i,j)*L[1].xelem(i2+1,j2-1))
             /(L[0].xelem(i2+1,j2-1));
   P[8].xelem(i,j)=-(L[4].xelem(i2-1,j2-1)+P[1].xelem(i,j)*L[3].xelem(i2-1,j2-1)+P[7].xelem(i,j)*L[5].xelem(i2-1,j2-1))
             /(L[0].xelem(i2-1,j2-1));
  }
 }
 
 // Refrain from interpolating to pads when m and/or n are even.
 if (mf % 2 == 0) {
  for (octave_idx_type j=1;j<=nc;j++) {
   P[4].xelem(mc,j) = 0;
   P[5].xelem(mc,j) = 0;
   P[6].xelem(mc,j) = 0;
  }
 }
 if (nf % 2 == 0) {
  for (octave_idx_type i=1;i<=mc;i++) {
    P[2].xelem(i,nc) = 0;
    P[3].xelem(i,nc) = 0;
    P[4].xelem(i,nc) = 0;
  }
 }

}
void Prolongation::coarsen(Matrix uf,Matrix* uc) {
 for (octave_idx_type i=1;i<=mc;i++) {
  for (octave_idx_type j=1;j<=nc;j++) {
   octave_idx_type i2=2*i;
   octave_idx_type j2=2*j;
   uc->xelem(i,j)=
      P[0].xelem(i,j)*uf.xelem(i2  ,j2  )
     +P[1].xelem(i,j)*uf.xelem(i2-1,j2  )
     +P[2].xelem(i,j)*uf.xelem(i2-1,j2+1)
     +P[3].xelem(i,j)*uf.xelem(i2  ,j2+1)
     +P[4].xelem(i,j)*uf.xelem(i2+1,j2+1)
     +P[5].xelem(i,j)*uf.xelem(i2+1,j2  )
     +P[6].xelem(i,j)*uf.xelem(i2+1,j2-1)
     +P[7].xelem(i,j)*uf.xelem(i2  ,j2-1)
     +P[8].xelem(i,j)*uf.xelem(i2-1,j2-1);
  }
 }
}
Matrix Prolongation::coarsen(Matrix uf) {
 Matrix uc(mc+2,nc+2,0.0);
 coarsen(uf,&uc);
 return uc;
}
void Prolongation::refine(Matrix uc, Matrix *uf) {
 for (octave_idx_type i=1;i<=mc;i++) {
  for (octave_idx_type j=1;j<=nc;j++) {
   octave_idx_type i2=2*i;
   octave_idx_type j2=2*j;
   uf->xelem(i2  ,j2  )+=P[0].xelem(i,j)*uc.xelem(i,j);
   uf->xelem(i2-1,j2  )+=P[1].xelem(i,j)*uc.xelem(i,j);
   uf->xelem(i2-1,j2+1)+=P[2].xelem(i,j)*uc.xelem(i,j);
   uf->xelem(i2  ,j2+1)+=P[3].xelem(i,j)*uc.xelem(i,j);
   uf->xelem(i2+1,j2+1)+=P[4].xelem(i,j)*uc.xelem(i,j);
   uf->xelem(i2+1,j2  )+=P[5].xelem(i,j)*uc.xelem(i,j);
   uf->xelem(i2+1,j2-1)+=P[6].xelem(i,j)*uc.xelem(i,j);
   uf->xelem(i2  ,j2-1)+=P[7].xelem(i,j)*uc.xelem(i,j);
   uf->xelem(i2-1,j2-1)+=P[8].xelem(i,j)*uc.xelem(i,j);
  }
 }
}
Matrix Prolongation::refine(Matrix uc) {
 Matrix uf(mf+2,nf+2,0.0);
 refine(uc,&uf);
 return uf;
} 
void Prolongation::Galerkin(FiniteDifferenceStencil Lf,FiniteDifferenceStencil Lc) {
 Lc[0].fill(1.0);
 for (int k=1;k<9;k++) 
  Lc[k].fill(0.0);

/* 
  We have to do this nine times to separate the contributions. We put 1 at
  some coarse-grid variables, prolong, apply L on the fine grid, and restrict.
  The coarse-grid values thus obtained are the required coefficients. 
*/
 octave_idx_type i0[9]={2,3,3,2,1,1,1,2,3};
 octave_idx_type j0[9]={2,2,1,1,1,2,3,3,3};
 Matrix uc(mc+2,nc+2);
 Matrix uf(mf+2,nf+2,0.0);
 Matrix Luf(mf+2,nf+2,0.0);
 for (int k=0;k<9;k++) {
  uc.fill(0.0);
  uf.fill(0.0);
  Luf.fill(0.0);
  
  for (octave_idx_type i=i0[k];i<=mc;i+=3) 
   for (octave_idx_type j=j0[k];j<=nc;j+=3) 
    uc.xelem(i,j)=1.0;
  refine(uc,&uf);
  Lf.apply(uf,&Luf);
  uc.fill(0.0);
  coarsen(Luf,&uc);
  for (int l=0;l<9;l++)
   for (octave_idx_type i=i0[k]+i0[l]-2;i<=mc;i+=3) 
    for (octave_idx_type j=j0[k]+j0[l]-2;j<=nc;j+=3) 
     Lc[l].xelem(i,j)=uc.xelem(i,j);

 }
 return;

}
ostream& operator<< (ostream& o, Prolongation& P) {
 for (int k=0;k<9;k++) o<<"k="<<k<<endl<<P[k]<<endl;
}
