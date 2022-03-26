// solve boundary value problem by relaxation method
// W. H. Press, et al, "Numerical Recipes" section 17.3

#include <iostream>
#include <iomanip>
#include <cmath>
#include "nr.h"
#include "Mat3D.h"
using namespace std;

void bksub(const int ne, const int nb, const int jf, const int k1,
           const int k2, Mat3D_IO_DP &c)
{
    int nbf,im,kp,k,j,i;
    DP xx;
    
    nbf=ne-nb;
    im=1;
    for (k=k2-1;k>=k1;k--) {
        if (k == k1) im=nbf+1;
        kp=k+1;
        for (j=0;j<nbf;j++) {
            xx=c[j][jf][kp];
            for (i=im-1;i<ne;i++)
                c[i][jf][k] -= c[i][j][k]*xx;
        }
    }
    for (k=k1;k<k2;k++) {
        kp=k+1;
        for (i=0;i<nb;i++) c[i][0][k]=c[i+nbf][jf][k];
        for (i=0;i<nbf;i++) c[i+nb][0][k]=c[i][jf][kp];
    }
}

void red(const int iz1, const int iz2, const int jz1, const int jz2,
         const int jm1, const int jm2, const int jmf, const int ic1,
         const int jc1, const int jcf, const int kc, Mat3D_I_DP &c,
         Mat_IO_DP &s)
{
    int loff,l,j,ic,i;
    DP vx;
    
    loff=jc1-jm1;
    ic=ic1;
    for (j=jz1;j<jz2;j++) {
        for (l=jm1;l<jm2;l++) {
            vx=c[ic][l+loff][kc-1];
            for (i=iz1;i<iz2;i++) s[i][l] -= s[i][j]*vx;
        }
        vx=c[ic][jcf][kc-1];
        for (i=iz1;i<iz2;i++) s[i][jmf] -= s[i][j]*vx;
        ic += 1;
    }
}

void pinvs(const int ie1, const int ie2, const int je1, const int jsf,
           const int jc1, const int k, Mat3D_O_DP &c, Mat_IO_DP &s)
{
    int jpiv,jp,je2,jcoff,j,irow,ipiv,id,icoff,i;
    DP pivinv,piv,dum,big;
    
    const int iesize=ie2-ie1;
    Vec_INT indxr(iesize);
    Vec_DP pscl(iesize);
    je2=je1+iesize;
    for (i=ie1;i<ie2;i++) {
        big=0.0;
        for (j=je1;j<je2;j++)
            if (fabs(s[i][j]) > big) big=fabs(s[i][j]);
        if (big == 0.0)
            nrerror("Singular matrix - row all 0, in pinvs");
        pscl[i-ie1]=1.0/big;
        indxr[i-ie1]=0;
    }
    for (id=0;id<iesize;id++) {
        piv=0.0;
        for (i=ie1;i<ie2;i++) {
            if (indxr[i-ie1] == 0) {
                big=0.0;
                for (j=je1;j<je2;j++) {
                    if (fabs(s[i][j]) > big) {
                        jp=j;
                        big=fabs(s[i][j]);
                    }
                }
                if (big*pscl[i-ie1] > piv) {
                    ipiv=i;
                    jpiv=jp;
                    piv=big*pscl[i-ie1];
                }
            }
        }
        if (s[ipiv][jpiv] == 0.0)
            nrerror("Singular matrix in routine pinvs");
        indxr[ipiv-ie1]=jpiv+1;
        pivinv=1.0/s[ipiv][jpiv];
        for (j=je1;j<=jsf;j++) s[ipiv][j] *= pivinv;
        s[ipiv][jpiv]=1.0;
        for (i=ie1;i<ie2;i++) {
            if (indxr[i-ie1] != jpiv+1) {
                if (s[i][jpiv] != 0.0) {
                    dum=s[i][jpiv];
                    for (j=je1;j<=jsf;j++)
                        s[i][j] -= dum*s[ipiv][j];
                    s[i][jpiv]=0.0;
                }
            }
        }
    }
    jcoff=jc1-je2;
    icoff=ie1-je1;
    for (i=ie1;i<ie2;i++) {
        irow=indxr[i-ie1]+icoff;
        for (j=je2;j<=jsf;j++) c[irow-1][j+jcoff][k]=s[i][j];
    }
}

void difeq(const int k, const int k1, const int k2, const int jsf,
           const int is1, const int isf, Vec_I_INT &indexv, Mat_O_DP &s,
           Mat_I_DP &y);

void solvde(const int itmax, const DP conv, const DP slowc,
    Vec_I_DP &scalv, Vec_I_INT &indexv, const int nb, Mat_IO_DP &y)
{
    int ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8,j9;
    int jc1,jcf,jv,k,k1,k2,km,kp,nvars;
    DP err,errj,fac,vmax,vz;

    int ne=y.nrows();
    int m=y.ncols();
    Vec_INT kmax(ne);
    Vec_DP ermax(ne);
    Mat3D_DP c(ne,ne-nb+1,m+1);
    Mat_DP s(ne,2*ne+1);
    k1=0; k2=m;
    nvars=ne*m;
    j1=0,j2=nb,j3=nb,j4=ne,j5=j4+j1;
    j6=j4+j2,j7=j4+j3,j8=j4+j4,j9=j8+j1;
    ic1=0,ic2=ne-nb,ic3=ic2,ic4=ne;
    jc1=0,jcf=ic3;
    for (it=0;it<itmax;it++) {
        k=k1;
        difeq(k,k1,k2,j9,ic3,ic4,indexv,s,y);
        pinvs(ic3,ic4,j5,j9,jc1,k1,c,s);
        for (k=k1+1;k<k2;k++) {
            kp=k;
            difeq(k,k1,k2,j9,ic1,ic4,indexv,s,y);
            red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,s);
            pinvs(ic1,ic4,j3,j9,jc1,k,c,s);
        }
        k=k2;
        difeq(k,k1,k2,j9,ic1,ic2,indexv,s,y);
        red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,c,s);
        pinvs(ic1,ic2,j7,j9,jcf,k2,c,s);
        bksub(ne,nb,jcf,k1,k2,c);
        err=0.0;
        for (j=0;j<ne;j++) {
            jv=indexv[j];
            errj=vmax=0.0;
            km=0;
            for (k=k1;k<k2;k++) {
                vz=fabs(c[jv][0][k]);
                if (vz > vmax) {
                    vmax=vz;
                    km=k+1;
                }
                errj += vz;
            }
            err += errj/scalv[j];
            ermax[j]=c[jv][0][km-1]/scalv[j];
            kmax[j]=km;
        }
        err /= nvars;
        fac=(err > slowc ? slowc/err : 1.0);
        for (j=0;j<ne;j++) {
            jv=indexv[j];
            for (k=k1;k<k2;k++)
            y[j][k] -= fac*c[jv][0][k];
        }
        cout << setw(8) << "Iter.";
        cout << setw(10) << "Error" << setw(10) <<  "FAC" << endl;
        cout << setw(6) << it;
        cout << fixed << setprecision(6) << setw(13) << err;
        cout << setw(12) << fac << endl;
        if (err < conv) return;
    }
    nrerror("Too many iterations in solvde");
}

void fdjac(Vec_IO_DP &x, Vec_I_DP &fvec, Mat_O_DP &df,
           void vecfunc(Vec_I_DP &, Vec_O_DP &));
void fdjac(double t, Vec_IO_DP &x, Vec_I_DP &fvec, Mat_O_DP &df,
           void vecfunc(double, Vec_I_DP &, Vec_O_DP &));

int solvde_nb;
const Vec_DP *solvde_xp;
void (*solvde_bc1)(const Vec_DP&, Vec_DP&);
void (*solvde_bc2)(const Vec_DP&, Vec_DP&);
void (*solvde_diff_eq)(double, const Vec_DP&, Vec_DP&);

void difeq(const int k, const int k1, const int k2, const int jsf,
       const int is1, const int isf, Vec_I_INT &indexv, Mat_O_DP &s,
       Mat_I_DP &y)
{
    int i,j,l,n=y.nrows(),m=y.ncols();
    Vec_DP v(n);
    if(k==k1) {// left boundary
        Vec_DP B(solvde_nb);
        Mat_DP J(solvde_nb,n);
        for(i=0; i<n; i++) v[i] = y[i][k1];
        solvde_bc1(v,B);
        fdjac(v,B,J,solvde_bc1);
        for(i=0; i<solvde_nb; i++) {
            l = n - solvde_nb + i;
            s[l][jsf] = B[i];
            for(j=0; j<n; j++)
                s[l][n+indexv[j]] = J[i][j];
        }
    }
    else if(k>=k2) {// right boundary
        l = n - solvde_nb;
        Vec_DP B(l);
        Mat_DP J(l,n);
        for(i=0; i<n; i++) v[i] = y[i][m-1];
        solvde_bc2(v,B);
        fdjac(v,B,J,solvde_bc2);
        for(i=0; i<l; i++) {
            s[i][jsf] = B[i];
            for(j=0; j<n; j++)
                s[i][n+indexv[j]] = J[i][j];
        }
    }
    else {// build finite difference equation
        Vec_DP E(n);
        Mat_DP J(n,n);
        const Vec_DP &x(*solvde_xp);
        double u(0.5*(x[k-1] + x[k])), h(x[k] - x[k-1]);
        for(i=0; i<n; i++) v[i] = 0.5*(y[i][k-1] + y[i][k]);
        solvde_diff_eq(u,v,E);
        fdjac(u,v,E,J,solvde_diff_eq);
        for(i=0; i<n; i++) {
            s[i][jsf] = y[i][k] - y[i][k-1] - h*E[i];
            for(j=0; j<n; j++) {
                s[i][indexv[j]] = -0.5*h*J[i][j];
                s[i][n+indexv[j]] = s[i][indexv[j]];
                if(i==j) {
                    s[i][indexv[j]] -= 1;
                    s[i][n+indexv[j]] += 1;
                }
            }
        }
    }
}

void solvde(Mat_DP &y, const Vec_DP& x, const Vec_DP& scalv,
        void diff_eq(double, const Vec_DP&, Vec_DP&),
        void bc1(const Vec_DP&, Vec_DP&),
        void bc2(const Vec_DP&, Vec_DP&),
        int nb, double slowc,
        double eps, int maxiter)
// solve boundaray value problem by relaxation method
// input:
//   y = initial guess for solution
//   x = mesh points of independent variable
//   scalv = typical scale of solution y
//   diff_eq = right hand side of differential equation dy/dx = f(x,y)
//   bc1 = left boundary condition B(y)=0 computed as bc1(y,B)
//   bc2 = right boundary condition C(y)=0 computed as bc2(y,C)
//   nb = number of boundary condition at left boundary at x[0]
//   slowc = controls step size (Numerical Recipes, section 17.3)
//   eps = error torelance
//   maxiter = maximum number of iterations
// output:
//   y = solution
{
    int i,n=y.nrows(),m=y.ncols();
    Vec_INT indexv(n);
    for(i=0; i<n; i++) indexv[i] = i;
    solvde_xp = &x;
    solvde_nb = nb;
    solvde_bc1 = bc1;
    solvde_bc2 = bc2;
    solvde_diff_eq = diff_eq;
    solvde(maxiter, eps, slowc, scalv, indexv, nb, y);
}

void solvde(Mat_DP &y, double a, double b, const Vec_DP& scalv,
            void diff_eq(double, const Vec_DP&, Vec_DP&),
            void bc1(const Vec_DP&, Vec_DP&),
            void bc2(const Vec_DP&, Vec_DP&),
            int nb, double slowc,
            double eps, int maxiter)
// solve boundaray value problem in range [a,b]
// mesh points on x-axis are a + i*h (i=0,1,...,m-1)
//   where h=(b-a)/(m-1) and m=y.ncols()
{
    int i,m=y.ncols();
    double dx((b-a)/(m-1));
    Vec_DP x(m);
    for(i=0; i<m; i++) x[i] = i*dx;
    solvde(y,x,scalv,diff_eq,bc1,bc2,nb,slowc,eps,maxiter);
}
