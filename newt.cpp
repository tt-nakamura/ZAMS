// root finding in multi-dimensional space
// W. H. Press, et al, "Numerical Recipes" section 9.7

#include <cmath>
#include <limits>
#include "nr.h"
using namespace std;

void fdjac(Vec_IO_DP &x, Vec_I_DP &fvec, Mat_O_DP &df,
           void vecfunc(Vec_I_DP &, Vec_O_DP &));

void lnsrch(Vec_I_DP &xold, const DP fold, Vec_I_DP &g, Vec_IO_DP &p,
        Vec_O_DP &x, DP &f, const DP stpmax, bool &check, DP func(Vec_I_DP &))
{
    const DP ALF=1.0e-4, TOLX=numeric_limits<DP>::epsilon();
    int i;
    DP a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
    DP rhs1,rhs2,slope,sum,temp,test,tmplam;
    
    int n=xold.size();
    check=false;
    sum=0.0;
    for (i=0;i<n;i++) sum += p[i]*p[i];
    sum=sqrt(sum);
    if (sum > stpmax)
        for (i=0;i<n;i++) p[i] *= stpmax/sum;
    slope=0.0;
    for (i=0;i<n;i++)
        slope += g[i]*p[i];
    if (slope >= 0.0) nrerror("Roundoff problem in lnsrch.");
    test=0.0;
    for (i=0;i<n;i++) {
        temp=fabs(p[i])/MAX(fabs(xold[i]),1.0);
        if (temp > test) test=temp;
    }
    alamin=TOLX/test;
    alam=1.0;
    for (;;) {
        for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
        f=func(x);
        if (alam < alamin) {
            for (i=0;i<n;i++) x[i]=xold[i];
            check=true;
            return;
        } else if (f <= fold+ALF*alam*slope) return;
        else {
            if (alam == 1.0)
                tmplam = -slope/(2.0*(f-fold-slope));
            else {
                rhs1=f-fold-alam*slope;
                rhs2=f2-fold-alam2*slope;
                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                if (a == 0.0) tmplam = -slope/(2.0*b);
                else {
                    disc=b*b-3.0*a*slope;
                    if (disc < 0.0) tmplam=0.5*alam;
                    else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
                    else tmplam=-slope/(b+sqrt(disc));
                }
                if (tmplam>0.5*alam)
                    tmplam=0.5*alam;
            }
        }
        alam2=alam;
        f2 = f;
        alam=MAX(tmplam,0.1*alam);
    }
}

Vec_DP *fvec_p;
void (*nrfuncv)(Vec_I_DP &v, Vec_O_DP &f);

DP fmin(Vec_I_DP &x)
{
    int i;
    DP sum;
    
    Vec_DP &fvec=*fvec_p;
    nrfuncv(x,fvec);
    int n=x.size();
    for (sum=0.0,i=0;i<n;i++) sum += SQR(fvec[i]);
    return 0.5*sum;
}

void ludcmp(Mat_IO_DP &a, Vec_O_INT &indx, DP &d);
void lubksb(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b);

void newt(Vec_IO_DP &x, bool &check, void vecfunc(Vec_I_DP &, Vec_O_DP &))
{
    const int MAXITS=500;
    const DP TOLF=1.0e-8,TOLMIN=1.0e-12,STPMX=1.0e-1;
    const DP TOLX=numeric_limits<DP>::epsilon();
    int i,j,its;
    DP d,den,f,fold,stpmax,sum,temp,test;
    
    int n=x.size();
    Vec_INT indx(n);
    Vec_DP g(n),p(n),xold(n);
    Mat_DP fjac(n,n);
    fvec_p=new Vec_DP(n);
    nrfuncv=vecfunc;
    Vec_DP &fvec=*fvec_p;
    f=fmin(x);
    test=0.0;
    for (i=0;i<n;i++)
        if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < 0.01*TOLF) {
        check=false;
        delete fvec_p;
        return;
    }
    sum=0.0;
    for (i=0;i<n;i++) sum += SQR(x[i]);
    stpmax=STPMX*MAX(sqrt(sum),DP(n));
    for (its=0;its<MAXITS;its++) {
        fdjac(x,fvec,fjac,vecfunc);
        for (i=0;i<n;i++) {
            sum=0.0;
            for (j=0;j<n;j++) sum += fjac[j][i]*fvec[j];
            g[i]=sum;
        }
        for (i=0;i<n;i++) xold[i]=x[i];
        fold=f;
        for (i=0;i<n;i++) p[i] = -fvec[i];
        ludcmp(fjac,indx,d);
        lubksb(fjac,indx,p);
        lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin);
        test=0.0;
        for (i=0;i<n;i++)
            if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
        if (test < TOLF) {
            check=false;
            delete fvec_p;
            return;
        }
        if (check) {
            test=0.0;
            den=MAX(f,0.5*n);
            for (i=0;i<n;i++) {
                temp=fabs(g[i])*MAX(fabs(x[i]),1.0)/den;
                if (temp > test) test=temp;
            }
            check=(test < TOLMIN);
            delete fvec_p;
            return;
        }
        test=0.0;
        for (i=0;i<n;i++) {
            temp=(fabs(x[i]-xold[i]))/MAX(fabs(x[i]),1.0);
            if (temp > test) test=temp;
        }
        if (test < TOLX) {
            delete fvec_p;
            return;
        }
    }
    nrerror("MAXITS exceeded in newt");
}

bool newt(Vec_IO_DP &x, void vecfunc(Vec_I_DP &, Vec_O_DP &))
// find a root x of f(x)=0; vecfunc(x,f) computes f(x)
{
    bool check;
    newt(x,check,vecfunc);
    return check;
}
