// finite difference jacobian
// W. H. Press, et al, "Numerical Recipes"

#include<cmath>
#include "nr.h"

void fdjac(Vec_IO_DP &x, Vec_I_DP &fvec, Mat_O_DP &df,
           void vecfunc(Vec_I_DP &, Vec_O_DP &))
{
    const DP EPS=1.0e-8;
    int i,j;
    DP h,temp;
    
    int n=x.size();
    int m=fvec.size();
    Vec_DP f(m);
    for (j=0;j<n;j++) {
        temp=x[j];
        h=EPS*fabs(temp);
        if (h == 0.0) h=EPS;
        x[j]=temp+h;
        h=x[j]-temp;
        vecfunc(x,f);
        x[j]=temp;
        for (i=0;i<m;i++)
            df[i][j]=(f[i]-fvec[i])/h;
    }
}

void fdjac(double t, Vec_IO_DP &x, Vec_I_DP &fvec, Mat_O_DP &df,
           void vecfunc(double, Vec_I_DP &, Vec_O_DP &))
{
    const DP EPS=1.0e-8;
    int i,j;
    DP h,temp;
    
    int n=x.size();
    int m=fvec.size();
    Vec_DP f(m);
    for (j=0;j<n;j++) {
        temp=x[j];
        h=EPS*fabs(temp);
        if (h == 0.0) h=EPS;
        x[j]=temp+h;
        h=x[j]-temp;
        vecfunc(t,x,f);
        x[j]=temp;
        for (i=0;i<m;i++)
            df[i][j]=(f[i]-fvec[i])/h;
    }
}