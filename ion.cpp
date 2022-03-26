// reference:
//   R. Kippenhahn, A. Weigert and A. Weiss
//    "Stellar Structure and Evolution" 2nd edition
//     chapter 14

#include<cmath>
#include "nr.h"
#include "constants.h"

double chi_H   = 13.598*eV;// ionization energy of H / erg
double chi_He0 = 24.587*eV;// ionization energy of neutral He / erg
double chi_He1 = 54.418*eV;// ionization energy of singly ionized He / erg

double KH, KH1, KH2, XH, XHe, xh, xhe1, xhe2;

double saha_eq(double x) {// saha equation
    double y = 1+1/x;
    xh = KH*y;
    xh = xh/(1+xh);
    xhe1 = KH1*y;
    xhe2 = KH2*y*xhe1;
    y = 1/(1 + xhe1 + xhe2);
    xhe1 *= y;
    xhe2 *= y;
    x -= XH*xh + XHe*(xhe1 + 2*xhe2);
    return x;
}

double zbrent(double func(double), double x1, double x2, double tol);
void ludcmp(Mat_DP &a, Vec_INT &indx, double &d);
void lubksb(const Mat_DP &a, const Vec_INT &indx, Vec_DP &b);

void eos_ion(double& rho, double& del_ad, double& delta, double& P_e,
             double P, double T, double X, double Y)
// equation of state,
//   taking partial ionization of H and He into account
// input:
//   P = gas pressure / dyne/cm^2
//   T = temperature / K
//   X,Y = mass fraction of H,He
// output:
//   rho = density / g/cm^3
//   del_ad = (dlnT/dlnP)_adiabatic
//   delta = -(dln(rho)/dlnT)_P
//   P_e = pressure of electron gas / dyne/cm^2
// reference: R. Kippenhahn and A. Weigert
//   "Stellar Structure and Evolution" chapter 14
{
    static double K(pow(sqrt(2*PI*me*kB)/hP, 3)*kB);
    static double err(1e-12), xmin(1e-100), dx(0.5);
    double kT = kB*T;
    double KT5 = K*pow(T, 2.5);
    double HT(chi_H/kT), He0T(chi_He0/kT), He1T(chi_He1/kT);
    double Z(1-X-Y), mu0(1/(X + 0.25*Y + Z/12.));
    // mu0 = number of nucleons / number of nucleus
    double x((X+0.5*Y)*mu0),x1,dum,y,z,E,E1;
    Mat_DP A(3,3);
    Vec_DP B(3);
    Vec_INT indx(3);
    KH  = KT5*exp(-HT)/P;
    KH1 = KT5*exp(-He0T)/P*4;
    KH2 = KT5*exp(-He1T)/P;
    XH  = X*mu0;
    XHe = Y*mu0*0.25;
    for(E=saha_eq(x); E*(E1=saha_eq(x1=x*dx))>0; E=E1,x=x1)
        if(x1<xmin) nrerror("Bad initial range in zbrac");
    E = zbrent(saha_eq, x1, x, err);
    // E = number of free electron / number of nucleus
    dum = 1/(E*(E+1));
    y = 1/xhe1/XHe;
    z = 1/(1-xhe1-xhe2)/XHe;
    A[0][0] = 1/(xh*(1-xh)*XH) + dum;
    A[0][1] = dum;
    A[0][2] = dum*2;
    A[1][0] = dum;
    A[1][1] = y + z + dum;
    A[1][2] = z + dum*2;
    A[2][0] = dum;
    A[2][1] = dum - y;
    A[2][2] = 1/xhe2/XHe + dum*2;
    B[0] = 2.5 + HT;
    B[1] = 2.5 + He0T;
    B[2] = 2.5 + He1T;
    ludcmp(A,indx,dum);
    lubksb(A,indx,B);
    dum = 1+E + B[0] + B[1] + 2*B[2];
    delta = dum/(1+E);
    dum = 2.5 + (HT*B[0] + He0T*B[1] + (He0T + He1T)*B[2])/dum;
    del_ad = 1/dum;
    P /= 1+E;
    P_e = P*E;
    rho = mu0*P*amu/kT;
}
