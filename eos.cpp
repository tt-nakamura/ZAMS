// reference:
//   R. Kippenhahn, A. Weigert and A. Weiss
//    "Stellar Structure and Evolution" 2nd edition
//     chapter 13

#include<cmath>
#include "constants.h"
#include "nr.h"

void eos_rad(double& rho, double& del_ad, double& delta, double& P_e,
     double P, double T, double X, double Y)
// equation of state, taking radiation pressure into account
// assume complete ionization
// input:
//   P = pressure / dyne/cm^2
//   T = temperature / K
//   X,Y = hydrogen and helium mass fraction
// output:
//   rho = density / g/cm^3
//   del_ad = (dlnT/dlnP)_adiabatic
//   delta = -(dln(rho)/dlnT)_P
//   P_e = electron pressure / dyne/cm^2
// reference: R. Kippenhahn and A. Weigert
//   "Stellar Structure and Evolution" chapter 13
{
    double mu, mu_e, P_rad, b, b1, b2;
    mu = 2./(1 + 3*X + 0.5*Y);// mean molecular weight
    mu_e = 2./(1+X);
    P_rad = a_rad*pow(T,4)/3.;// radiation pressure
    b1 = P_rad/P;
    b = 1-b1;
    b1 *= 4+b;
    b2 = b*b;
    del_ad = (b2+b1)/(2.5*b2 + 4*b1);
    delta = 4/b-3;
    P -= P_rad;
    P *= mu;
    P_e = P/mu_e;
    rho = amu*P/(kB*T);
}

void eos_ion(double& rho, double& del_ad, double& delta, double& P_e,
         double P, double T, double X, double Y);
void opacity(double& KAPPA, double Y, double Z, double T, double RHO, double PE);

void eos(double& rho, double& del_ad, double& delta, double& kappa,
         double P, double T, double X, double Y)
// equation of state for stellar matter
//   if T>=T0, assume complete ionization
//   if T<=T0, take partial ionization into account
// input:
//   P = pressure / dyne/cm^2
//   T = temperature / K
//   X,Y = hydrogen and helium mass fraction
// output:
//   rho = density / g/cm^3
//   del_ad = (dlnT/dlnP)_adiabatic
//   delta = -(dln(rho)/dlnT)_P
//   P_e = electron pressure / dyne/cm^2
//   kappa = opacity / cm^2/g
{
    static double T1(.9e6), T2(1.1e6);
    static double T0(.5*(T1+T2)), Delta_T(.05*(T2-T1));
    double P_e;
    if(T<=T1)      eos_ion(rho, del_ad, delta, P_e, P,T,X,Y);
    else if(T>=T2) eos_rad(rho, del_ad, delta, P_e, P,T,X,Y);
    else {// interpolate to erase discontinuity
        double rho1, del_ad1, delta1, P_e1, f_step;
        eos_rad(rho, del_ad, delta, P_e, P,T,X,Y);
        eos_ion(rho1, del_ad1, delta1, P_e1, P,T,X,Y);
        f_step = 1/(1 + exp((T-T0)/Delta_T));
        rho += (rho1 - rho)*f_step;
        del_ad += (del_ad1 - del_ad)*f_step;
        delta += (delta1 - delta)*f_step;
        P_e += (P_e1 - P_e)*f_step;
    }
    opacity(kappa, Y, 1-X-Y, T, rho, P_e);
}