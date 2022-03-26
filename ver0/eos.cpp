#include<cmath>
#include "constants.h"

void opacity(double& KAPPA, double Y, double Z, double T, double RHO, double PE);

void eos(double& rho, double& kappa, double P, double T, double X, double Y)
// assume complete ionization
// P = pressure / dyne/cm^2
// T = temperature / K
// X,Y = mass fraction of H,He
// rho = density / g/cm^3
// kappa = opacity / cm^2/g
{
    double mu, mu_e, P_rad, P_e, b1, b2;
    mu = 2./(1 + 3*X + 0.5*Y);// mean molecular weight
    mu_e = 2./(1+X);
    P *= mu;
    P_e = P/mu_e;// electron pressure
    rho = P*amu/(kB*T);
    opacity(kappa, Y, 1-X-Y, T, rho, P_e);
}
