// reference:
//   A. Weiss, et al.
//    "Cox and Giuli's Principles of Stellar Structure"

#include<cmath>

void pp_chain(double& e_pp, double rho, double T, double X, double Y) {
    double T6,T63,T633,w,alpha,gamma,psi;
    T6 = T*1.e-6;
    T63 = pow(T6, 1./3);
    T633 = T63*T63;
    w = 1.22e16*exp(-102.6/T63)/sqrt(T63)*X/(1+X); // Cox&Giuli(17.309)
    alpha = 5.48e17*exp(-100/T63)*pow(0.25*Y/X, 2);// Cox&Giuli(17.307)
    gamma = alpha*(sqrt(1 + 2/alpha) - 1);         // Cox&Giuli(17.298)
    psi = 1 + gamma*(0.959 + 0.47*w)/(1+w);        // Cox&Giuli(17.306)
    e_pp = 2.06e6*rho/T633*X*X*exp(-33.81/T63);    // Cox&Giuli(17.310)
    e_pp *= 1 + 0.25*sqrt(rho)*pow(T6, -1.5);      // Cox&Giuli(17.311)
    e_pp *= 1 + 1.2e-3*T63 + 7.8e-3*T633 + 6e-4*T6;// Cox&Giuli(17.312)
    e_pp *= psi;
}

void CNO_cycle(double& e_cn, double rho, double T, double X, double Y) {
    double Z(1-X-Y),T6,T63,T633;
    T6 = T*1.e-6;
    T63 = pow(T6, 1./3);
    T633 = T63*T63;
    e_cn = 7.94e27*rho/T633*X*Z*exp(-152.313/T63); // Cox&Giuli(17.280)
    e_cn *= 1 + 1.75*sqrt(rho)*pow(T6, -1.5);      // Cox&Giuli(17.282)
    e_cn *= 1 + 2.7e-3*T63 - 3.7e-3*T633 - 7e-4*T6;// Cox&Giuli(17.283)
}

void nuclear(double& epsilon, double rho, double T, double X, double Y)
// rho = density / g/cm^3
// T = temperature / K
// X,Y = mass fraction of H,He
// epsilon = nuclear energy generation rate / erg/s/g
{
    double e_pp, e_cn;
    pp_chain(e_pp, rho, T, X, Y);
    CNO_cycle(e_cn, rho, T, X, Y);
    epsilon = e_pp + e_cn;
}
