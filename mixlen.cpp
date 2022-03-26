// reference:
//   R. Kippenhahn, A. Weigert and A. Weiss
//    "Stellar Structure and Evolution" 2nd edition
//     chapter 7

#include<cmath>
#include "constants.h"

double mix_alpha = 1.5;// mixing length / scale height

void mixlen(double& del, double del_ad, double del_rad, double delta,
            double P, double T, double rho, double kappa, double g)
// mixing-length theory
// input:
//   del_ad = (dlnT/dlnP)_adiabatic
//   del_rad = (dlnT/dlnP)_radiative
//   delta = -(dln(rho)/dlnT)_P
//   P = pressure / dyne/cm^2
//   T = temperature / K
//   rho = density / g/cm^3
//   kappa = opacity / cm^2/g
//   g = grav acceleration / cm/s^2
// output:
//   del = dlnT/dlnP for convection
// reference: R. Kippenhahn and A. Weigert
//   "Stellar Structure and Evolution" section 7.2
{
    static double eps(1e-12);
    if(g==0) { del=del_ad; return; }
    double H_P = P/(rho*g);// scale height
    double l_m = mix_alpha*H_P;// mixing length
    double c_P = P*delta/(rho*T*del_ad);
    double U = 12*Stef*pow(T,3)/(c_P*kappa*pow(rho*l_m, 2))
                 *sqrt(8*H_P/(g*delta));
    double W = 8*U*(del_rad - del_ad);
    double U2(16*U*U), x(0), dx;
    do {
        dx = -(((9*x + 8*U)*x + U2)*x - W)
              /((27*x + 16*U)*x + U2);
        x += dx;
    } while(fabs(dx) > eps*fabs(x));
    del = x*(x+2*U) + del_ad;
}