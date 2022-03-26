#include "zams_shoot.h"
#include "constants.h"
#include<cmath>
#include<fstream>

main() {
    double X(0.71), Y(0.27);
    ZAMS_shoot sun(M_solar, X, Y);
    std::ofstream ofs("fig5.txt");
    int i,n(1<<8);
    Vec_DP r(n+1);
    Mat_DP y(4,n+1);
    double x(1.e-6), dx(pow(x, 1./n));
    for(i=0; i<=n; i++) r[i] = sun.radius()*(1-pow(dx,i));
    sun.get(y,r);
    double rho, kappa, del, del_ad, del_rad, delta, g;
    double *m(y[0]), *l(y[1]), *P(y[2]), *T(y[3]);
    double PSG = 3./64./PI/Stef/Grav;
    for(i=1; i<=n; i++) {
        eos(rho, del_ad, delta, kappa, P[i], T[i], X, Y);
        del_rad = l[i]/m[i]*PSG*kappa*P[i]/pow(T[i],4);
        g = Grav*m[i]/(r[i]*r[i]);
        if(del_rad <= del_ad) del = del_rad;
        else mixlen(del, del_ad, del_rad, delta, P[i], T[i], rho, kappa, g);
        ofs << P[i] << '\t';
        ofs << T[i] << '\t';
        ofs << del_ad << '\t';
        ofs << del_rad << '\t';
        ofs << del << '\n';
    }
}