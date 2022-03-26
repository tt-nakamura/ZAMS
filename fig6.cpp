#include "zams_shoot.h"
#include "constants.h"
#include<cmath>
#include<fstream>

main() {
    double X(0.71), Y(0.27);
    ZAMS_shoot sun(M_solar, X, Y);
    std::ofstream ofs("fig6.txt");
    int i,n(1<<8);
    Vec_DP r(n+1);
    Mat_DP y(4,n+1);
    double x(1.e-6), dx(pow(x, 1./n));
    for(i=0; i<=n; i++) r[i] = sun.radius()*(1-pow(dx,i));
    sun.get(y,r);
    double rho, P_e, del_ad, delta;
    double *P(y[2]), *T(y[3]);
    extern double xh, xhe1, xhe2;
    for(i=0; i<=n; i++) {
        eos_ion(rho, del_ad, delta, P_e, P[i], T[i], X, Y);
        ofs << P[i] << '\t';
        ofs << T[i] << '\t';
        ofs << xh   << '\t';
        ofs << xhe1 << '\t';
        ofs << xhe2 << '\n';
    }
}