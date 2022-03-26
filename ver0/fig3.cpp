#include "zams.h"
#include "constants.h"
#include<fstream>

main() {
    ZAMS sun(M_solar, 0.71, 0.27);
    std::ofstream ofs("fig3.txt");
    int i,n(1<<8);
    double dr(sun.radius()/n);
    Mat_DP y(4,n+1);
    Vec_DP r(n+1);
    for(i=0; i<=n; i++) r[i] = i*dr;
    sun.get(y,r);
    for(i=0; i<=n; i++) {
        ofs << r[i]/R_solar << ' ';
        ofs << y[0][i]/M_solar << '\t';
        ofs << y[1][i]/L_solar << '\t';
        ofs << y[2][i]/P_solar << '\t';
        ofs << y[3][i]/T_solar << '\n';
    }
}
