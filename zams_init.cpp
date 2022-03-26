#include "zams.h"
#include "zams_shoot.h"
#include "constants.h"
#include<fstream>

main() {
    int i,j,N(1<<10);
    double X(0.71), Y(0.27), A(1), dr;
    const char *fname = "zams_M1X71Y27.dat";
    Mat_DP y(4,N+1);// (m,l,P,T) in cgs
    Vec_DP r(N+1);
    ZAMS_shoot sun(M_solar, X, Y);
    dr = sun.radius()/N;
    for(i=0; i<=N; i++) r[i] = i*dr;
    sun.get(y,r);
    // save initial guess
    std::ofstream s(fname, std::ofstream::binary);
    s.write((const char*)&X, sizeof(X));
    s.write((const char*)&Y, sizeof(Y));
    s.write((const char*)&N, sizeof(N));
    s.write((const char*)&A, sizeof(A));
    for(j=0; j<=N; j++) {
        s.write((const char*)&r[j], sizeof(double));
        for(i=0; i<4; i++)
            s.write((const char*)&y[i][j], sizeof(double));
    }
    s.close();
    ZAMS z(fname);// load initial guess
    // slowc = 1e-3 to avoid singular matrix
    z.build(M_solar, X, Y, 1e-3);
    z.save(fname);
}
