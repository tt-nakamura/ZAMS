#include "zams.h"
#include "constants.h"
#include<cmath>
#include<fstream>

main() {
    std::ofstream ofs("fig7-9.txt");
    int i,n(1<<8);
    double X(0.71), Y(0.27);
    double M1(0.1*M_solar), M2(30*M_solar), dM(pow(M2/M1, 1./n));
    double M[n+1], R[n+1], L[n+1], T[n+1];
    for(i=0; i<=n; i++) M[i] = M1*pow(dM,i);
    ZAMS s,t;
    for(i=0; i<=n; i++) {
        if(M[i] <= M_solar) continue;        
        s.build(M[i],X,Y);
        R[i] = s.radius();
        L[i] = s.luminosity();
        T[i] = s.surface_temperature();
    }
    for(i=n; i>=0; i--) {
        if(M[i] > M_solar) continue;
        t.build(M[i],X,Y);
        R[i] = t.radius();
        L[i] = t.luminosity();
        T[i] = t.surface_temperature();
    }
    for(i=0; i<=n; i++) {
        ofs << M[i] << '\t';// mass
        ofs << R[i] << '\t';// radius
        ofs << L[i] << '\t';// luminosity
        ofs << T[i] << '\n';// surface temperature
    }
}
