#include<fstream>
#include<cmath>

void eos(double& rho, double& kappa, double P, double T, double X, double Y);

main() {
    std::ofstream ofs("fig2.txt");
    int i,j,n(100),m(4);
    double X(0.71), Y(0.27), T, kappa, rho;
    double T1(1.e2), T2(1.e8), dT(pow(T2/T1, 1./n));
    double P[] = { 1.e16, 1.e12, 1.e8, 1.e4 };
    for(i=0; i<=n; i++) {
        T = T1*pow(dT,i);
        ofs << T;
        for(j=0; j<m; j++) {
            eos(rho, kappa, P[j], T, X, Y);
            ofs << ' ' << kappa;
        }
        ofs << std::endl;
    }
}
