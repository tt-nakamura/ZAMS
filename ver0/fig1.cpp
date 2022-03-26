#include<fstream>
#include<cmath>

void pp_chain(double& epsilon, double rho, double T, double X, double Y);
void CNO_cycle(double& epsilon, double rho, double T, double X, double Y);

main() {
    std::ofstream ofs("fig1.txt");
    int i,n(100);
    double X(0.71), Y(0.27), rho(80), T, e_pp, e_cn;
    double T1(1.e6), T2(1.e8), dT(pow(T2/T1, 1./n));
    for(i=0; i<=n; i++) {
        T = T1*pow(dT,i);
        pp_chain(e_pp, rho, T, X, Y);
        CNO_cycle(e_cn, rho, T, X, Y);
        ofs << T << ' ' << e_pp << ' ' << e_cn << std::endl;
    }
}
