// reference:
//   R. Kippenhahn, A. Weigert and A. Weiss
//    "Stellar Structure and Evolution" 2nd edition
//     chapter 12

#include<cmath>
#include "zams.h"
#include "constants.h"

const char* ZAMS::init_file = "zams_M1X71Y27.dat";

ZAMS::ZAMS(double M,// set mass / gram
           double X,// set hydrogen mass fraction
           double Y)// set helium mass fraction
{
    load(init_file);
    build(M,X,Y);
}

ZAMS::ZAMS(const char *f) { load(f ? f : init_file); }

static ZAMS *zm;
static void zm_center(const Vec_DP& y, Vec_DP& f) { zm->center(y,f); }
static void zm_surface(const Vec_DP& y, Vec_DP& f) { zm->surface(y,f); }
static void zm_diff_eq(double x, const Vec_DP& y, Vec_DP& f) { zm->diff_eq(y,f); }

void ZAMS::build(double M, double X1, double Y1, double slowc)
// slowc = control parameter of relaxation method (0<slowc<=1)
//   (slowc must be small when initial guess is bad)
{
    static double dM(1.01);// gradual change in M
    Vec_DP scale(1.,6);
    X=X1; Y=Y1;
    zm = this;
    do {// adjust mass gradually
        if(mass<M) { if((mass*=dM)>M) mass=M; }
        else       { if((mass/=dM)<M) mass=M; }
        scale[r_idx] = radius();
        scale[m_idx] = mass;
        scale[l_idx] = luminosity();
        scale[P_idx] = central_pressure();
        scale[T_idx] = central_temperature();
        solvde(data, 0, 1, scale, zm_diff_eq, zm_center, zm_surface, 3, slowc);
    } while(mass!=M);
}

void ZAMS::diff_eq(const Vec_DP& y, Vec_DP& f)
// y = dependent variable (r,m,l,P,T,A)
// ouput: f = dy/dx
// reference: P. P. Eggleton
//   Monthly Notices of Royal Astronomical Society 151 (1971) 351
{
    const double& r = y[r_idx];// distance from center / cm
    const double& m = y[m_idx];// mass within r / g
    const double& l = y[l_idx];// luminosity at r / erg/s
    const double& P = y[P_idx];// pressure at r / dyne/cm^2
    const double& T = y[T_idx];// temperature at r / K
    const double& A = y[A_idx];
    double& dr_dx = f[r_idx];
    double& dm_dx = f[m_idx];
    double& dl_dx = f[l_idx];
    double& dP_dx = f[P_idx];
    double& dT_dx = f[T_idx];
    double& dA_dx = f[A_idx];
    double r2(r*r), rho, kappa, epsilon;
    double delta, del, del_rad, del_ad, g;
    double c1(1/radius()), c2(1/luminosity());
    double c3(1/log(central_pressure()/surface_pressure()));
    static double PSG = 3./64./PI/Stef/Grav;
    eos(rho, del_ad, delta, kappa, P,T,X,Y);
    nuclear(epsilon,rho,T,X,Y);
    g = Grav*m/r2;
    del_rad = l/m*PSG*kappa*P/pow(T,4);
    if(del_rad <= del_ad) del = del_rad;
    else mixlen(del, del_ad, del_rad, delta, P, T, rho, kappa, g);
    double dm_dr = 4*PI*rho*r2;
    double dP_dr = -g*rho;
    double dl_dr = dm_dr*epsilon;
    double dlnP_dr = dP_dr/P;
    double dT_dr = dlnP_dr*T*del;
    dA_dx = 0;
    dr_dx = A/(c1 + c2*dl_dr - c3*dlnP_dr);// adaptive mesh
    dm_dx = dm_dr*dr_dx;
    dl_dx = dl_dr*dr_dx;
    dP_dx = dP_dr*dr_dx;
    dT_dx = dT_dr*dr_dx;
}

void ZAMS::center(const Vec_DP& y, Vec_DP& f)
// boundary condition at r=0
{
    f[0] = y[r_idx];// r=0 at x=0
    f[1] = y[m_idx];// m=0 at x=0
    f[2] = y[l_idx];// l=0 at x=0
}

void ZAMS::surface(const Vec_DP& y, Vec_DP& f)
// boundary condition at r=R
{
    const double& R = y[r_idx];// stellar radius / cm
    const double& P = y[P_idx];// pressure at R / dyne/cm^2
    const double& T = y[T_idx];// temperature at R / K
    double R2(R*R), kappa, dum, dum1, dum2, L1, P1;
    static double PS4(4*PI*Stef);
    eos(dum,dum1,dum2,kappa, P,T,X,Y);
    L1 = PS4*R2*pow(T,4);
    P1 = 2/3.*Grav*mass/R2/kappa;
    f[0] = y[m_idx] - mass;// m=M at x=1
    f[1] = y[l_idx] - L1;// l=L at x=1
    f[2] = y[P_idx] - P1;// atmospheric equilibrium
}

#include<fstream>

void ZAMS::save(const char *fname) {
    int i,j;
    double A(data[A_idx][0]);
    std::ofstream s(fname, std::ofstream::binary);
    s.write((const char*)&X, sizeof(X));
    s.write((const char*)&Y, sizeof(Y));
    s.write((const char*)&N, sizeof(N));
    s.write((const char*)&A, sizeof(A));
    for(j=0; j<=N; j++)
        for(i=0; i<5; i++)
            s.write((const char*)&data[i][j], sizeof(double));
}

void ZAMS::load(const char *fname) {
    int i,j;
    double A;
    std::ifstream s(fname, std::ifstream::binary);
    s.read((char *)&X, sizeof(X));
    s.read((char *)&Y, sizeof(Y));
    s.read((char *)&N, sizeof(N));
    s.read((char *)&A, sizeof(A));
    data = Mat_DP(6,N+1);
    for(j=0; j<=N; j++) {
        for(i=0; i<5; i++)
            s.read((char *)&data[i][j], sizeof(double));
        data[A_idx][j] = A;
    }
    mass = data[m_idx][N];
}
