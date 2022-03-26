#include<cmath>
#include "zams.h"
#include "constants.h"

double ZAMS::xf = 0.2;

ZAMS::ZAMS(double M, // set mass / gram
           double X1,// set hydrogen mass fraction
           double Y1 // set helium mass fraction
) : param(4) {
    mass = M_solar;
    param[R_shoot] = R_solar;
    param[L_shoot] = L_solar;
    param[P_shoot] = P_solar;
    param[T_shoot] = T_solar;
    set(M,X1,Y1);
}

ZAMS *zm;
void zm_shoot(const Vec_DP& x, Vec_DP& y) { zm->shoot(x,y); }
void zm_diff_eq(double x, const Vec_DP& y, Vec_DP& f) { zm->diff_eq(x,y,f); }
double scale[4];

void ZAMS::set(double M, double X1, double Y1) {
    static double dM(1.01);
    X=X1; Y=Y1;
    zm = this;
    do {
        if(mass<M) { if((mass*=dM)>M) mass=M; }
        else       { if((mass/=dM)<M) mass=M; }
        scale[m_difeq] = mass;
        scale[l_difeq] = luminosity();
        scale[P_difeq] = central_pressure();
        scale[T_difeq] = central_temperature();
        if(newt(param, zm_shoot))
            nrerror("shoot failed: bad initial guess");
    } while(mass!=M);
}

void ZAMS::diff_eq(double r, const Vec_DP& y, Vec_DP& f)
// r = distance from center / cm
{
    const double& m = y[m_difeq];// mass within r / g
    const double& l = y[l_difeq];// luminosity at r / erg/s
    const double& P = y[P_difeq];// pressure at r / dyne/cm^2
    const double& T = y[T_difeq];// temperature at r / K
    double& dm_dr = f[m_difeq];
    double& dl_dr = f[l_difeq];
    double& dP_dr = f[P_difeq];
    double& dT_dr = f[T_difeq];
    double r2(r*r), rho, kappa, epsilon, del_rad;
    static double PSG = 3./64./PI/Stef/Grav;
    static double del_ad(0.4);// adiabatic dlnT/dlnP
    eos(rho,kappa,P,T,X,Y);
    nuclear(epsilon,rho,T,X,Y);
    del_rad = (m==0 ? epsilon : l/m)*PSG*kappa*P/pow(T,4);
    dm_dr = 4*PI*rho*r2;
    dP_dr = (r==0 ? 0 : -Grav*m*rho/r2);
    dl_dr = dm_dr*epsilon;
    dT_dr = dP_dr/P*T*MIN(del_rad, del_ad);
}

void ZAMS::shoot(const Vec_DP& x, Vec_DP& y) {
    int i;
    Vec_DP w(4);
    const double& R = x[R_shoot];
    double rf = xf*R;
    center(y,x);// boundary condition at center
    odeint(y, zm_diff_eq, 0, rf);// integrate outward
    surface(w,x);// boundary condition at surface
    odeint(w, zm_diff_eq, R, rf);// integrate inward
    for(i=0; i<4; i++) y[i] -= w[i];
    for(i=0; i<4; i++) y[i] /= scale[i];
    for(int i=0; i<4; i++) std::cout << y[i] << ' ';
    std::cout << std::endl;
}

void ZAMS::center(Vec_DP& y, const Vec_DP& x) {
    y[m_difeq] = 0;         // mass
    y[P_difeq] = x[P_shoot];// pressure
    y[T_difeq] = x[T_shoot];// temperature
    y[l_difeq] = 0;         // luminosity
}

double zm_T, zm_P;
double zm_surface_eq(double x) { return zm->surface_eq(x); }

void ZAMS::surface(Vec_DP& y, const Vec_DP& x) {
    const double& R = x[R_shoot];
    const double& L = x[L_shoot];
    double R2(R*R), x1, x2, dx(2), err(1.e-12);
    zm_T = pow(L/(4*PI*R2*Stef), 0.25);
    zm_P = 2./3.*Grav*mass/R2;
    for(x1=2, x2=x1*dx;
        surface_eq(x1)*surface_eq(x2) > 0;
        x1=x2, x2*=dx);
    y[P_difeq] = zbrent(zm_surface_eq, x1, x2, err);
    y[T_difeq] = zm_T;
    y[m_difeq] = mass;
    y[l_difeq] = L;
}

double ZAMS::surface_eq(double P) {
    double dum,kappa;
    eos(dum,kappa, P, zm_T, X,Y);
    return P - zm_P/kappa;
}

void ZAMS::get(Mat_DP& y_mat, const Vec_DP& r_vec)
// r_vec = radial coordinates at which data are evaluated
// data = (m,l,P,T) in cgs (shape(4, len(r_vec)))
{
    int i, j, n(r_vec.size());
    double r(0), rf(xf*radius());
    Vec_DP y(4);
    center(y,param);
    for(j=0; r_vec[j] <= rf; j++) {
        odeint(y, zm_diff_eq, r, r_vec[j]);
        for(i=0; i<4; i++) y_mat[i][j] = y[i];
        r = r_vec[j];
    }
    surface(y,param);
    r = radius();
    for(j=n-1; r_vec[j] > rf; j--) {
        odeint(y, zm_diff_eq, r, r_vec[j]);
        for(i=0; i<4; i++) y_mat[i][j] = y[i];
        r = r_vec[j];
    }
}
