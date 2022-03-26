// reference:
//   R. Kippenhahn, A. Weigert and A. Weiss
//    "Stellar Structure and Evolution" 2nd edition
//     chapter 12

#include<cmath>
#include "zams_shoot.h"
#include "constants.h"

double ZAMS_shoot::xf = 0.2;

ZAMS_shoot::ZAMS_shoot(
    double M, // set mass / gram
    double X1,// set hydrogen mass fraction
    double Y1 // set helium mass fraction
) : param(4) {
    mass = M_solar;
    // initial guess of parameters
    param[R_shoot] = 0.758317*R_solar;
    param[L_shoot] = 0.499495*L_solar;
    param[P_shoot] = 0.524567*P_solar;
    param[T_shoot] = 0.836627*T_solar;
    build(M,X1,Y1);
}

static ZAMS_shoot *zs;
static void zs_shoot(const Vec_DP& x, Vec_DP& y) { zs->shoot(x,y); }
static void zs_diff_eq(double x, const Vec_DP& y, Vec_DP& f) { zs->diff_eq(x,y,f); }
static double scale[4];

void ZAMS_shoot::build(double M, double X1, double Y1) {
    static double dM(1.01);// gradual change in M
    X=X1; Y=Y1;
    zs = this;
    do {// adjust mass gradually
        if(mass<M) { if((mass*=dM)>M) mass=M; }
        else       { if((mass/=dM)<M) mass=M; }
        scale[m_difeq] = mass;
        scale[l_difeq] = luminosity();
        scale[P_difeq] = central_pressure();
        scale[T_difeq] = central_temperature();
        if(newt(param, zs_shoot))
            nrerror("shoot failed: bad initial guess");
    } while(mass!=M);
}

void ZAMS_shoot::diff_eq(double r, const Vec_DP& y, Vec_DP& f)
// r = distance from center / cm
// y = dependent variables [m,P,T]
// output: f = dy/dr
{
    const double& m = y[m_difeq];// mass within r / g
    const double& l = y[l_difeq];// luminosity at r / erg/s
    const double& P = y[P_difeq];// pressure at r / dyne/cm^2
    const double& T = y[T_difeq];// temperature at r / K
    double& dm_dr = f[m_difeq];
    double& dl_dr = f[l_difeq];
    double& dP_dr = f[P_difeq];
    double& dT_dr = f[T_difeq];
    double r2(r*r), rho, kappa, epsilon;
    double delta, del, del_rad, del_ad, g;
    static double PSG = 3./64./PI/Stef/Grav;
    eos(rho, del_ad, delta, kappa, P,T,X,Y);
    nuclear(epsilon,rho,T,X,Y);
    g = (r==0 ? 0 : Grav*m/r2);
    del_rad = (m==0 ? epsilon : l/m)*PSG*kappa*P/pow(T,4);
    if(del_rad <= del_ad) del = del_rad;
    else mixlen(del, del_ad, del_rad, delta, P, T, rho, kappa, g);
    dm_dr = 4*PI*rho*r2;
    dP_dr = -g*rho;
    dl_dr = dm_dr*epsilon;
    dT_dr = dP_dr*T/P*del;
}

void ZAMS_shoot::shoot(const Vec_DP& x, Vec_DP& y) {
    int i;
    Vec_DP w(4);
    const double& R = x[R_shoot];
    double rf = xf*R;
    center(y,x);// boundary condition at center
    odeint(y, zs_diff_eq, 0, rf);// integrate outward
    surface(w,x);// boundary condition at surface
    odeint(w, zs_diff_eq, R, rf);// integrate inward
    for(i=0; i<4; i++) y[i] -= w[i];
    for(i=0; i<4; i++) y[i] /= scale[i];
    for(i=0; i<4; i++) std::cout << y[i] << ' ';
    std::cout << std::endl;
}

void ZAMS_shoot::center(Vec_DP& y, const Vec_DP& x) {
    y[m_difeq] = 0;            // mass
    y[P_difeq] = x[P_shoot];   // pressure
    y[T_difeq] = x[T_shoot];   // temperature
    y[l_difeq] = 0;            // luminosity
}

static double zs_P, zs_T;
static double zs_surface_eq(double x) { return zs->surface_eq(x); }

void ZAMS_shoot::surface(Vec_DP& y, const Vec_DP& x) {
    const double& R = x[R_shoot];
    const double& L = x[L_shoot];
    double R2(R*R), x1, x2, dx(2), err(1.e-12);
    zs_T = pow(L/(4*PI*R2*Stef), 0.25);
    zs_P = 2./3.*Grav*mass/R2;
    for(x1=2, x2=x1*dx;
        surface_eq(x1)*surface_eq(x2) > 0;
        x1=x2, x2*=dx);
    y[P_difeq] = zbrent(zs_surface_eq, x1, x2, err);
    y[T_difeq] = zs_T;
    y[m_difeq] = mass;
    y[l_difeq] = L;
}

double ZAMS_shoot::surface_eq(double P) {
    double dum,dum1,dum2,kappa;
    eos(dum,dum1,dum2,kappa, P, zs_T, X,Y);
    return P - zs_P/kappa;
}

void ZAMS_shoot::get(Mat_DP& data, const Vec_DP& r_vec)
// r_vec = radial coordinates at which data are evaluated
// data = (m,l,P,T) in cgs (shape(4, len(r_vec)))
{
    int i, j, n(r_vec.size());
    double r(0), rf(xf*radius());
    Vec_DP y(4);
    center(y,param);
    for(j=0; r_vec[j] <= rf; j++) {
        odeint(y, zs_diff_eq, r, r_vec[j]);
        for(i=0; i<4; i++) data[i][j] = y[i];
        r = r_vec[j];
    }
    surface(y,param);
    r = radius();
    for(j=n-1; r_vec[j] > rf; j--) {
        odeint(y, zs_diff_eq, r, r_vec[j]);
        for(i=0; i<4; i++) data[i][j] = y[i];
        r = r_vec[j];
    }
}
