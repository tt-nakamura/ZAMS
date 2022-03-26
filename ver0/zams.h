#ifndef __zams_h__
#define __zams_h__

#include "nr.h"

struct ZAMS {// Zero Age Main Sequence (homogeneous composition)
    static double xf; // fitting point / stellar radius
    double mass;// gram
    double X;   // hydrogen mass fraction
    double Y;   // helium mass fraction
    Vec_DP param;// R, L, P_c, T_c / cgs
    enum { R_shoot, L_shoot, P_shoot, T_shoot };
    enum { m_difeq, l_difeq, P_difeq, T_difeq };
    ZAMS(double M, // set mass / gram
         double X, // set hydrogen mass fraction
         double Y);// set helium mass fraction
    void diff_eq(double, const Vec_DP&, Vec_DP&);
    void shoot(const Vec_DP&, Vec_DP&);
    void center(Vec_DP&, const Vec_DP&);
    void surface(Vec_DP&, const Vec_DP&);
    double surface_eq(double);
    void set(double M, double X, double Y);
    void get(Mat_DP& y, const Vec_DP& r);
    inline double& radius() { return param[R_shoot]; }
    inline double& luminosity() { return param[L_shoot]; }
    inline double& central_pressure() { return param[P_shoot]; }
    inline double& central_temperature() { return param[T_shoot]; }
};

bool newt(Vec_DP&, void(const Vec_DP&, Vec_DP&));
void odeint(Vec_DP&, void(double, const Vec_DP&, Vec_DP&),
            double, double, double=1e-9);
double zbrent(double(double), double, double, double);
void eos(double& rho, double& kappa, double P, double T, double X, double Y);
void nuclear(double& epsilon, double rho, double T, double X, double Y);

#endif // __zams_h__