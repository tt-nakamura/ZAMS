#ifndef __zams_h__
#define __zams_h__

#include "nr.h"

struct ZAMS {// Zero Age Main Sequence (homogeneous composition)
    static const char *init_file;
    int N;// number of mesh points - 1
    double mass;// gram
    double X,Y; // mass fraction of H,He
    enum { r_idx, m_idx, l_idx, P_idx, T_idx, A_idx };
    Mat_DP data;// (r, m, L, P, T, A) in cgs
    ZAMS(double M, // set mass / gram
         double X, // set hydrogen mass fraction
         double Y);// set helium mass fraction
    ZAMS(const char* =0);
    void diff_eq(const Vec_DP&, Vec_DP&);
    void center(const Vec_DP&, Vec_DP&);
    void surface(const Vec_DP&, Vec_DP&);
    void build(double M, double X, double Y, double slowc=1);
    void save(const char*);
    void load(const char*);
    inline double& radius() { return data[r_idx][N]; }
    inline double& luminosity() { return data[l_idx][N]; }
    inline double& central_pressure() { return data[P_idx][0]; }
    inline double& surface_pressure() { return data[P_idx][N]; }
    inline double& central_temperature() { return data[T_idx][0]; }
    inline double& surface_temperature() { return data[T_idx][N]; }
};

void solvde(Mat_DP &y, double a, double b, const Vec_DP& scalv,
            void diff_eq(double, const Vec_DP&, Vec_DP&),
            void boundary1(const Vec_DP&, Vec_DP&),
            void boundary2(const Vec_DP&, Vec_DP&),
            int nb, double slowc=1,
            double eps=1.e-8, int maxiter=10000);

void eos(double& rho, double& del_ad, double& delta, double& kappa,
         double P, double T, double X, double Y);
void nuclear(double& epsilon, double rho, double T, double X, double Y);
void mixlen(double& del, double del_ad, double del_rad, double delta,
            double P, double T, double rho, double kappa, double g);

#endif // __zams_h__
