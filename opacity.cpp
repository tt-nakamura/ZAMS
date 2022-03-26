#include<cmath>

void opacity(double& KAPPA, double Y, double Z, double T,
	     double RHO, double PE) {
    /***********************************************
     ! This routine computes the opacity (OP) given the
     ! density (RHO), temperature (T), Hydrogen mass fraction
     ! (X), Helium mass fraction (Y), and electron pressure (Pe
     ! from the equation of state routine called EOS).
     ! Taken from Stellingwerf 1975, ApJ, 195,
     ! 441., with corrections given in Stellingwerf 1975,
     ! ApJ, 199, 705.
     ! This routine is to be used for 0.6<X<0.8, 0.2<Y<0.4,
     ! and 0.001<Z<0.02.
     ! ***********************************************
     ! Temperature and volume in these fits are given in 10^4 K
     ! and 1/RHO.
     ! ***********************************************/
    //! Set up variables.*/
    double T4 = T/1.0e4;
    double V = 1./RHO;
    double V1 = pow(V, 0.35);
    double V2 = sqrt(V1);
    double Y1 = 6.294e-5-(6.0e-5*Y);
    double Y2 = 3.53e6*Y-3.0447e5;
    double Z1 = 21.0*Z+0.979;
    double Z2 = 105.0*Z+0.895;
    /***********************************************
     !  Compute kappa  from equation D3 of Stellingwerf
     !  1975 which involves some continued fractions
     !  (temp1, temp2).
     ! ***********************************************/
    double TEMP1,TEMP2;
    TEMP1 = Y1*V1*sqrt(T4)*pow(T4,3)+1.0/(760.0*pow(T4,5) + 316.0/V2);
    TEMP1 = 1.0/(10.0*pow(T4,6)+1.0/TEMP1);
    TEMP2 = 1780.0*sqrt(T4)*T4*T4/Z1+1.0/(Z1*Y2/pow(T4,10)
            + 2.13e-3*Z2*V2/sqrt(T4)/pow(T4,4));
    TEMP2 = 47.3/pow(T4,8) + 1.0/TEMP2;
    TEMP2 = 1.0/(4.0e3+1.0/TEMP2);
    KAPPA = PE*(4.819e-13*V/T4+TEMP1+TEMP2);
}
