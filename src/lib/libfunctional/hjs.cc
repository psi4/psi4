#include <cmath>

double hjs_F(double rho, double s, double omega, double* Fhjs, double *Fhjs_rho, double* Fhjs_s)
{
    // => Parameters <= //

    double A, B, C, D, E;
    double s0, k0;
    
    double Ha[10];
    double Hb[10];

    A =  0.757211;
    B = -0.106364;
    C = -0.118649;
    D =  0.609650;
    E = -0.0477963;

    s0 = 2.0;
    k0 = pow(3.0 * M_PI * M_PI, 1.0/3.0);

    Ha[0] =  0.0000000;
    Ha[1] =  0.0000000;
    Ha[2] =  0.0159941;
    Ha[3] =  0.0852995;
    Ha[4] = -0.160368;
    Ha[5] =  0.152645; 
    Ha[6] = -0.0971263;
    Ha[7] =  0.0422061;
    Ha[8] =  0.0000000;
    Ha[9] =  0.0000000;

    Hb[0] =  1.0000000;
    Hb[1] =  5.3331900;
    Hb[2] = -12.478000;
    Hb[3] =  11.098800;
    Hb[4] = -5.110130;
    Hb[5] =  1.714680; 
    Hb[6] = -0.610380; 
    Hb[7] =  0.307555; 
    Hb[8] = -0.0770547;
    Hb[9] =  0.0334840;

    // => Convenience <= //
    
    double rho13 = pow(rho, 1.0/3.0);
    double rho43 = rho * rho13;
    
    // => Zero-th Order Partials <= //
    
    // => H <= //
    
    double H;
    {
        double HN = 0.0;
        double HD = 0.0;
        
        double temp1 = 1.0;
        for (int i = 0; i <= 8; i++) {
            HN += Ha[i] * temp1;
            HD += Hb[i] * temp1;
            temp1 *= s;
        } 

        H = HN / HD;
    }

    // => Zeta <= //
    
    double zeta = s * s * H;

    // => Eta <= //

    double eta = zeta + A;

    // => Lambda <= //

    double lambda = zeta + D;

    // => F <= //
    
    double F = 1.0 - 1.0 / (27.0 * C) * s * s / (1.0 + s * s / (s0 * s0)) - 1.0 / (2.0 * C) * zeta;

    // => G <= //
    
    double G = 1.0/E * (-2.0/5.0 * C * F * lambda - 4.0/15.0 * B * lambda * lambda - 6.0/5.0 * A * lambda * lambda * lambda - 4.0/5.0 * sqrt(M_PI) * pow(lambda,7.0/2.0) - 12.0/5.0 * pow(lambda,7.0/2.0) * (sqrt(zeta) - sqrt(eta)));

    // => Nu <= //

    double nu = omega / (k0 * rho13);

    // => Chi <= //
    
    double chi = nu / sqrt(lambda + nu * nu);
    
    // => Q (Cancellation Hazard for large nu) <= //

    double Q = 2.0 * nu * (sqrt(zeta + nu * nu) - sqrt(eta + nu * nu));

    // => Factor <= //
    
    double Factor = A - 4.0/9.0 * B / lambda * (1.0 - chi) - 4.0/9.0 * C * F / (lambda * lambda) * (1.0 - 3.0/2.0 * chi + 1.0/2.0 * chi * chi) - 8.0/9.0 * E * G / (lambda * lambda * lambda) * (1 - 15.0 / 8.0 * chi + 5.0 / 4.0 * chi * chi *chi - 3.0 / 8.0 * chi * chi * chi * chi * chi) + Q + 2.0 * zeta * log((nu + sqrt(zeta + nu * nu)) / (nu + sqrt(lambda + nu * nu))) + 2.0 * eta * log((nu + sqrt(eta + nu * nu)) / (nu + sqrt(lambda + nu * nu)));
    
    *Fhjs = Factor;
    *Fhjs_rho = 0.0;
    *Fhjs_s = 0.0;
}
