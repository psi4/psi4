/**********************************************************
* dispersion.cc: definitions for -D for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/

#include <libmints/molecule.h>
#include <libciomr/libciomr.h>
#include "dispersion.h"
#include "dispersion_defines.h"
#include <boost/algorithm/string.hpp>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace psi;
using namespace boost;
using namespace std;

namespace psi { namespace functional {


Dispersion::Dispersion()
{
}
Dispersion::~Dispersion()
{
}
std::string Dispersion::printEnergy(boost::shared_ptr<Molecule> m)
{
    double e = computeEnergy(m);
    std::stringstream s;
    s.setf(ios::scientific);
    s.precision(11);

    s << "   " << name_ << " Dispersion Energy: " << e << " [H]" << endl;

    return s.str();
}
std::string Dispersion::printGradient(boost::shared_ptr<Molecule> m)
{
    double* g = computeGradient(m);
    std::stringstream s;
    s.setf(ios::scientific);
    s.precision(11);

    s << "   " << name_ << " Dispersion Gradient ([a.u.]): " << endl << endl;
    s << "    Atom #:       E_x                E_y                 E_z" << endl;
    s << "   -----------------------------------------------------------------" << endl;

    for (int k = 1; k <= m->natom(); k++) {
        s << "  " << setw(5) << k <<  \
            setw(20) << g[(k-1)*3 + 0] << \
            setw(20) << g[(k-1)*3 + 1] << \
            setw(20) << g[(k-1)*3 + 2] << endl;
    }
    return s.str();
}
std::string Dispersion::printHessian(boost::shared_ptr<Molecule> m)
{
    double** h = computeHessian(m);

    std::stringstream s;
    s.setf(ios::scientific);
    s.precision(11);

    //print_mat(h, 3*m->natom(), 3*m->natom(), outfile);

    s << "   " << name_ << " Dispersion Hessian ([a.u.]): " << endl << endl;
    for (int k = 1; k <= m->natom(); k++) {
        for (int j = 1; j <= m->natom(); j++) {

            s << "    Atom Pair A = " << k << " B = "  << j << ":" << endl << endl;
            s << "                   xB                 yB                  zB" << endl;
            s << "   -----------------------------------------------------------------" << endl;

            s << "  " << setw(5) << "xA" <<  \
                setw(20) << h[(k-1)*3 + 0][(j-1)*3 + 0] << \
                setw(20) << h[(k-1)*3 + 0][(j-1)*3 + 1] << \
                setw(20) << h[(k-1)*3 + 0][(j-1)*3 + 2] << endl;
            s << "  " << setw(5) << "yA" <<  \
                setw(20) << h[(k-1)*3 + 1][(j-1)*3 + 0] << \
                setw(20) << h[(k-1)*3 + 1][(j-1)*3 + 1] << \
                setw(20) << h[(k-1)*3 + 1][(j-1)*3 + 2] << endl;
            s << "  " << setw(5) << "zA" <<  \
                setw(20) << h[(k-1)*3 + 2][(j-1)*3 + 0] << \
                setw(20) << h[(k-1)*3 + 2][(j-1)*3 + 1] << \
                setw(20) << h[(k-1)*3 + 2][(j-1)*3 + 2] << endl;
            s << endl;
        }
    }
    return s.str();
}
boost::shared_ptr<Dispersion> Dispersion::createDispersion(const std::string & name, double s6, double s8, double sr6, double sr8)
{
    if (boost::to_upper_copy(name) == "") {
        boost::shared_ptr<Dispersion> disp;
        return disp;
    } else if (boost::to_upper_copy(name) == "-D1") {
        return boost::shared_ptr<Dispersion> (new D1(s6));
    } else if (boost::to_upper_copy(name) == "-D2") {
        return boost::shared_ptr<Dispersion> (new D2(s6));
    } else if (boost::to_upper_copy(name) == "-D3") {
        return boost::shared_ptr<Dispersion> (new D3(s6, s8, sr6, sr8));
    }
}
D1::D1(double s6) : Dispersion()
{
    name_ = "-D1";
    description_ = "Grimme's -D1 Dispersion Correction";
    citation_ = "Grimme, S. (2004), J. Comp. Chem., 25: 1463-1473";

    d_ = 23.0;
    RvdW_ = RvdW_D1_;
    C6_ = C6_D1_;
    s6_ = s6;
}
D1::~D1()
{
}
double D1::computeEnergy(boost::shared_ptr<Molecule> mol)
{
    double energy = 0.0;

    // Build Z, x, y, and z
    int natom = 0;

    int *Z = init_int_array(mol->natom());
    double *x = init_array(mol->natom());
    double *y = init_array(mol->natom());
    double *z = init_array(mol->natom());

    for (int i = 0; i < mol->natom(); i++) {
        if (mol->Z(i) == 0)
        continue;
        Z[natom] = mol->Z(i);
        x[natom] = mol->x(i);
        y[natom] = mol->y(i);
        z[natom] = mol->z(i);
        natom ++;
    }

    // Build C6 (C6_ij = C6_i*C6_j/(C6_i+C6_j) in -D1 )
    double **C6 = block_matrix(natom, natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            C6[i][j] = 2.0 * C6_[Z[i]] * C6_[Z[j]] / (C6_[Z[i]] + C6_[Z[j]]);
            C6[j][i] = C6[i][j];
        }
    }

    // Build RvdW (sum of vdW radii)
    double **RvdW = block_matrix(natom, natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            RvdW[i][j] = RvdW_[Z[i]] + RvdW_[Z[j]];
            RvdW[j][i] = RvdW[i][j];
        }
    }

    // Build r
    double **r = block_matrix(natom,natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            r[i][j] = sqrt((x[i]-x[j])*(x[i]-x[j]) + \
                           (y[i]-y[j])*(y[i]-y[j]) + \
                           (z[i]-z[j])*(z[i]-z[j]));
            r[j][i] = r[i][j];

        }
    }

    // Compute energy
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
              double t3960 = x[i]-x[j];
              double t3961 = y[i]-y[j];
              double t3962 = z[i]-z[j];
              double t3963 = t3960*t3960;
              double t3964 = t3961*t3961;
              double t3965 = t3962*t3962;
              double t3966 = t3963+t3964+t3965;
              energy += (C6[i][j]*1/(t3966*t3966*t3966))/(exp(-d_*(sqrt(t3966)/RvdW[i][j]-1.0))+1.0);
        }
    }

    // Scale energy by s6_
    energy *= s6_;

    // Dispersion energy is always stabilizing
    energy *= -1.0;

    // Free stuff
    free_block(C6);
    free_block(RvdW);
    free_block(r);
    delete []Z;
    delete []x;
    delete []y;
    delete []z;

    return energy;
}
double* D1::computeGradient(boost::shared_ptr<Molecule> mol)
{
    // Build Z, x, y, and z
    int natom = 0;

    int *Z = init_int_array(mol->natom());
    double *x = init_array(mol->natom());
    double *y = init_array(mol->natom());
    double *z = init_array(mol->natom());

    for (int i = 0; i < mol->natom(); i++) {
        if (mol->Z(i) == 0)
        continue;
        Z[natom] = mol->Z(i);
        x[natom] = mol->x(i);
        y[natom] = mol->y(i);
        z[natom] = mol->z(i);
        natom ++;
    }

    // gradient array in [x1 y1 z1 x2 y2 z2 ... ]
    double* grad = init_array(3*natom);

    // Build C6 (C6_ij = C6_i*C6_j/(C6_i+C6_j) in -D1 )
    double **C6 = block_matrix(natom, natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            C6[i][j] = 2.0 * C6_[Z[i]] * C6_[Z[j]] / (C6_[Z[i]] + C6_[Z[j]]);
            C6[j][i] = C6[i][j];
        }
    }

    // Build RvdW (sum of vdW radii)
    double **RvdW = block_matrix(natom, natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            RvdW[i][j] = RvdW_[Z[i]] + RvdW_[Z[j]];
            RvdW[j][i] = RvdW[i][j];
        }
    }

    // Build r
    double **r = block_matrix(natom,natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            r[i][j] = sqrt((x[i]-x[j])*(x[i]-x[j]) + \
                           (y[i]-y[j])*(y[i]-y[j]) + \
                           (z[i]-z[j])*(z[i]-z[j]));
            r[j][i] = r[i][j];

        }
    }

    // Compute gradient by perturbing nucleus i
    int Ax, Ay, Az;
    double f;
    for (int i = 0; i < natom; i++) {
        Ax = 3*i;
        Ay = 3*i + 1;
        Az = 3*i + 2;

        for (int j = 0; j < natom; j++) {
            if (i == j)
                continue;
              double t3968 = x[i]-x[j];
              double t3969 = y[i]-y[j];
              double t3970 = z[i]-z[j];
              double t3971 = t3968*t3968;
              double t3972 = t3969*t3969;
              double t3973 = t3970*t3970;
              double t3974 = t3971+t3972+t3973;
              double t3975 = 1/RvdW[i][j];
              double t3976 = sqrt(t3974);
              double t3977 = t3975*t3976;
              double t3978 = t3977-1.0;
              double t3982 = d_*t3978;
              double t3979 = exp(-t3982);
              double t3980 = x[i]*2.0;
              double t3981 = t3980-x[j]*2.0;
              double t3983 = t3979+1.0;
              grad[Ax] += (C6[i][j]*t3981*1/(t3974*t3974*t3974*t3974)*-3.0)/t3983+C6[i][j]*d_*t3981*1/pow(t3974,7.0/ \
                   2.0)*1/(t3983*t3983)*t3975*t3979*(1.0/2.0);
              double t3985 = x[i]-x[j];
              double t3986 = y[i]-y[j];
              double t3987 = z[i]-z[j];
              double t3988 = t3985*t3985;
              double t3989 = t3986*t3986;
              double t3990 = t3987*t3987;
              double t3991 = t3990+t3988+t3989;
              double t3992 = 1/RvdW[i][j];
              double t3993 = sqrt(t3991);
              double t3994 = t3992*t3993;
              double t3995 = t3994-1.0;
              double t3999 = d_*t3995;
              double t3996 = exp(-t3999);
              double t3997 = y[i]*2.0;
              double t3998 = t3997-y[j]*2.0;
              double t4000 = t3996+1.0;
              grad[Ay] += (C6[i][j]*1/(t3991*t3991*t3991*t3991)*t3998*-3.0)/t4000+C6[i][j]*d_*1/pow(t3991,7.0/2.0) \
                   *t3992*t3996*t3998*1/(t4000*t4000)*(1.0/2.0);
              double t4002 = x[i]-x[j];
              double t4003 = y[i]-y[j];
              double t4004 = z[i]-z[j];
              double t4005 = t4002*t4002;
              double t4006 = t4003*t4003;
              double t4007 = t4004*t4004;
              double t4008 = t4005+t4006+t4007;
              double t4009 = 1/RvdW[i][j];
              double t4010 = sqrt(t4008);
              double t4011 = t4010*t4009;
              double t4012 = t4011-1.0;
              double t4016 = d_*t4012;
              double t4013 = exp(-t4016);
              double t4014 = z[i]*2.0;
              double t4015 = t4014-z[j]*2.0;
              double t4017 = t4013+1.0;
              grad[Az] += (C6[i][j]*t4015*1/(t4008*t4008*t4008*t4008)*-3.0)/t4017+C6[i][j]*d_*t4013*t4015*1/pow(t4008,7.0/ \
                   2.0)*1/(t4017*t4017)*t4009*(1.0/2.0);
        }
    }

    // Scale gradient by -s6_
    // Dispersion energy is always stabilizing
    for (int A = 0; A < 3*natom; A++) {
        grad[A] *= -s6_;
    }

    // Free stuff
    free_block(C6);
    free_block(RvdW);
    free_block(r);
    delete []Z;
    delete []x;
    delete []y;
    delete []z;
    return grad;
}
double** D1::computeHessian(boost::shared_ptr<Molecule> mol)
{
    // Build Z, x, y, and z
    int natom = 0;

    int *Z = init_int_array(mol->natom());
    double *x = init_array(mol->natom());
    double *y = init_array(mol->natom());
    double *z = init_array(mol->natom());

    for (int i = 0; i < mol->natom(); i++) {
        if (mol->Z(i) == 0)
        continue;
        Z[natom] = mol->Z(i);
        x[natom] = mol->x(i);
        y[natom] = mol->y(i);
        z[natom] = mol->z(i);
        natom ++;
    }

    // gradient array in [x1 y1 z1 x2 y2 z2 ... ]
    double** hess = block_matrix(3*natom,3*natom);

    // Build C6 (C6_ij = C6_i*C6_j/(C6_i+C6_j) in -D1 )
    double **C6 = block_matrix(natom, natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            C6[i][j] = 2.0 * C6_[Z[i]] * C6_[Z[j]] / (C6_[Z[i]] + C6_[Z[j]]);
            C6[j][i] = C6[i][j];
        }
    }

    // Build RvdW (sum of vdW radii)
    double **RvdW = block_matrix(natom, natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            RvdW[i][j] = RvdW_[Z[i]] + RvdW_[Z[j]];
            RvdW[j][i] = RvdW[i][j];
        }
    }

    // Build r
    double **r = block_matrix(natom,natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            r[i][j] = sqrt((x[i]-x[j])*(x[i]-x[j]) + \
                           (y[i]-y[j])*(y[i]-y[j]) + \
                           (z[i]-z[j])*(z[i]-z[j]));
            r[j][i] = r[i][j];

        }
    }

    int Ax, Ay, Az, Bx, By, Bz;
    double f;

    // Case 1: A A
    for (int i = 0; i < natom; i++) {
        Ax = 3*i;
        Ay = 3*i + 1;
        Az = 3*i + 2;

        for (int j = 0; j < natom; j++) {
            if (i == j)
                continue;

              double t4019 = x[i]-x[j];
              double t4020 = y[i]-y[j];
              double t4021 = z[i]-z[j];
              double t4022 = t4019*t4019;
              double t4023 = t4020*t4020;
              double t4024 = t4021*t4021;
              double t4025 = t4022+t4023+t4024;
              double t4035 = x[i]*2.0;
              double t4036 = x[j]*2.0;
              double t4026 = t4035-t4036;
              double t4027 = 1/RvdW[i][j];
              double t4028 = sqrt(t4025);
              double t4029 = t4027*t4028;
              double t4030 = t4029-1.0;
              double t4034 = d_*t4030;
              double t4031 = exp(-t4034);
              double t4032 = t4031+1.0;
              double t4033 = 1/t4032;
              double t4037 = t4026*t4026;
              double t4038 = 1/(t4032*t4032);
              double t4039 = 1/(t4025*t4025*t4025*t4025);
              double t4040 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4041 = d_*d_;
              hess[Ax][Ax] += C6[i][j]*t4033*t4039*-6.0+C6[i][j]*t4033*1/(t4025*t4025*t4025*t4025*t4025)*t4037*1.2E1+ \
                   C6[i][j]*d_*t4031*1/pow(t4025,7.0/2.0)*t4027*t4038-C6[i][j]*t4031*t4040*t4041*t4037*t4038*t4039*(1.0/ \
                   4.0)-C6[i][j]*d_*t4031*1/pow(t4025,9.0/2.0)*t4027*t4037*t4038*(1.3E1/4.0)+C6[i][j]*t4040*1/(t4032* \
                   t4032*t4032)*t4041*t4037*t4039*exp(d_*t4030*-2.0)*(1.0/2.0);
              double t4043 = x[i]-x[j];
              double t4044 = y[i]-y[j];
              double t4045 = z[i]-z[j];
              double t4046 = t4043*t4043;
              double t4047 = t4044*t4044;
              double t4048 = t4045*t4045;
              double t4049 = t4046+t4047+t4048;
              double t4050 = 1/RvdW[i][j];
              double t4051 = sqrt(t4049);
              double t4052 = t4050*t4051;
              double t4053 = t4052-1.0;
              double t4059 = d_*t4053;
              double t4054 = exp(-t4059);
              double t4055 = x[i]*2.0;
              double t4063 = x[j]*2.0;
              double t4056 = -t4063+t4055;
              double t4057 = y[i]*2.0;
              double t4064 = y[j]*2.0;
              double t4058 = -t4064+t4057;
              double t4060 = t4054+1.0;
              double t4061 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4062 = d_*d_;
              double t4065 = 1/(t4049*t4049*t4049*t4049);
              double t4066 = 1/(t4060*t4060);
              hess[Ax][Ay] += (C6[i][j]*t4056*1/(t4049*t4049*t4049*t4049*t4049)*t4058*1.2E1)/t4060-C6[i][j]*d_*t4050* \
                   t4054*t4056*t4066*1/pow(t4049,9.0/2.0)*t4058*(1.3E1/4.0)+C6[i][j]*1/(t4060*t4060*t4060)*t4061*t4062* \
                   t4056*t4065*t4058*exp(d_*t4053*-2.0)*(1.0/2.0)-C6[i][j]*t4061*t4062*t4054*t4056*t4065*t4066*t4058* \
                   (1.0/4.0);
              double t4068 = x[i]-x[j];
              double t4069 = y[i]-y[j];
              double t4070 = z[i]-z[j];
              double t4071 = t4068*t4068;
              double t4072 = t4069*t4069;
              double t4073 = t4070*t4070;
              double t4074 = t4071+t4072+t4073;
              double t4075 = 1/RvdW[i][j];
              double t4076 = sqrt(t4074);
              double t4077 = t4075*t4076;
              double t4078 = t4077-1.0;
              double t4084 = d_*t4078;
              double t4079 = exp(-t4084);
              double t4080 = x[i]*2.0;
              double t4088 = x[j]*2.0;
              double t4081 = t4080-t4088;
              double t4082 = z[i]*2.0;
              double t4089 = z[j]*2.0;
              double t4083 = t4082-t4089;
              double t4085 = t4079+1.0;
              double t4086 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4087 = d_*d_;
              double t4090 = 1/(t4074*t4074*t4074*t4074);
              double t4091 = 1/(t4085*t4085);
              hess[Ax][Az] += (C6[i][j]*t4081*1/(t4074*t4074*t4074*t4074*t4074)*t4083*1.2E1)/t4085-C6[i][j]*d_*t4081* \
                   t4091*1/pow(t4074,9.0/2.0)*t4083*t4075*t4079*(1.3E1/4.0)+C6[i][j]*t4081*t4090*t4083*1/(t4085*t4085* \
                   t4085)*t4086*t4087*exp(d_*t4078*-2.0)*(1.0/2.0)-C6[i][j]*t4081*t4090*t4091*t4083*t4086*t4087*t4079* \
                   (1.0/4.0);
              double t4093 = x[i]-x[j];
              double t4094 = y[i]-y[j];
              double t4095 = z[i]-z[j];
              double t4096 = t4093*t4093;
              double t4097 = t4094*t4094;
              double t4098 = t4095*t4095;
              double t4099 = t4096+t4097+t4098;
              double t4109 = y[i]*2.0;
              double t4110 = y[j]*2.0;
              double t4100 = -t4110+t4109;
              double t4101 = 1/RvdW[i][j];
              double t4102 = sqrt(t4099);
              double t4103 = t4101*t4102;
              double t4104 = t4103-1.0;
              double t4108 = d_*t4104;
              double t4105 = exp(-t4108);
              double t4106 = t4105+1.0;
              double t4107 = 1/t4106;
              double t4111 = t4100*t4100;
              double t4112 = 1/(t4106*t4106);
              double t4113 = 1/(t4099*t4099*t4099*t4099);
              double t4114 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4115 = d_*d_;
              hess[Ay][Ay] += C6[i][j]*t4113*t4107*-6.0+C6[i][j]*1/(t4099*t4099*t4099*t4099*t4099)*t4111*t4107*1.2E1+ \
                   C6[i][j]*d_*1/pow(t4099,7.0/2.0)*t4101*t4112*t4105-C6[i][j]*t4111*t4112*t4113*t4105*t4114*t4115*(1.0/ \
                   4.0)-C6[i][j]*d_*1/pow(t4099,9.0/2.0)*t4101*t4111*t4112*t4105*(1.3E1/4.0)+C6[i][j]*t4111*t4113*t4114* \
                   1/(t4106*t4106*t4106)*t4115*exp(d_*t4104*-2.0)*(1.0/2.0);
              double t4117 = x[i]-x[j];
              double t4118 = y[i]-y[j];
              double t4119 = z[i]-z[j];
              double t4120 = t4117*t4117;
              double t4121 = t4118*t4118;
              double t4122 = t4119*t4119;
              double t4123 = t4120+t4121+t4122;
              double t4124 = 1/RvdW[i][j];
              double t4125 = sqrt(t4123);
              double t4126 = t4124*t4125;
              double t4127 = t4126-1.0;
              double t4133 = d_*t4127;
              double t4128 = exp(-t4133);
              double t4129 = y[i]*2.0;
              double t4137 = y[j]*2.0;
              double t4130 = -t4137+t4129;
              double t4131 = z[i]*2.0;
              double t4138 = z[j]*2.0;
              double t4132 = t4131-t4138;
              double t4134 = t4128+1.0;
              double t4135 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4136 = d_*d_;
              double t4139 = 1/(t4123*t4123*t4123*t4123);
              double t4140 = 1/(t4134*t4134);
              hess[Ay][Az] += (C6[i][j]*t4130*1/(t4123*t4123*t4123*t4123*t4123)*t4132*1.2E1)/t4134-C6[i][j]*d_*t4130* \
                   t4140*1/pow(t4123,9.0/2.0)*t4132*t4124*t4128*(1.3E1/4.0)+C6[i][j]*t4130*t4132*1/(t4134*t4134*t4134) \
                   *t4135*t4136*t4139*exp(d_*t4127*-2.0)*(1.0/2.0)-C6[i][j]*t4130*t4140*t4132*t4135*t4136*t4128*t4139* \
                   (1.0/4.0);
              double t4142 = x[i]-x[j];
              double t4143 = y[i]-y[j];
              double t4144 = z[i]-z[j];
              double t4145 = t4142*t4142;
              double t4146 = t4143*t4143;
              double t4147 = t4144*t4144;
              double t4148 = t4145+t4146+t4147;
              double t4158 = z[i]*2.0;
              double t4159 = z[j]*2.0;
              double t4149 = t4158-t4159;
              double t4150 = 1/RvdW[i][j];
              double t4151 = sqrt(t4148);
              double t4152 = t4150*t4151;
              double t4153 = t4152-1.0;
              double t4157 = d_*t4153;
              double t4154 = exp(-t4157);
              double t4155 = t4154+1.0;
              double t4156 = 1/t4155;
              double t4160 = t4149*t4149;
              double t4161 = 1/(t4155*t4155);
              double t4162 = 1/(t4148*t4148*t4148*t4148);
              double t4163 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4164 = d_*d_;
              hess[Az][Az] += C6[i][j]*t4162*t4156*-6.0+C6[i][j]*t4160*t4156*1/(t4148*t4148*t4148*t4148*t4148)*1.2E1+ \
                   C6[i][j]*d_*t4150*t4161*t4154*1/pow(t4148,7.0/2.0)-C6[i][j]*t4160*t4161*t4162*t4154*t4163*t4164*(1.0/ \
                   4.0)-C6[i][j]*d_*t4150*t4160*t4161*t4154*1/pow(t4148,9.0/2.0)*(1.3E1/4.0)+C6[i][j]*t4160*t4162*t4163* \
                   1/(t4155*t4155*t4155)*t4164*exp(d_*t4153*-2.0)*(1.0/2.0);
        }
        hess[Ay][Ax] = hess[Ax][Ay];
        hess[Az][Ax] = hess[Ax][Az];
        hess[Az][Ay] = hess[Ay][Az];

        //fprintf(outfile," A = %d B = %d\n",i,i);
        //print_mat(hess,3*mol->natom(),3*mol->natom(),outfile);
    }
    // Case 2: A B
    for (int i = 1; i < natom; i++) {
        Ax = 3*i;
        Ay = 3*i + 1;
        Az = 3*i + 2;

        for (int j = 0; j < i; j++) {
            if (i == j)
                continue;

            Bx = 3*j;
            By = 3*j + 1;
            Bz = 3*j + 2;

              double t4166 = x[i]-x[j];
              double t4167 = y[i]-y[j];
              double t4168 = z[i]-z[j];
              double t4169 = t4166*t4166;
              double t4170 = t4167*t4167;
              double t4171 = t4168*t4168;
              double t4172 = t4170+t4171+t4169;
              double t4182 = x[i]*2.0;
              double t4183 = x[j]*2.0;
              double t4173 = t4182-t4183;
              double t4174 = 1/RvdW[i][j];
              double t4175 = sqrt(t4172);
              double t4176 = t4174*t4175;
              double t4177 = t4176-1.0;
              double t4181 = d_*t4177;
              double t4178 = exp(-t4181);
              double t4179 = t4178+1.0;
              double t4180 = 1/t4179;
              double t4184 = t4173*t4173;
              double t4185 = 1/(t4179*t4179);
              double t4186 = 1/(t4172*t4172*t4172*t4172);
              double t4187 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4188 = d_*d_;
              hess[Ax][Bx] += C6[i][j]*t4180*t4186*6.0-C6[i][j]*t4180*1/(t4172*t4172*t4172*t4172*t4172)*t4184*1.2E1- \
                   C6[i][j]*d_*1/pow(t4172,7.0/2.0)*t4174*t4185*t4178+C6[i][j]*t4184*t4185*t4186*t4178*t4187*t4188*(1.0/ \
                   4.0)+C6[i][j]*d_*1/pow(t4172,9.0/2.0)*t4174*t4184*t4185*t4178*(1.3E1/4.0)-C6[i][j]*t4184*t4186*t4187* \
                   1/(t4179*t4179*t4179)*t4188*exp(d_*t4177*-2.0)*(1.0/2.0);
              double t4190 = x[i]-x[j];
              double t4191 = y[i]-y[j];
              double t4192 = z[i]-z[j];
              double t4193 = t4190*t4190;
              double t4194 = t4191*t4191;
              double t4195 = t4192*t4192;
              double t4196 = t4193+t4194+t4195;
              double t4197 = 1/RvdW[i][j];
              double t4198 = sqrt(t4196);
              double t4199 = t4197*t4198;
              double t4200 = t4199-1.0;
              double t4206 = d_*t4200;
              double t4201 = exp(-t4206);
              double t4202 = x[i]*2.0;
              double t4210 = x[j]*2.0;
              double t4203 = -t4210+t4202;
              double t4204 = y[i]*2.0;
              double t4211 = y[j]*2.0;
              double t4205 = -t4211+t4204;
              double t4207 = t4201+1.0;
              double t4208 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4209 = d_*d_;
              double t4212 = 1/(t4196*t4196*t4196*t4196);
              double t4213 = 1/(t4207*t4207);
              hess[Ax][By] += (C6[i][j]*1/(t4196*t4196*t4196*t4196*t4196)*t4203*t4205*-1.2E1)/t4207+C6[i][j]*d_*1/ \
                   pow(t4196,9.0/2.0)*t4197*t4201*t4203*t4213*t4205*(1.3E1/4.0)-C6[i][j]*t4203*t4212*t4205*1/(t4207*t4207* \
                   t4207)*t4208*t4209*exp(d_*t4200*-2.0)*(1.0/2.0)+C6[i][j]*t4201*t4203*t4212*t4213*t4205*t4208*t4209* \
                   (1.0/4.0);
              double t4215 = x[i]-x[j];
              double t4216 = y[i]-y[j];
              double t4217 = z[i]-z[j];
              double t4218 = t4215*t4215;
              double t4219 = t4216*t4216;
              double t4220 = t4217*t4217;
              double t4221 = t4220+t4218+t4219;
              double t4222 = 1/RvdW[i][j];
              double t4223 = sqrt(t4221);
              double t4224 = t4222*t4223;
              double t4225 = t4224-1.0;
              double t4231 = d_*t4225;
              double t4226 = exp(-t4231);
              double t4227 = x[i]*2.0;
              double t4235 = x[j]*2.0;
              double t4228 = -t4235+t4227;
              double t4229 = z[i]*2.0;
              double t4236 = z[j]*2.0;
              double t4230 = -t4236+t4229;
              double t4232 = t4226+1.0;
              double t4233 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4234 = d_*d_;
              double t4237 = 1/(t4221*t4221*t4221*t4221);
              double t4238 = 1/(t4232*t4232);
              hess[Ax][Bz] += (C6[i][j]*1/(t4221*t4221*t4221*t4221*t4221)*t4230*t4228*-1.2E1)/t4232+C6[i][j]*d_*1/ \
                   pow(t4221,9.0/2.0)*t4230*t4222*t4226*t4228*t4238*(1.3E1/4.0)-C6[i][j]*t4230*1/(t4232*t4232*t4232)* \
                   t4233*t4234*t4228*t4237*exp(d_*t4225*-2.0)*(1.0/2.0)+C6[i][j]*t4230*t4233*t4234*t4226*t4228*t4237* \
                   t4238*(1.0/4.0);
              double t4240 = x[i]-x[j];
              double t4241 = y[i]-y[j];
              double t4242 = z[i]-z[j];
              double t4243 = t4240*t4240;
              double t4244 = t4241*t4241;
              double t4245 = t4242*t4242;
              double t4246 = t4243+t4244+t4245;
              double t4247 = 1/RvdW[i][j];
              double t4248 = sqrt(t4246);
              double t4249 = t4247*t4248;
              double t4250 = t4249-1.0;
              double t4256 = d_*t4250;
              double t4251 = exp(-t4256);
              double t4252 = x[i]*2.0;
              double t4260 = x[j]*2.0;
              double t4253 = -t4260+t4252;
              double t4254 = y[i]*2.0;
              double t4261 = y[j]*2.0;
              double t4255 = -t4261+t4254;
              double t4257 = t4251+1.0;
              double t4258 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4259 = d_*d_;
              double t4262 = 1/(t4246*t4246*t4246*t4246);
              double t4263 = 1/(t4257*t4257);
              hess[Ay][Bx] += (C6[i][j]*t4253*1/(t4246*t4246*t4246*t4246*t4246)*t4255*-1.2E1)/t4257+C6[i][j]*d_*t4251* \
                   t4253*t4263*1/pow(t4246,9.0/2.0)*t4255*t4247*(1.3E1/4.0)-C6[i][j]*t4253*t4262*t4255*1/(t4257*t4257* \
                   t4257)*t4258*t4259*exp(d_*t4250*-2.0)*(1.0/2.0)+C6[i][j]*t4251*t4253*t4262*t4263*t4255*t4258*t4259* \
                   (1.0/4.0);
              double t4265 = x[i]-x[j];
              double t4266 = y[i]-y[j];
              double t4267 = z[i]-z[j];
              double t4268 = t4265*t4265;
              double t4269 = t4266*t4266;
              double t4270 = t4267*t4267;
              double t4271 = t4270+t4268+t4269;
              double t4281 = y[i]*2.0;
              double t4282 = y[j]*2.0;
              double t4272 = t4281-t4282;
              double t4273 = 1/RvdW[i][j];
              double t4274 = sqrt(t4271);
              double t4275 = t4273*t4274;
              double t4276 = t4275-1.0;
              double t4280 = d_*t4276;
              double t4277 = exp(-t4280);
              double t4278 = t4277+1.0;
              double t4279 = 1/t4278;
              double t4283 = t4272*t4272;
              double t4284 = 1/(t4278*t4278);
              double t4285 = 1/(t4271*t4271*t4271*t4271);
              double t4286 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4287 = d_*d_;
              hess[Ay][By] += C6[i][j]*t4285*t4279*6.0-C6[i][j]*1/(t4271*t4271*t4271*t4271*t4271)*t4283*t4279*1.2E1- \
                   C6[i][j]*d_*1/pow(t4271,7.0/2.0)*t4273*t4284*t4277+C6[i][j]*t4283*t4284*t4285*t4277*t4286*t4287*(1.0/ \
                   4.0)+C6[i][j]*d_*1/pow(t4271,9.0/2.0)*t4273*t4283*t4284*t4277*(1.3E1/4.0)-C6[i][j]*t4283*t4285*t4286* \
                   1/(t4278*t4278*t4278)*t4287*exp(d_*t4276*-2.0)*(1.0/2.0);
              double t4289 = x[i]-x[j];
              double t4290 = y[i]-y[j];
              double t4291 = z[i]-z[j];
              double t4292 = t4289*t4289;
              double t4293 = t4290*t4290;
              double t4294 = t4291*t4291;
              double t4295 = t4292+t4293+t4294;
              double t4296 = 1/RvdW[i][j];
              double t4297 = sqrt(t4295);
              double t4298 = t4296*t4297;
              double t4299 = t4298-1.0;
              double t4305 = d_*t4299;
              double t4300 = exp(-t4305);
              double t4301 = y[i]*2.0;
              double t4309 = y[j]*2.0;
              double t4302 = t4301-t4309;
              double t4303 = z[i]*2.0;
              double t4310 = z[j]*2.0;
              double t4304 = -t4310+t4303;
              double t4306 = t4300+1.0;
              double t4307 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4308 = d_*d_;
              double t4311 = 1/(t4295*t4295*t4295*t4295);
              double t4312 = 1/(t4306*t4306);
              hess[Ay][Bz] += (C6[i][j]*1/(t4295*t4295*t4295*t4295*t4295)*t4302*t4304*-1.2E1)/t4306+C6[i][j]*d_*1/ \
                   pow(t4295,9.0/2.0)*t4296*t4300*t4302*t4312*t4304*(1.3E1/4.0)-C6[i][j]*t4302*t4311*t4304*1/(t4306*t4306* \
                   t4306)*t4307*t4308*exp(d_*t4299*-2.0)*(1.0/2.0)+C6[i][j]*t4300*t4302*t4311*t4312*t4304*t4307*t4308* \
                   (1.0/4.0);
              double t4314 = x[i]-x[j];
              double t4315 = y[i]-y[j];
              double t4316 = z[i]-z[j];
              double t4317 = t4314*t4314;
              double t4318 = t4315*t4315;
              double t4319 = t4316*t4316;
              double t4320 = t4317+t4318+t4319;
              double t4321 = 1/RvdW[i][j];
              double t4322 = sqrt(t4320);
              double t4323 = t4321*t4322;
              double t4324 = t4323-1.0;
              double t4330 = d_*t4324;
              double t4325 = exp(-t4330);
              double t4326 = x[i]*2.0;
              double t4334 = x[j]*2.0;
              double t4327 = -t4334+t4326;
              double t4328 = z[i]*2.0;
              double t4335 = z[j]*2.0;
              double t4329 = -t4335+t4328;
              double t4331 = t4325+1.0;
              double t4332 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4333 = d_*d_;
              double t4336 = 1/(t4320*t4320*t4320*t4320);
              double t4337 = 1/(t4331*t4331);
              hess[Az][Bx] += (C6[i][j]*1/(t4320*t4320*t4320*t4320*t4320)*t4327*t4329*-1.2E1)/t4331+C6[i][j]*d_*1/ \
                   pow(t4320,9.0/2.0)*t4321*t4325*t4327*t4337*t4329*(1.3E1/4.0)-C6[i][j]*1/(t4331*t4331*t4331)*t4332* \
                   t4333*t4327*t4336*t4329*exp(d_*t4324*-2.0)*(1.0/2.0)+C6[i][j]*t4332*t4333*t4325*t4327*t4336*t4337* \
                   t4329*(1.0/4.0);
              double t4339 = x[i]-x[j];
              double t4340 = y[i]-y[j];
              double t4341 = z[i]-z[j];
              double t4342 = t4339*t4339;
              double t4343 = t4340*t4340;
              double t4344 = t4341*t4341;
              double t4345 = t4342+t4343+t4344;
              double t4346 = 1/RvdW[i][j];
              double t4347 = sqrt(t4345);
              double t4348 = t4346*t4347;
              double t4349 = t4348-1.0;
              double t4355 = d_*t4349;
              double t4350 = exp(-t4355);
              double t4351 = y[i]*2.0;
              double t4359 = y[j]*2.0;
              double t4352 = t4351-t4359;
              double t4353 = z[i]*2.0;
              double t4360 = z[j]*2.0;
              double t4354 = -t4360+t4353;
              double t4356 = t4350+1.0;
              double t4357 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4358 = d_*d_;
              double t4361 = 1/(t4345*t4345*t4345*t4345);
              double t4362 = 1/(t4356*t4356);
              hess[Az][By] += (C6[i][j]*t4352*1/(t4345*t4345*t4345*t4345*t4345)*t4354*-1.2E1)/t4356+C6[i][j]*d_*t4350* \
                   t4352*t4362*1/pow(t4345,9.0/2.0)*t4354*t4346*(1.3E1/4.0)-C6[i][j]*t4352*t4361*t4354*1/(t4356*t4356* \
                   t4356)*t4357*t4358*exp(d_*t4349*-2.0)*(1.0/2.0)+C6[i][j]*t4350*t4352*t4361*t4362*t4354*t4357*t4358* \
                   (1.0/4.0);
              double t4364 = x[i]-x[j];
              double t4365 = y[i]-y[j];
              double t4366 = z[i]-z[j];
              double t4367 = t4364*t4364;
              double t4368 = t4365*t4365;
              double t4369 = t4366*t4366;
              double t4370 = t4367+t4368+t4369;
              double t4380 = z[i]*2.0;
              double t4381 = z[j]*2.0;
              double t4371 = t4380-t4381;
              double t4372 = 1/RvdW[i][j];
              double t4373 = sqrt(t4370);
              double t4374 = t4372*t4373;
              double t4375 = t4374-1.0;
              double t4379 = d_*t4375;
              double t4376 = exp(-t4379);
              double t4377 = t4376+1.0;
              double t4378 = 1/t4377;
              double t4382 = t4371*t4371;
              double t4383 = 1/(t4377*t4377);
              double t4384 = 1/(t4370*t4370*t4370*t4370);
              double t4385 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4386 = d_*d_;
              hess[Az][Bz] += C6[i][j]*t4384*t4378*6.0-C6[i][j]*1/(t4370*t4370*t4370*t4370*t4370)*t4382*t4378*1.2E1- \
                   C6[i][j]*d_*1/pow(t4370,7.0/2.0)*t4372*t4383*t4376+C6[i][j]*t4382*t4383*t4384*t4376*t4385*t4386*(1.0/ \
                   4.0)+C6[i][j]*d_*1/pow(t4370,9.0/2.0)*t4372*t4382*t4383*t4376*(1.3E1/4.0)-C6[i][j]*t4382*t4384*t4385* \
                   1/(t4377*t4377*t4377)*t4386*exp(d_*t4375*-2.0)*(1.0/2.0);

            hess[Bx][Ax] = hess[Ax][Bx];
            hess[Bx][Ay] = hess[Ay][Bx];
            hess[Bx][Az] = hess[Az][Bx];
            hess[By][Ax] = hess[Ax][By];
            hess[By][Ay] = hess[Ay][By];
            hess[By][Az] = hess[Az][By];
            hess[Bz][Ax] = hess[Ax][Bz];
            hess[Bz][Ay] = hess[Ay][Bz];
            hess[Bz][Az] = hess[Az][Bz];
        }
    }

    // Scale gradient by -s6_
    // Dispersion energy is always stabilizing
    for (int A = 0; A < 3*natom; A++) {
        for (int B = 0; B < 3*natom; B++) {
            hess[A][B] *= -s6_;
        }
    }

    // Free stuff
    free_block(C6);
    free_block(RvdW);
    free_block(r);
    delete []Z;
    delete []x;
    delete []y;
    delete []z;
    return hess;
}
D2::D2(double s6) : Dispersion()
{
    name_ = "-D2";
    description_ = "Grimme's -D2 Dispersion Correction";
    citation_ = "Grimme, S. (2006),  J. Comp. Chem., 27: 1787-1799";

    d_ = 20.0;
    RvdW_ = RvdW_D2_;
    C6_ = C6_D2_;
    s6_ = s6;
}
D2::~D2()
{
}
double D2::computeEnergy(boost::shared_ptr<Molecule> mol)
{
    double energy = 0.0;

    // Build Z, x, y, and z
    int natom = 0;

    int *Z = init_int_array(mol->natom());
    double *x = init_array(mol->natom());
    double *y = init_array(mol->natom());
    double *z = init_array(mol->natom());

    for (int i = 0; i < mol->natom(); i++) {
        if (mol->Z(i) == 0)
        continue;
        Z[natom] = mol->Z(i);
        x[natom] = mol->x(i);
        y[natom] = mol->y(i);
        z[natom] = mol->z(i);
        natom ++;
    }

    // Build C6 (C6_ij = C6_i*C6_j/(C6_i+C6_j) in -D1 )
    double **C6 = block_matrix(natom, natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            C6[i][j] = sqrt( C6_[Z[i]] * C6_[Z[j]] );
            C6[j][i] = C6[i][j];
        }
    }

    // Build RvdW (sum of vdW radii)
    double **RvdW = block_matrix(natom, natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            RvdW[i][j] = RvdW_[Z[i]] + RvdW_[Z[j]];
            RvdW[j][i] = RvdW[i][j];
        }
    }

    // Build r
    double **r = block_matrix(natom,natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            r[i][j] = sqrt((x[i]-x[j])*(x[i]-x[j]) + \
                           (y[i]-y[j])*(y[i]-y[j]) + \
                           (z[i]-z[j])*(z[i]-z[j]));
            r[j][i] = r[i][j];

        }
    }

    // Compute energy
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
              double t3960 = x[i]-x[j];
              double t3961 = y[i]-y[j];
              double t3962 = z[i]-z[j];
              double t3963 = t3960*t3960;
              double t3964 = t3961*t3961;
              double t3965 = t3962*t3962;
              double t3966 = t3963+t3964+t3965;
              energy += (C6[i][j]*1/(t3966*t3966*t3966))/(exp(-d_*(sqrt(t3966)/RvdW[i][j]-1.0))+1.0);
        }
    }

    // Scale energy by s6_
    energy *= s6_;

    // Dispersion energy is always stabilizing
    energy *= -1.0;

    // Free stuff
    free_block(C6);
    free_block(RvdW);
    free_block(r);
    delete []Z;
    delete []x;
    delete []y;
    delete []z;

    return energy;
}
double* D2::computeGradient(boost::shared_ptr<Molecule> mol)
{
    // Build Z, x, y, and z
    int natom = 0;

    int *Z = init_int_array(mol->natom());
    double *x = init_array(mol->natom());
    double *y = init_array(mol->natom());
    double *z = init_array(mol->natom());

    for (int i = 0; i < mol->natom(); i++) {
        if (mol->Z(i) == 0)
        continue;
        Z[natom] = mol->Z(i);
        x[natom] = mol->x(i);
        y[natom] = mol->y(i);
        z[natom] = mol->z(i);
        natom ++;
    }

    // gradient array in [x1 y1 z1 x2 y2 z2 ... ]
    double* grad = init_array(3*natom);

    // Build C6 (C6_ij = C6_i*C6_j/(C6_i+C6_j) in -D1 )
    double **C6 = block_matrix(natom, natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            C6[i][j] = sqrt( C6_[Z[i]] * C6_[Z[j]] );
            C6[j][i] = C6[i][j];
        }
    }

    // Build RvdW (sum of vdW radii)
    double **RvdW = block_matrix(natom, natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            RvdW[i][j] = RvdW_[Z[i]] + RvdW_[Z[j]];
            RvdW[j][i] = RvdW[i][j];
        }
    }

    // Build r
    double **r = block_matrix(natom,natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            r[i][j] = sqrt((x[i]-x[j])*(x[i]-x[j]) + \
                           (y[i]-y[j])*(y[i]-y[j]) + \
                           (z[i]-z[j])*(z[i]-z[j]));
            r[j][i] = r[i][j];

        }
    }

    // Compute gradient by perturbing nucleus i
    int Ax, Ay, Az;
    double f;
    for (int i = 0; i < natom; i++) {
        Ax = 3*i;
        Ay = 3*i + 1;
        Az = 3*i + 2;

        for (int j = 0; j < natom; j++) {
            if (i == j)
                continue;
              double t3968 = x[i]-x[j];
              double t3969 = y[i]-y[j];
              double t3970 = z[i]-z[j];
              double t3971 = t3968*t3968;
              double t3972 = t3969*t3969;
              double t3973 = t3970*t3970;
              double t3974 = t3971+t3972+t3973;
              double t3975 = 1/RvdW[i][j];
              double t3976 = sqrt(t3974);
              double t3977 = t3975*t3976;
              double t3978 = t3977-1.0;
              double t3982 = d_*t3978;
              double t3979 = exp(-t3982);
              double t3980 = x[i]*2.0;
              double t3981 = t3980-x[j]*2.0;
              double t3983 = t3979+1.0;
              grad[Ax] += (C6[i][j]*t3981*1/(t3974*t3974*t3974*t3974)*-3.0)/t3983+C6[i][j]*d_*t3981*1/pow(t3974,7.0/ \
                   2.0)*1/(t3983*t3983)*t3975*t3979*(1.0/2.0);
              double t3985 = x[i]-x[j];
              double t3986 = y[i]-y[j];
              double t3987 = z[i]-z[j];
              double t3988 = t3985*t3985;
              double t3989 = t3986*t3986;
              double t3990 = t3987*t3987;
              double t3991 = t3990+t3988+t3989;
              double t3992 = 1/RvdW[i][j];
              double t3993 = sqrt(t3991);
              double t3994 = t3992*t3993;
              double t3995 = t3994-1.0;
              double t3999 = d_*t3995;
              double t3996 = exp(-t3999);
              double t3997 = y[i]*2.0;
              double t3998 = t3997-y[j]*2.0;
              double t4000 = t3996+1.0;
              grad[Ay] += (C6[i][j]*1/(t3991*t3991*t3991*t3991)*t3998*-3.0)/t4000+C6[i][j]*d_*1/pow(t3991,7.0/2.0) \
                   *t3992*t3996*t3998*1/(t4000*t4000)*(1.0/2.0);
              double t4002 = x[i]-x[j];
              double t4003 = y[i]-y[j];
              double t4004 = z[i]-z[j];
              double t4005 = t4002*t4002;
              double t4006 = t4003*t4003;
              double t4007 = t4004*t4004;
              double t4008 = t4005+t4006+t4007;
              double t4009 = 1/RvdW[i][j];
              double t4010 = sqrt(t4008);
              double t4011 = t4010*t4009;
              double t4012 = t4011-1.0;
              double t4016 = d_*t4012;
              double t4013 = exp(-t4016);
              double t4014 = z[i]*2.0;
              double t4015 = t4014-z[j]*2.0;
              double t4017 = t4013+1.0;
              grad[Az] += (C6[i][j]*t4015*1/(t4008*t4008*t4008*t4008)*-3.0)/t4017+C6[i][j]*d_*t4013*t4015*1/pow(t4008,7.0/ \
                   2.0)*1/(t4017*t4017)*t4009*(1.0/2.0);
        }
    }

    // Scale gradient by -s6_
    // Dispersion energy is always stabilizing
    for (int A = 0; A < 3*natom; A++) {
        grad[A] *= -s6_;
    }

    // Free stuff
    free_block(C6);
    free_block(RvdW);
    free_block(r);
    delete []Z;
    delete []x;
    delete []y;
    delete []z;
    return grad;
}
double** D2::computeHessian(boost::shared_ptr<Molecule> mol)
{
    // Build Z, x, y, and z
    int natom = 0;

    int *Z = init_int_array(mol->natom());
    double *x = init_array(mol->natom());
    double *y = init_array(mol->natom());
    double *z = init_array(mol->natom());

    for (int i = 0; i < mol->natom(); i++) {
        if (mol->Z(i) == 0)
        continue;
        Z[natom] = mol->Z(i);
        x[natom] = mol->x(i);
        y[natom] = mol->y(i);
        z[natom] = mol->z(i);
        natom ++;
    }

    // gradient array in [x1 y1 z1 x2 y2 z2 ... ]
    double** hess = block_matrix(3*natom,3*natom);

    // Build C6 (C6_ij = C6_i*C6_j/(C6_i+C6_j) in -D1 )
    double **C6 = block_matrix(natom, natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            C6[i][j] = sqrt( C6_[Z[i]] * C6_[Z[j]] );
            C6[j][i] = C6[i][j];
        }
    }

    // Build RvdW (sum of vdW radii)
    double **RvdW = block_matrix(natom, natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            RvdW[i][j] = RvdW_[Z[i]] + RvdW_[Z[j]];
            RvdW[j][i] = RvdW[i][j];
        }
    }

    // Build r
    double **r = block_matrix(natom,natom);
    for (int i = 1; i < natom; i++) {
        for (int j = 0; j < i; j++) {
            r[i][j] = sqrt((x[i]-x[j])*(x[i]-x[j]) + \
                           (y[i]-y[j])*(y[i]-y[j]) + \
                           (z[i]-z[j])*(z[i]-z[j]));
            r[j][i] = r[i][j];

        }
    }

    int Ax, Ay, Az, Bx, By, Bz;
    double f;

    // Case 1: A A
    for (int i = 0; i < natom; i++) {
        Ax = 3*i;
        Ay = 3*i + 1;
        Az = 3*i + 2;

        for (int j = 0; j < natom; j++) {
            if (i == j)
                continue;

              double t4019 = x[i]-x[j];
              double t4020 = y[i]-y[j];
              double t4021 = z[i]-z[j];
              double t4022 = t4019*t4019;
              double t4023 = t4020*t4020;
              double t4024 = t4021*t4021;
              double t4025 = t4022+t4023+t4024;
              double t4035 = x[i]*2.0;
              double t4036 = x[j]*2.0;
              double t4026 = t4035-t4036;
              double t4027 = 1/RvdW[i][j];
              double t4028 = sqrt(t4025);
              double t4029 = t4027*t4028;
              double t4030 = t4029-1.0;
              double t4034 = d_*t4030;
              double t4031 = exp(-t4034);
              double t4032 = t4031+1.0;
              double t4033 = 1/t4032;
              double t4037 = t4026*t4026;
              double t4038 = 1/(t4032*t4032);
              double t4039 = 1/(t4025*t4025*t4025*t4025);
              double t4040 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4041 = d_*d_;
              hess[Ax][Ax] += C6[i][j]*t4033*t4039*-6.0+C6[i][j]*t4033*1/(t4025*t4025*t4025*t4025*t4025)*t4037*1.2E1+ \
                   C6[i][j]*d_*t4031*1/pow(t4025,7.0/2.0)*t4027*t4038-C6[i][j]*t4031*t4040*t4041*t4037*t4038*t4039*(1.0/ \
                   4.0)-C6[i][j]*d_*t4031*1/pow(t4025,9.0/2.0)*t4027*t4037*t4038*(1.3E1/4.0)+C6[i][j]*t4040*1/(t4032* \
                   t4032*t4032)*t4041*t4037*t4039*exp(d_*t4030*-2.0)*(1.0/2.0);
              double t4043 = x[i]-x[j];
              double t4044 = y[i]-y[j];
              double t4045 = z[i]-z[j];
              double t4046 = t4043*t4043;
              double t4047 = t4044*t4044;
              double t4048 = t4045*t4045;
              double t4049 = t4046+t4047+t4048;
              double t4050 = 1/RvdW[i][j];
              double t4051 = sqrt(t4049);
              double t4052 = t4050*t4051;
              double t4053 = t4052-1.0;
              double t4059 = d_*t4053;
              double t4054 = exp(-t4059);
              double t4055 = x[i]*2.0;
              double t4063 = x[j]*2.0;
              double t4056 = -t4063+t4055;
              double t4057 = y[i]*2.0;
              double t4064 = y[j]*2.0;
              double t4058 = -t4064+t4057;
              double t4060 = t4054+1.0;
              double t4061 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4062 = d_*d_;
              double t4065 = 1/(t4049*t4049*t4049*t4049);
              double t4066 = 1/(t4060*t4060);
              hess[Ax][Ay] += (C6[i][j]*t4056*1/(t4049*t4049*t4049*t4049*t4049)*t4058*1.2E1)/t4060-C6[i][j]*d_*t4050* \
                   t4054*t4056*t4066*1/pow(t4049,9.0/2.0)*t4058*(1.3E1/4.0)+C6[i][j]*1/(t4060*t4060*t4060)*t4061*t4062* \
                   t4056*t4065*t4058*exp(d_*t4053*-2.0)*(1.0/2.0)-C6[i][j]*t4061*t4062*t4054*t4056*t4065*t4066*t4058* \
                   (1.0/4.0);
              double t4068 = x[i]-x[j];
              double t4069 = y[i]-y[j];
              double t4070 = z[i]-z[j];
              double t4071 = t4068*t4068;
              double t4072 = t4069*t4069;
              double t4073 = t4070*t4070;
              double t4074 = t4071+t4072+t4073;
              double t4075 = 1/RvdW[i][j];
              double t4076 = sqrt(t4074);
              double t4077 = t4075*t4076;
              double t4078 = t4077-1.0;
              double t4084 = d_*t4078;
              double t4079 = exp(-t4084);
              double t4080 = x[i]*2.0;
              double t4088 = x[j]*2.0;
              double t4081 = t4080-t4088;
              double t4082 = z[i]*2.0;
              double t4089 = z[j]*2.0;
              double t4083 = t4082-t4089;
              double t4085 = t4079+1.0;
              double t4086 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4087 = d_*d_;
              double t4090 = 1/(t4074*t4074*t4074*t4074);
              double t4091 = 1/(t4085*t4085);
              hess[Ax][Az] += (C6[i][j]*t4081*1/(t4074*t4074*t4074*t4074*t4074)*t4083*1.2E1)/t4085-C6[i][j]*d_*t4081* \
                   t4091*1/pow(t4074,9.0/2.0)*t4083*t4075*t4079*(1.3E1/4.0)+C6[i][j]*t4081*t4090*t4083*1/(t4085*t4085* \
                   t4085)*t4086*t4087*exp(d_*t4078*-2.0)*(1.0/2.0)-C6[i][j]*t4081*t4090*t4091*t4083*t4086*t4087*t4079* \
                   (1.0/4.0);
              double t4093 = x[i]-x[j];
              double t4094 = y[i]-y[j];
              double t4095 = z[i]-z[j];
              double t4096 = t4093*t4093;
              double t4097 = t4094*t4094;
              double t4098 = t4095*t4095;
              double t4099 = t4096+t4097+t4098;
              double t4109 = y[i]*2.0;
              double t4110 = y[j]*2.0;
              double t4100 = -t4110+t4109;
              double t4101 = 1/RvdW[i][j];
              double t4102 = sqrt(t4099);
              double t4103 = t4101*t4102;
              double t4104 = t4103-1.0;
              double t4108 = d_*t4104;
              double t4105 = exp(-t4108);
              double t4106 = t4105+1.0;
              double t4107 = 1/t4106;
              double t4111 = t4100*t4100;
              double t4112 = 1/(t4106*t4106);
              double t4113 = 1/(t4099*t4099*t4099*t4099);
              double t4114 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4115 = d_*d_;
              hess[Ay][Ay] += C6[i][j]*t4113*t4107*-6.0+C6[i][j]*1/(t4099*t4099*t4099*t4099*t4099)*t4111*t4107*1.2E1+ \
                   C6[i][j]*d_*1/pow(t4099,7.0/2.0)*t4101*t4112*t4105-C6[i][j]*t4111*t4112*t4113*t4105*t4114*t4115*(1.0/ \
                   4.0)-C6[i][j]*d_*1/pow(t4099,9.0/2.0)*t4101*t4111*t4112*t4105*(1.3E1/4.0)+C6[i][j]*t4111*t4113*t4114* \
                   1/(t4106*t4106*t4106)*t4115*exp(d_*t4104*-2.0)*(1.0/2.0);
              double t4117 = x[i]-x[j];
              double t4118 = y[i]-y[j];
              double t4119 = z[i]-z[j];
              double t4120 = t4117*t4117;
              double t4121 = t4118*t4118;
              double t4122 = t4119*t4119;
              double t4123 = t4120+t4121+t4122;
              double t4124 = 1/RvdW[i][j];
              double t4125 = sqrt(t4123);
              double t4126 = t4124*t4125;
              double t4127 = t4126-1.0;
              double t4133 = d_*t4127;
              double t4128 = exp(-t4133);
              double t4129 = y[i]*2.0;
              double t4137 = y[j]*2.0;
              double t4130 = -t4137+t4129;
              double t4131 = z[i]*2.0;
              double t4138 = z[j]*2.0;
              double t4132 = t4131-t4138;
              double t4134 = t4128+1.0;
              double t4135 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4136 = d_*d_;
              double t4139 = 1/(t4123*t4123*t4123*t4123);
              double t4140 = 1/(t4134*t4134);
              hess[Ay][Az] += (C6[i][j]*t4130*1/(t4123*t4123*t4123*t4123*t4123)*t4132*1.2E1)/t4134-C6[i][j]*d_*t4130* \
                   t4140*1/pow(t4123,9.0/2.0)*t4132*t4124*t4128*(1.3E1/4.0)+C6[i][j]*t4130*t4132*1/(t4134*t4134*t4134) \
                   *t4135*t4136*t4139*exp(d_*t4127*-2.0)*(1.0/2.0)-C6[i][j]*t4130*t4140*t4132*t4135*t4136*t4128*t4139* \
                   (1.0/4.0);
              double t4142 = x[i]-x[j];
              double t4143 = y[i]-y[j];
              double t4144 = z[i]-z[j];
              double t4145 = t4142*t4142;
              double t4146 = t4143*t4143;
              double t4147 = t4144*t4144;
              double t4148 = t4145+t4146+t4147;
              double t4158 = z[i]*2.0;
              double t4159 = z[j]*2.0;
              double t4149 = t4158-t4159;
              double t4150 = 1/RvdW[i][j];
              double t4151 = sqrt(t4148);
              double t4152 = t4150*t4151;
              double t4153 = t4152-1.0;
              double t4157 = d_*t4153;
              double t4154 = exp(-t4157);
              double t4155 = t4154+1.0;
              double t4156 = 1/t4155;
              double t4160 = t4149*t4149;
              double t4161 = 1/(t4155*t4155);
              double t4162 = 1/(t4148*t4148*t4148*t4148);
              double t4163 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4164 = d_*d_;
              hess[Az][Az] += C6[i][j]*t4162*t4156*-6.0+C6[i][j]*t4160*t4156*1/(t4148*t4148*t4148*t4148*t4148)*1.2E1+ \
                   C6[i][j]*d_*t4150*t4161*t4154*1/pow(t4148,7.0/2.0)-C6[i][j]*t4160*t4161*t4162*t4154*t4163*t4164*(1.0/ \
                   4.0)-C6[i][j]*d_*t4150*t4160*t4161*t4154*1/pow(t4148,9.0/2.0)*(1.3E1/4.0)+C6[i][j]*t4160*t4162*t4163* \
                   1/(t4155*t4155*t4155)*t4164*exp(d_*t4153*-2.0)*(1.0/2.0);
        }
        hess[Ay][Ax] = hess[Ax][Ay];
        hess[Az][Ax] = hess[Ax][Az];
        hess[Az][Ay] = hess[Ay][Az];

        //fprintf(outfile," A = %d B = %d\n",i,i);
        //print_mat(hess,3*mol->natom(),3*mol->natom(),outfile);
    }
    // Case 2: A B
    for (int i = 1; i < natom; i++) {
        Ax = 3*i;
        Ay = 3*i + 1;
        Az = 3*i + 2;

        for (int j = 0; j < i; j++) {
            if (i == j)
                continue;

            Bx = 3*j;
            By = 3*j + 1;
            Bz = 3*j + 2;

              double t4166 = x[i]-x[j];
              double t4167 = y[i]-y[j];
              double t4168 = z[i]-z[j];
              double t4169 = t4166*t4166;
              double t4170 = t4167*t4167;
              double t4171 = t4168*t4168;
              double t4172 = t4170+t4171+t4169;
              double t4182 = x[i]*2.0;
              double t4183 = x[j]*2.0;
              double t4173 = t4182-t4183;
              double t4174 = 1/RvdW[i][j];
              double t4175 = sqrt(t4172);
              double t4176 = t4174*t4175;
              double t4177 = t4176-1.0;
              double t4181 = d_*t4177;
              double t4178 = exp(-t4181);
              double t4179 = t4178+1.0;
              double t4180 = 1/t4179;
              double t4184 = t4173*t4173;
              double t4185 = 1/(t4179*t4179);
              double t4186 = 1/(t4172*t4172*t4172*t4172);
              double t4187 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4188 = d_*d_;
              hess[Ax][Bx] += C6[i][j]*t4180*t4186*6.0-C6[i][j]*t4180*1/(t4172*t4172*t4172*t4172*t4172)*t4184*1.2E1- \
                   C6[i][j]*d_*1/pow(t4172,7.0/2.0)*t4174*t4185*t4178+C6[i][j]*t4184*t4185*t4186*t4178*t4187*t4188*(1.0/ \
                   4.0)+C6[i][j]*d_*1/pow(t4172,9.0/2.0)*t4174*t4184*t4185*t4178*(1.3E1/4.0)-C6[i][j]*t4184*t4186*t4187* \
                   1/(t4179*t4179*t4179)*t4188*exp(d_*t4177*-2.0)*(1.0/2.0);
              double t4190 = x[i]-x[j];
              double t4191 = y[i]-y[j];
              double t4192 = z[i]-z[j];
              double t4193 = t4190*t4190;
              double t4194 = t4191*t4191;
              double t4195 = t4192*t4192;
              double t4196 = t4193+t4194+t4195;
              double t4197 = 1/RvdW[i][j];
              double t4198 = sqrt(t4196);
              double t4199 = t4197*t4198;
              double t4200 = t4199-1.0;
              double t4206 = d_*t4200;
              double t4201 = exp(-t4206);
              double t4202 = x[i]*2.0;
              double t4210 = x[j]*2.0;
              double t4203 = -t4210+t4202;
              double t4204 = y[i]*2.0;
              double t4211 = y[j]*2.0;
              double t4205 = -t4211+t4204;
              double t4207 = t4201+1.0;
              double t4208 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4209 = d_*d_;
              double t4212 = 1/(t4196*t4196*t4196*t4196);
              double t4213 = 1/(t4207*t4207);
              hess[Ax][By] += (C6[i][j]*1/(t4196*t4196*t4196*t4196*t4196)*t4203*t4205*-1.2E1)/t4207+C6[i][j]*d_*1/ \
                   pow(t4196,9.0/2.0)*t4197*t4201*t4203*t4213*t4205*(1.3E1/4.0)-C6[i][j]*t4203*t4212*t4205*1/(t4207*t4207* \
                   t4207)*t4208*t4209*exp(d_*t4200*-2.0)*(1.0/2.0)+C6[i][j]*t4201*t4203*t4212*t4213*t4205*t4208*t4209* \
                   (1.0/4.0);
              double t4215 = x[i]-x[j];
              double t4216 = y[i]-y[j];
              double t4217 = z[i]-z[j];
              double t4218 = t4215*t4215;
              double t4219 = t4216*t4216;
              double t4220 = t4217*t4217;
              double t4221 = t4220+t4218+t4219;
              double t4222 = 1/RvdW[i][j];
              double t4223 = sqrt(t4221);
              double t4224 = t4222*t4223;
              double t4225 = t4224-1.0;
              double t4231 = d_*t4225;
              double t4226 = exp(-t4231);
              double t4227 = x[i]*2.0;
              double t4235 = x[j]*2.0;
              double t4228 = -t4235+t4227;
              double t4229 = z[i]*2.0;
              double t4236 = z[j]*2.0;
              double t4230 = -t4236+t4229;
              double t4232 = t4226+1.0;
              double t4233 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4234 = d_*d_;
              double t4237 = 1/(t4221*t4221*t4221*t4221);
              double t4238 = 1/(t4232*t4232);
              hess[Ax][Bz] += (C6[i][j]*1/(t4221*t4221*t4221*t4221*t4221)*t4230*t4228*-1.2E1)/t4232+C6[i][j]*d_*1/ \
                   pow(t4221,9.0/2.0)*t4230*t4222*t4226*t4228*t4238*(1.3E1/4.0)-C6[i][j]*t4230*1/(t4232*t4232*t4232)* \
                   t4233*t4234*t4228*t4237*exp(d_*t4225*-2.0)*(1.0/2.0)+C6[i][j]*t4230*t4233*t4234*t4226*t4228*t4237* \
                   t4238*(1.0/4.0);
              double t4240 = x[i]-x[j];
              double t4241 = y[i]-y[j];
              double t4242 = z[i]-z[j];
              double t4243 = t4240*t4240;
              double t4244 = t4241*t4241;
              double t4245 = t4242*t4242;
              double t4246 = t4243+t4244+t4245;
              double t4247 = 1/RvdW[i][j];
              double t4248 = sqrt(t4246);
              double t4249 = t4247*t4248;
              double t4250 = t4249-1.0;
              double t4256 = d_*t4250;
              double t4251 = exp(-t4256);
              double t4252 = x[i]*2.0;
              double t4260 = x[j]*2.0;
              double t4253 = -t4260+t4252;
              double t4254 = y[i]*2.0;
              double t4261 = y[j]*2.0;
              double t4255 = -t4261+t4254;
              double t4257 = t4251+1.0;
              double t4258 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4259 = d_*d_;
              double t4262 = 1/(t4246*t4246*t4246*t4246);
              double t4263 = 1/(t4257*t4257);
              hess[Ay][Bx] += (C6[i][j]*t4253*1/(t4246*t4246*t4246*t4246*t4246)*t4255*-1.2E1)/t4257+C6[i][j]*d_*t4251* \
                   t4253*t4263*1/pow(t4246,9.0/2.0)*t4255*t4247*(1.3E1/4.0)-C6[i][j]*t4253*t4262*t4255*1/(t4257*t4257* \
                   t4257)*t4258*t4259*exp(d_*t4250*-2.0)*(1.0/2.0)+C6[i][j]*t4251*t4253*t4262*t4263*t4255*t4258*t4259* \
                   (1.0/4.0);
              double t4265 = x[i]-x[j];
              double t4266 = y[i]-y[j];
              double t4267 = z[i]-z[j];
              double t4268 = t4265*t4265;
              double t4269 = t4266*t4266;
              double t4270 = t4267*t4267;
              double t4271 = t4270+t4268+t4269;
              double t4281 = y[i]*2.0;
              double t4282 = y[j]*2.0;
              double t4272 = t4281-t4282;
              double t4273 = 1/RvdW[i][j];
              double t4274 = sqrt(t4271);
              double t4275 = t4273*t4274;
              double t4276 = t4275-1.0;
              double t4280 = d_*t4276;
              double t4277 = exp(-t4280);
              double t4278 = t4277+1.0;
              double t4279 = 1/t4278;
              double t4283 = t4272*t4272;
              double t4284 = 1/(t4278*t4278);
              double t4285 = 1/(t4271*t4271*t4271*t4271);
              double t4286 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4287 = d_*d_;
              hess[Ay][By] += C6[i][j]*t4285*t4279*6.0-C6[i][j]*1/(t4271*t4271*t4271*t4271*t4271)*t4283*t4279*1.2E1- \
                   C6[i][j]*d_*1/pow(t4271,7.0/2.0)*t4273*t4284*t4277+C6[i][j]*t4283*t4284*t4285*t4277*t4286*t4287*(1.0/ \
                   4.0)+C6[i][j]*d_*1/pow(t4271,9.0/2.0)*t4273*t4283*t4284*t4277*(1.3E1/4.0)-C6[i][j]*t4283*t4285*t4286* \
                   1/(t4278*t4278*t4278)*t4287*exp(d_*t4276*-2.0)*(1.0/2.0);
              double t4289 = x[i]-x[j];
              double t4290 = y[i]-y[j];
              double t4291 = z[i]-z[j];
              double t4292 = t4289*t4289;
              double t4293 = t4290*t4290;
              double t4294 = t4291*t4291;
              double t4295 = t4292+t4293+t4294;
              double t4296 = 1/RvdW[i][j];
              double t4297 = sqrt(t4295);
              double t4298 = t4296*t4297;
              double t4299 = t4298-1.0;
              double t4305 = d_*t4299;
              double t4300 = exp(-t4305);
              double t4301 = y[i]*2.0;
              double t4309 = y[j]*2.0;
              double t4302 = t4301-t4309;
              double t4303 = z[i]*2.0;
              double t4310 = z[j]*2.0;
              double t4304 = -t4310+t4303;
              double t4306 = t4300+1.0;
              double t4307 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4308 = d_*d_;
              double t4311 = 1/(t4295*t4295*t4295*t4295);
              double t4312 = 1/(t4306*t4306);
              hess[Ay][Bz] += (C6[i][j]*1/(t4295*t4295*t4295*t4295*t4295)*t4302*t4304*-1.2E1)/t4306+C6[i][j]*d_*1/ \
                   pow(t4295,9.0/2.0)*t4296*t4300*t4302*t4312*t4304*(1.3E1/4.0)-C6[i][j]*t4302*t4311*t4304*1/(t4306*t4306* \
                   t4306)*t4307*t4308*exp(d_*t4299*-2.0)*(1.0/2.0)+C6[i][j]*t4300*t4302*t4311*t4312*t4304*t4307*t4308* \
                   (1.0/4.0);
              double t4314 = x[i]-x[j];
              double t4315 = y[i]-y[j];
              double t4316 = z[i]-z[j];
              double t4317 = t4314*t4314;
              double t4318 = t4315*t4315;
              double t4319 = t4316*t4316;
              double t4320 = t4317+t4318+t4319;
              double t4321 = 1/RvdW[i][j];
              double t4322 = sqrt(t4320);
              double t4323 = t4321*t4322;
              double t4324 = t4323-1.0;
              double t4330 = d_*t4324;
              double t4325 = exp(-t4330);
              double t4326 = x[i]*2.0;
              double t4334 = x[j]*2.0;
              double t4327 = -t4334+t4326;
              double t4328 = z[i]*2.0;
              double t4335 = z[j]*2.0;
              double t4329 = -t4335+t4328;
              double t4331 = t4325+1.0;
              double t4332 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4333 = d_*d_;
              double t4336 = 1/(t4320*t4320*t4320*t4320);
              double t4337 = 1/(t4331*t4331);
              hess[Az][Bx] += (C6[i][j]*1/(t4320*t4320*t4320*t4320*t4320)*t4327*t4329*-1.2E1)/t4331+C6[i][j]*d_*1/ \
                   pow(t4320,9.0/2.0)*t4321*t4325*t4327*t4337*t4329*(1.3E1/4.0)-C6[i][j]*1/(t4331*t4331*t4331)*t4332* \
                   t4333*t4327*t4336*t4329*exp(d_*t4324*-2.0)*(1.0/2.0)+C6[i][j]*t4332*t4333*t4325*t4327*t4336*t4337* \
                   t4329*(1.0/4.0);
              double t4339 = x[i]-x[j];
              double t4340 = y[i]-y[j];
              double t4341 = z[i]-z[j];
              double t4342 = t4339*t4339;
              double t4343 = t4340*t4340;
              double t4344 = t4341*t4341;
              double t4345 = t4342+t4343+t4344;
              double t4346 = 1/RvdW[i][j];
              double t4347 = sqrt(t4345);
              double t4348 = t4346*t4347;
              double t4349 = t4348-1.0;
              double t4355 = d_*t4349;
              double t4350 = exp(-t4355);
              double t4351 = y[i]*2.0;
              double t4359 = y[j]*2.0;
              double t4352 = t4351-t4359;
              double t4353 = z[i]*2.0;
              double t4360 = z[j]*2.0;
              double t4354 = -t4360+t4353;
              double t4356 = t4350+1.0;
              double t4357 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4358 = d_*d_;
              double t4361 = 1/(t4345*t4345*t4345*t4345);
              double t4362 = 1/(t4356*t4356);
              hess[Az][By] += (C6[i][j]*t4352*1/(t4345*t4345*t4345*t4345*t4345)*t4354*-1.2E1)/t4356+C6[i][j]*d_*t4350* \
                   t4352*t4362*1/pow(t4345,9.0/2.0)*t4354*t4346*(1.3E1/4.0)-C6[i][j]*t4352*t4361*t4354*1/(t4356*t4356* \
                   t4356)*t4357*t4358*exp(d_*t4349*-2.0)*(1.0/2.0)+C6[i][j]*t4350*t4352*t4361*t4362*t4354*t4357*t4358* \
                   (1.0/4.0);
              double t4364 = x[i]-x[j];
              double t4365 = y[i]-y[j];
              double t4366 = z[i]-z[j];
              double t4367 = t4364*t4364;
              double t4368 = t4365*t4365;
              double t4369 = t4366*t4366;
              double t4370 = t4367+t4368+t4369;
              double t4380 = z[i]*2.0;
              double t4381 = z[j]*2.0;
              double t4371 = t4380-t4381;
              double t4372 = 1/RvdW[i][j];
              double t4373 = sqrt(t4370);
              double t4374 = t4372*t4373;
              double t4375 = t4374-1.0;
              double t4379 = d_*t4375;
              double t4376 = exp(-t4379);
              double t4377 = t4376+1.0;
              double t4378 = 1/t4377;
              double t4382 = t4371*t4371;
              double t4383 = 1/(t4377*t4377);
              double t4384 = 1/(t4370*t4370*t4370*t4370);
              double t4385 = 1/(RvdW[i][j]*RvdW[i][j]);
              double t4386 = d_*d_;
              hess[Az][Bz] += C6[i][j]*t4384*t4378*6.0-C6[i][j]*1/(t4370*t4370*t4370*t4370*t4370)*t4382*t4378*1.2E1- \
                   C6[i][j]*d_*1/pow(t4370,7.0/2.0)*t4372*t4383*t4376+C6[i][j]*t4382*t4383*t4384*t4376*t4385*t4386*(1.0/ \
                   4.0)+C6[i][j]*d_*1/pow(t4370,9.0/2.0)*t4372*t4382*t4383*t4376*(1.3E1/4.0)-C6[i][j]*t4382*t4384*t4385* \
                   1/(t4377*t4377*t4377)*t4386*exp(d_*t4375*-2.0)*(1.0/2.0);

            hess[Bx][Ax] = hess[Ax][Bx];
            hess[Bx][Ay] = hess[Ay][Bx];
            hess[Bx][Az] = hess[Az][Bx];
            hess[By][Ax] = hess[Ax][By];
            hess[By][Ay] = hess[Ay][By];
            hess[By][Az] = hess[Az][By];
            hess[Bz][Ax] = hess[Ax][Bz];
            hess[Bz][Ay] = hess[Ay][Bz];
            hess[Bz][Az] = hess[Az][Bz];
        }
    }

    // Scale gradient by -s6_
    // Dispersion energy is always stabilizing
    for (int A = 0; A < 3*natom; A++) {
        for (int B = 0; B < 3*natom; B++) {
            hess[A][B] *= -s6_;
        }
    }

    // Free stuff
    free_block(C6);
    free_block(RvdW);
    free_block(r);
    delete []Z;
    delete []x;
    delete []y;
    delete []z;
    return hess;
}
D3::D3(double s6, double s8, double sr6, double sr8)
{
    name_ = "-D3";
    description_ = "Grimme's -D3 Dispersion Correction";
    citation_ = "Grimme, S. (2010),  J.C.P., 132: 154104";

    s6_ = s6;
    s8_ = s8;
    sr6_ = sr6;
    sr8_ = sr8;
    alpha6_ = 14.0;
    alpha6_ = 16.0;
    k1_ = 16.0;
    k2_ = 4.0/3.0;
    k3_ = 4.0;
}
D3::~D3()
{
}
double D3::computeEnergy(boost::shared_ptr<Molecule> mol)
{
    double energy = 0.0;
    return energy;
}
double* D3::computeGradient(boost::shared_ptr<Molecule> mol)
{
    double* grad;
    return grad;
}
double** D3::computeHessian(boost::shared_ptr<Molecule> mol)
{
    double** hess;
    return hess;
}

}}
