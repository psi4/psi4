/**********************************************************
* dispersion.cc: definitions for -D for KS-DFT
* Robert Parrish, robparrish@gmail.com
* 09/01/2010
*
***********************************************************/

#include <libmints/matrix.h>
#include <libmints/molecule.h>
#include <libciomr/libciomr.h>
#include "dispersion.h"
#include "dispersion_defines.h"
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>

using namespace boost;
using namespace std;

namespace psi {

Dispersion::Dispersion()
{
}
Dispersion::~Dispersion()
{
}
boost::shared_ptr<Dispersion> Dispersion::build(const std::string & name, double s6)
{
    if (boost::to_upper_copy(name) == "-D1") {
        boost::shared_ptr<Dispersion> disp(new Dispersion());
        disp->name_ = "-D1";
        disp->description_ = "    Grimme's -D1 Dispersion Correction\n";
        disp->citation_ = "    Grimme, S. (2004), J. Comp. Chem., 25: 1463-1473\n";
        disp->s6_ = s6;
        disp->d_ = 23.0;
        disp->C6_ = C6_D1_;
        disp->RvdW_ = RvdW_D1_;
        disp->C6_type_ = C6_arit;
        disp->Damping_type_ = Damping_D1;
        return disp;
    } else if (boost::to_upper_copy(name) == "-D2") {
        boost::shared_ptr<Dispersion> disp(new Dispersion());
        disp->name_ = "-D2";
        disp->description_ = "    Grimme's -D2 Dispersion Correction\n";
        disp->citation_ = "    Grimme, S. (2006),  J. Comp. Chem., 27: 1787-1799\n";
        disp->s6_ = s6;
        disp->d_ = 20.0;
        disp->C6_ = C6_D2_;
        disp->RvdW_ = RvdW_D2_;
        disp->C6_type_ = C6_geom;
        disp->Damping_type_ = Damping_D1;
        return disp;
    } else if (boost::to_upper_copy(name) == "-CHG") {
        boost::shared_ptr<Dispersion> disp(new Dispersion());
        disp->name_ = "-CHG";
        disp->description_ = "    Chai and Head-Gordon Dispersion Correction\n";
        disp->citation_ = "    Chai, J.-D.; Head-Gordon, M. (2010), J. Chem. Phys., 132: 6615-6620\n";
        disp->s6_ = s6;
        disp->d_ = 6.0;
        disp->C6_ = C6_D2_;
        disp->RvdW_ = RvdW_D2_;
        disp->C6_type_ = C6_geom;
        disp->Damping_type_ = Damping_CHG;
        return disp;
    } else {
        throw PSIEXCEPTION("Dispersion: Unknown -D type specified");
    }
}
void Dispersion::print(FILE* out, int level) const
{
    if (level < 1) return;

    fprintf(out, "   => %s: Empirical Dispersion <=\n\n", name_.c_str());

    fprintf(out, "%s", description_.c_str());
    fprintf(out, "\n");

    fprintf(out, "%s", citation_.c_str());
    fprintf(out, "\n");

    fprintf(out, "    S6 = %14.6E\n", s6_);
    fprintf(out, "\n");
}
std::string Dispersion::print_energy(boost::shared_ptr<Molecule> m)
{
    double e = compute_energy(m);
    std::stringstream s;
    s.setf(ios::scientific);
    s.precision(11);

    s << "   " << name_ << " Dispersion Energy: " << e << " [H]" << endl;

    return s.str();
}
std::string Dispersion::print_gradient(boost::shared_ptr<Molecule> m)
{
    SharedMatrix G = compute_gradient(m);
    double* g = G->pointer()[0];
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
std::string Dispersion::print_hessian(boost::shared_ptr<Molecule> m)
{
    SharedMatrix H = compute_hessian(m);
    double** h = H->pointer();;

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
double Dispersion::compute_energy(boost::shared_ptr<Molecule> m)
{
    double E = 0.0;

    for (int i = 0; i < m->natom(); i++) {
        for (int j = 0; j < i; j++) {

            double C6, Rm6, f;

            double dx = m->x(j) - m->x(i);
            double dy = m->y(j) - m->y(i);
            double dz = m->z(j) - m->z(i);
    
            double R2 = dx * dx + dy * dy + dz * dz;
            double R = sqrt(R2);
            double R6 = R2 * R2 * R2;
            Rm6 = 1.0 / R6;

            double RvdW = RvdW_[(int)m->Z(i)] + RvdW_[(int)m->Z(j)];

            if (C6_type_ == C6_arit) {
                C6 = 2.0 * C6_[(int)m->Z(i)] * C6_[(int)m->Z(j)] / (C6_[(int)m->Z(i)] + C6_[(int)m->Z(j)]);
            } else if (C6_type_ == C6_geom) {
                C6 = sqrt(C6_[(int)m->Z(i)] * C6_[(int)m->Z(j)]);
            } else {
                throw PSIEXCEPTION("Unrecognized C6 Type");
            } 

            if (Damping_type_ == Damping_D1) {
                f = 1.0 / (1.0 + exp(-d_ * (R / RvdW - 1)));
            } else if (Damping_type_ == Damping_CHG) {
                f = 1.0 / (1.0 + d_ * pow((R / RvdW),-12.0));
            } else {
                throw PSIEXCEPTION("Unrecognized Damping Function");
            }

            E += C6 * Rm6 * f;
        }
    } 

    E *= - s6_;
    
    return E;
}
SharedMatrix Dispersion::compute_gradient(boost::shared_ptr<Molecule> m)
{
    SharedMatrix G(new Matrix("Dispersion Gradient", m->natom(), 3));
    double** Gp = G->pointer();

    for (int i = 0; i < m->natom(); i++) {
        for (int j = 0; j < i; j++) {

            double C6, Rm6, f;
            double C6_R, Rm6_R, f_R;
    
            double R; 
            double R_xi, R_yi, R_zi;
            double R_xj, R_yj, R_zj;

            double dx = m->x(j) - m->x(i);
            double dy = m->y(j) - m->y(i);
            double dz = m->z(j) - m->z(i);
    
            double R2 = dx * dx + dy * dy + dz * dz;
            R = sqrt(R2);

            R_xi = - dx / R;
            R_xj =   dx / R;
            R_yi = - dy / R;
            R_yj =   dy / R;
            R_zi = - dz / R;
            R_zj =   dz / R;

            double R6 = R2 * R2 * R2;
            Rm6 = 1.0 / R6;
            Rm6_R = -6.0 * Rm6 / R;

            double RvdW = RvdW_[(int)m->Z(i)] + RvdW_[(int)m->Z(j)];

            if (C6_type_ == C6_arit) {
                C6 = 2.0 * C6_[(int)m->Z(i)] * C6_[(int)m->Z(j)] / (C6_[(int)m->Z(i)] + C6_[(int)m->Z(j)]);
                C6_R = 0.0;
            } else if (C6_type_ == C6_geom) {
                C6 = sqrt(C6_[(int)m->Z(i)] * C6_[(int)m->Z(j)]);
                C6_R = 0.0;
            } else {
                throw PSIEXCEPTION("Unrecognized C6 Type");
            } 

            if (Damping_type_ == Damping_D1) {
                f = 1.0 / (1.0 + exp(-d_ * (R / RvdW - 1.0)));
                f_R = - f * f * exp(-d_ * (R / RvdW - 1.0)) * (-d_ / RvdW); 
            } else if (Damping_type_ == Damping_CHG) {
                f = 1.0 / (1.0 + d_ * pow((R / RvdW),-12.0));
                f_R = - f * f * d_ * (-12.0) * pow((R / RvdW), -13.0) * (1.0 / RvdW);
            } else {
                throw PSIEXCEPTION("Unrecognized Damping Function");
            }

            double E_R = C6_R * Rm6 * f + C6 * Rm6_R * f + C6 * Rm6 * f_R;

            Gp[i][0] += E_R * R_xi; 
            Gp[i][1] += E_R * R_yi; 
            Gp[i][2] += E_R * R_zi; 
            Gp[j][0] += E_R * R_xj; 
            Gp[j][1] += E_R * R_yj; 
            Gp[j][2] += E_R * R_zj; 
        }
    } 

    G->scale(-s6_);
    
    return G;
}
SharedMatrix Dispersion::compute_hessian(boost::shared_ptr<Molecule> m)
{
    throw PSIEXCEPTION("Dispersion: Hessians not implemented");
}

} // end namespace
