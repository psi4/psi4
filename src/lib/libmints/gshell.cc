/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <cstdlib>
#include <cmath>
#include "vector3.h"
#include "integral.h"
#include "gshell.h"

#include <libmints/wavefunction.h>

#include <psiconfig.h>
#include <psi4-dec.h>

using namespace psi;

ShellInfo::ShellInfo(int am, const std::vector<double> &c,
                             const std::vector<double> &e, GaussianType pure,
                             int nc, const Vector3 &center, int start,
                             PrimitiveType pt)
    : l_(am), puream_(pure), exp_(e), coef_(c),
      nc_(nc), center_(center), start_(start)
{
    for(int n = 0; n < c.size(); ++n)
        original_coef_.push_back(c[n]);

    ncartesian_ = INT_NCART(l_);
    nfunction_  = INT_NFUNC(puream_, l_);

    // Compute the normalization constants
    if (pt == Unnormalized){
        normalize_shell();
        erd_normalize_shell();
    }
}

ShellInfo ShellInfo::copy()
{
    return ShellInfo(l_, original_coef_, exp_,
                         GaussianType(puream_),
                         nc_, center_, start_, Unnormalized);
}

ShellInfo ShellInfo::copy(int nc, const Vector3& c)
{
    return ShellInfo(l_, original_coef_, exp_,
                         GaussianType(puream_),
                         nc, c, start_, Unnormalized);
}

double ShellInfo::primitive_normalization(int p)
{
    double tmp1 = l_ + 1.5;
    double g = 2.0 * exp_[p];
    double z = pow(g, tmp1);
    double normg = sqrt( (pow(2.0, l_) * z) / (M_PI * sqrt(M_PI) * df[2*l_]));
    return normg;
}

void ShellInfo::erd_normalize_shell()
{
    erd_coef_.clear();
    double sum = 0.0;
    for(int j = 0; j < nprimitive(); j++){
        for(int k = 0; k <= j; k++){
            double a1 = exp_[j];
            double a2 = exp_[k];
            double temp = (original_coef(j) * original_coef(k));
            double temp2 = ((double) l_ + 1.5);
            double temp3 = (2.0 * sqrt(a1 * a2) / (a1 + a2));
            temp3 = pow(temp3, temp2);
            temp = temp * temp3;
            sum = sum + temp;
            if(j != k)
                sum = sum + temp;
        }
    }
    double prefac = 1.0;
    if(l_ > 1)
        prefac = pow(2.0, 2*l_) / df[2*l_];
    double norm = sqrt(prefac / sum);
    for(int j = 0; j < nprimitive(); j++){
        erd_coef_.push_back(original_coef_[j] * norm);
    }
}

void ShellInfo::contraction_normalization()
{
    int i, j;
    double e_sum = 0.0, g, z;

    for (i=0; i<nprimitive(); ++i) {
        for (j=0; j<nprimitive(); ++j) {
            g = exp_[i] + exp_[j];
            z = pow(g, l_+1.5);
            e_sum += coef_[i] * coef_[j] / z;
        }
    }

    double tmp = ((2.0*M_PI/M_2_SQRTPI) * df[2*l_])/pow(2.0, l_);
    double norm = sqrt(1.0 / (tmp*e_sum));

    // Set the normalization
    for (i=0; i<nprimitive(); ++i)
        coef_[i] *= norm;

    if (norm != norm)
        for (i=0; i<nprimitive(); ++i)
            coef_[i] = 1.0;
}

void ShellInfo::normalize_shell()
{
    int i;

    for (i = 0; i < nprimitive(); ++i) {
        double normalization = primitive_normalization(i);
        coef_[i] *= normalization;
    }
    contraction_normalization();
}

int ShellInfo::nfunction() const
{
    return INT_NFUNC(puream_, l_);
}

int ShellInfo::nprimitive() const
{
    return exp_.size();
}

void ShellInfo::print(std::string out) const
{
   boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         boost::shared_ptr<OutFile>(new OutFile(out)));

   printer->Printf( "    %c %3d 1.00\n", AMCHAR(), nprimitive());

    for (int K = 0; K < nprimitive(); K++) {

            printer->Printf( "               %20.8f %20.8f\n",exp_[K], original_coef_[K]);

    }
}

double ShellInfo::normalize(int l, int m, int n)
{
    return 1.0;
}

const Vector3& ShellInfo::center() const
{
    return center_;
}

const char *ShellInfo::amtypes = "spdfghiklmnopqrtuvwxyz";
const char *ShellInfo::AMTYPES = "SPDFGHIKLMNOPQRTUVWXYZ";


GaussianShell::GaussianShell(int am, int nprimitive, const double *oc, const double *c, const double *ec,
                             const double *e, GaussianType pure,
                             int nc, const double *center, int start)
    : l_(am), nprimitive_(nprimitive), puream_(pure), exp_(e), original_coef_(oc), coef_(c), erd_coef_(ec),
      nc_(nc), center_(center), start_(start)
{
    ncartesian_ = INT_NCART(l_);
    nfunction_  = INT_NFUNC(puream_, l_);
}

int GaussianShell::nfunction() const
{
    return INT_NFUNC(puream_, l_);
}

int GaussianShell::nprimitive() const
{
    return nprimitive_;
}

void GaussianShell::print(std::string out) const
{
   boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         boost::shared_ptr<OutFile>(new OutFile(out)));

        printer->Printf( "    %c %3d 1.00\n", AMCHAR(), nprimitive());

    for (int K = 0; K < nprimitive(); K++) {

            printer->Printf( "               %20.8f %20.8f\n",exp_[K], original_coef_[K]);

    }
}


const double* GaussianShell::center() const
{
    return center_;
}

const char *GaussianShell::amtypes = "spdfghiklmnopqrtuvwxyz";
const char *GaussianShell::AMTYPES = "SPDFGHIKLMNOPQRTUVWXYZ";
