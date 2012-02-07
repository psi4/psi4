#include <cstdlib>
#include <cmath>
#include "vector3.h"
#include "integral.h"
#include "gshell.h"

#include <libmints/wavefunction.h>

#include <psiconfig.h>
#include <psi4-dec.h>

using namespace psi;

GaussianShell::GaussianShell(int am, const std::vector<double> &c,
                             const std::vector<double> &e, GaussianType pure,
                             int nc, const Vector3 &center, int start,
                             PrimitiveType pt)
    : l_(am), puream_(pure), exp_(e), coef_(c),
      nc_(nc), center_(center), start_(start)
{
    ncartesian_ = INT_NCART(l_);
    nfunction_  = INT_NFUNC(puream_, l_);

    // Compute the normalization constants
    if (pt == Unnormalized)
        normalize_shell();
}

GaussianShell GaussianShell::copy()
{
    return GaussianShell(l_, coef_, exp_,
                         GaussianType(puream_),
                         nc_, center_, start_);
}

GaussianShell GaussianShell::copy(int nc, const Vector3& c)
{
    return GaussianShell(l_, coef_, exp_,
                         GaussianType(puream_),
                         nc, c, start_);
}

double GaussianShell::primitive_normalization(int p)
{
    double tmp1 = l_ + 1.5;
    double g = 2.0 * exp_[p];
    double z = pow(g, tmp1);
    double normg = sqrt( (pow(2.0, l_) * z) / (M_PI * sqrt(M_PI) * df[2*l_]));
    return normg;
}

void GaussianShell::contraction_normalization()
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

void GaussianShell::normalize_shell()
{
    int i;

    for (i = 0; i < nprimitive(); ++i) {
        double normalization = primitive_normalization(i);
        coef_[i] *= normalization;
    }
    contraction_normalization();
}

int GaussianShell::nfunction() const
{
    return INT_NFUNC(puream_, l_);
}

int GaussianShell::nprimitive() const
{
    return exp_.size();
}

void GaussianShell::print(FILE *out) const
{
    fprintf(out, "      Number of primitives: %d\n", nprimitive());
    fprintf(out, "      Number of Cartesian Gaussians: %d\n", ncartesian());
    fprintf(out, "      Spherical Harmonics?: %s\n", is_pure() ? "true" : "false");
    fprintf(out, "      Angular momentum: %d\n", am());
    fprintf(out, "      Center: %d\n", nc_);
    fprintf(out, "      Start index: %d\n", start_);
    fprintf(out, "      # Cartesians: %d\n", ncartesian_);
    fprintf(out, "      # functions: %d\n", nfunction_);
    fprintf(out, "      Exponent       ");
    fprintf(out, "\n");
    for(int p=0; p<nprimitive(); p++) {
        fprintf(out, "      %15.10f",exp_[p]);
        fprintf(out, " %12.9f",coef_[p]);
        fprintf(out, "\n");
    }
    fprintf(out, "\n");
}

double GaussianShell::normalize(int l, int m, int n)
{
    return 1.0;
}

const Vector3& GaussianShell::center() const
{
    return center_;
}

const char *GaussianShell::amtypes = "spdfghiklmn";
const char *GaussianShell::AMTYPES = "SPDFGHIKLMN";
