#include <cstdlib>
#include <cmath>
#include "vector3.h"
#include "integral.h"
#include "gshell.h"

#include <libmints/wavefunction.h>

#include <psiconfig.h>
#include <psi4-dec.h>

using namespace psi;

void GaussianShell::init(int nprm, double* e, int am, GaussianType pure,
    double* c, int nc, Vector3& center, int start, PrimitiveType pt)
{
    nprimitive_ = nprm;
    nc_ = nc;
    center_ = center;
    start_ = start;

    puream_ = pure;

    // Directly copy the arrays over
    copy_data(am, e, c);
    // Compute the number of basis functions in this shell
    init_data();

    // Compute the normalization constants
    if (pt == Unnormalized)
        normalize_shell();
}

GaussianShell::~GaussianShell()
{
    delete[] exp_;
    delete[] coef_;
}

// expects coef to be in primitive x contraction format. the data is transposed here
void GaussianShell::copy_data(int l, double *exp, double *coef)
{
    l_ = l;

    exp_ = new double[nprimitive()];
    coef_ = new double[nprimitive()];
    for (int p=0; p<nprimitive(); ++p) {
        exp_[p] = exp[p];
        coef_[p] = coef[p];
    }
}

GaussianShell* GaussianShell::copy(int nc, Vector3& c)
{
    GaussianShell* temp = new GaussianShell();
    temp->init(nprimitive_, exp_, l_, GaussianType(puream_), coef_, nc, c, start_);
    return temp;
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

    for (i=0; i<nprimitive_; ++i) {
        for (j=0; j<nprimitive_; ++j) {
            g = exp_[i] + exp_[j];
            z = pow(g, l_+1.5);
            e_sum += coef_[i] * coef_[j] / z;
        }
    }

    double tmp = ((2.0*M_PI/M_2_SQRTPI) * df[2*l_])/pow(2.0, l_);
    double norm = sqrt(1.0 / (tmp*e_sum));

    // Set the normalization
    for (i=0; i<nprimitive_; ++i)
        coef_[i] *= norm;

    //if (std::isnan(norm))
    if (norm != norm) 
        for (i=0; i<nprimitive_; ++i)
            coef_[i] = 1.0;
}

void GaussianShell::normalize_shell()
{
    int i;

    for (i = 0; i < nprimitive_; ++i) {
        double normalization = primitive_normalization(i);
        coef_[i] *= normalization;
    }
    contraction_normalization();
}

int GaussianShell::nfunction() const
{
    return INT_NFUNC(puream_, l_);
}

void GaussianShell::init_data()
{
    ncartesian_ = INT_NCART(l_);
    nfunction_ = INT_NFUNC(puream_, l_);
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
    static int use_cca_integrals_standard = (PSI_INTEGRALS_STANDARD == 1);
    if (use_cca_integrals_standard) {
        return 1.0;
    } else {
        double numer = df[2*l_];
        double denom = df[2*l] * df[2*m] * df[2*n];
        // printf("l_=%d, l=%d, m=%d, n=%d, norm=%15.10f\n", l_[0], l, m, n, sqrt(numer/denom));
        return sqrt(numer/denom);
    }
}

Vector3 GaussianShell::center() const
{
    return center_;
}
const char *GaussianShell::amtypes = "spdfghiklmn";
const char *GaussianShell::AMTYPES = "SPDFGHIKLMN";

