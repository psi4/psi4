#include <cstdlib>
#include <cmath>
#include <libmints/gshell.h>

#include <libmints/wavefunction.h>

#include <psiconfig.h>
#include <psi4-dec.h>

using namespace psi;

void GaussianShell::init(int ncn, int nprm, double* e, int* am, GaussianType pure,
    double** c, int nc, Vector3& center, int start, PrimitiveType pt)
{
    nprimitives_ = nprm;
    ncontractions_ = ncn;
    nc_ = nc;
    center_ = center;
    start_ = start;

    puream_ = new int[ncontraction()];
    for (int i=0; i<ncontraction(); ++i) {
        puream_[i] = (pure == Pure);
    }
 
    copy_data(am, e, c);
    // Compute the number of basis functions in this shell
    init_data();

    // Compute the normalization constants
    if (pt == Unnormalized)
        normalize_shell();
}

GaussianShell::~GaussianShell()
{
    delete[] l_;
    delete[] puream_;
    delete[] exp_;
    
    for (int i=0; i<ncontractions_; ++i)
        delete[] coef_[i];
    
    delete[] coef_;
}

// expects coef to be in primitive x contraction format. the data is transposed here
void GaussianShell::copy_data(int *l, double *exp, double **coef)
{
    l_ = new int[ncontraction()];
    for (int c=0; c<ncontraction(); ++c)
        l_[c] = l[c];
    
    exp_ = new double[nprimitive()];
    coef_ = new double*[ncontraction()];
    for (int p=0; p<nprimitive(); ++p) {
        exp_[p] = exp[p];
    }
    for (int c=0; c<ncontraction(); ++c) {
        coef_[c] = new double[nprimitive()];
        for (int p=0; p<nprimitive(); ++p) {
            // Yes, I want it stored c x p, but it came in from chkpt as p x c
            coef_[c][p] = coef[p][c];
        }
    }
}

double GaussianShell::primitive_normalization(int p)
{
    int am = l_[0];
    double tmp1 = am + 1.5;
    double g = 2.0 * exp_[p];
    double z = pow(g, tmp1);
    double normg = sqrt( (pow(2.0, am) * z) / (M_PI * sqrt(M_PI) * df[2*am]));
    return normg;
}

void GaussianShell::contraction_normalization(int gs)
{
    int i, j;
    double e_sum = 0.0, g, z;
    
    for (i=0; i<nprimitives_; ++i) {
        for (j=0; j<nprimitives_; ++j) {
            g = exp_[i] + exp_[j];
            z = pow(g, l_[0]+1.5);
            e_sum += coef_[gs][i] * coef_[gs][j] / z;
        }
    }
    
    double tmp = ((2.0*M_PI/M_2_SQRTPI) * df[2*l_[0]])/pow(2.0, l_[0]);
    double norm = sqrt(1.0 / (tmp*e_sum));
    
    // Set the normalization
    for (i=0; i<nprimitives_; ++i)
        coef_[gs][i] *= norm;
        
    if (std::isnan(norm))
        for (i=0; i<nprimitives_; ++i)
            coef_[gs][i] = 1.0;
}

void GaussianShell::normalize_shell()
{
    int i, gs;
    
    for (gs = 0; gs < ncontractions_; ++gs) {
        for (i = 0; i < nprimitives_; ++i) {
            double normalization = primitive_normalization(i);
            coef_[gs][i] *= normalization;
        }
        contraction_normalization(gs);
    }
}

int GaussianShell::nfunction(int c) const
{
    return INT_NFUNC(puream_[c], l_[c]);
}

void GaussianShell::init_data()
{
    int max = 0;
    int min = 0;
    int nc = 0;
    int nf = 0;
    has_pure_ = false;
    
    for (int i=0; i<ncontraction(); ++i) {
        int maxi = l_[i];
        if (max < maxi)
            max = maxi;
            
        int mini = l_[i];
        if (min > mini || i == 0)
            min = mini;
            
        nc += ncartesian(i);
        nf += nfunction(i);
        
        if (is_pure(i))
            has_pure_ = true;
    }
    
    max_am_ = max;
    min_am_ = min;
    ncartesians_ = nc;
    nfunctions_ = nf;
}

void GaussianShell::print(FILE *out) const
{
    fprintf(out, "      Number of contractions: %d\n", ncontraction());
    fprintf(out, "      Number of primitives: %d\n", nprimitive());
    fprintf(out, "      Number of Cartesian Gaussians: %d\n", ncartesian());
    fprintf(out, "      Spherical Harmonics?: %s\n", has_pure() ? "true" : "false");
    if (max_am() == min_am())
        fprintf(out, "      Angular momentum: %d\n", max_am());
    else {
        fprintf(out, "      Max angular momentum: %d\n", max_am());
        fprintf(out, "      Min angular momentum: %d\n", min_am());
    }
    fprintf(out, "      Center: %d\n", nc_);
    fprintf(out, "      Start index: %d\n", start_);
    fprintf(out, "      # Cartesians: %d\n", ncartesians_);
    fprintf(out, "      # functions: %d\n", nfunctions_);
    fprintf(out, "      Exponent       ");
    for(int c=0; c<ncontraction(); c++)
        fprintf(out, " Contr. %3d  ",c);
    fprintf(out, "\n");
    for(int p=0; p<nprimitive(); p++) {
        fprintf(out, "      %15.10f",exp_[p]);
        for(int c=0; c<ncontraction(); c++)
            fprintf(out, " %12.9f",coef_[c][p]);
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
        double numer = df[2*l_[0]];
        double denom = df[2*l] * df[2*m] * df[2*n];
        // printf("l_=%d, l=%d, m=%d, n=%d, norm=%15.10f\n", l_[0], l, m, n, sqrt(numer/denom));
        return sqrt(numer/denom);        
    }
}

