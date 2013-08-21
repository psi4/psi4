#include "psiconfig.h"

#ifdef HAVE_ERD

#include "psi4-dec.h"
#include "erd_eri.h"
#include "libmints/basisset.h"
#include "libmints/integral.h"
#include "libmints/wavefunction.h"

#define DEBUG 0

// Name mangling
#define FC_SYMBOL 2

#if FC_SYMBOL == 1
#define C_ERD__GENER_ERI_BATCH erd__gener_eri_batch
#define C_ERD__MEMORY_CSGTO erd__memory_csgto
#define C_ERD__MEMORY_ERI_BATCH erd__memory_eri_batch
#elif FC_SYMBOL == 2
#define C_ERD__GENER_ERI_BATCH erd__gener_eri_batch_
#define C_ERD__MEMORY_CSGTO erd__memory_csgto_
#define C_ERD__MEMORY_ERI_BATCH erd__memory_eri_batch_
#elif FC_SYMBOL == 3
#define C_ERD__GENER_ERI_BATCH ERD__GENER_ERI_BATCH
#define C_ERD__MEMORY_CSGTO ERD__MEMORY_CSGTO
#define C_ERD__MEMORY_ERI_BATCH ERD__MEMORY_ERI_BATCH
#elif FC_SYMBOL == 4
#define C_ERD__GENER_ERI_BATCH ERD__GENER_ERI_BATCH_
#define C_ERD__MEMORY_CSGTO ERD__MEMORY_CSGTO_
#define C_ERD__MEMORY_ERI_BATCH ERD__MEMORY_ERI_BATCH_
#else
#error FC_SYMBOL is not defined
#endif

typedef int F_INT;
typedef int F_BOOL;

extern "C" {
void C_ERD__GENER_ERI_BATCH(const F_INT &imax, const F_INT &zmax, const F_INT &nalpha, const F_INT &ncoeff,
                            const F_INT &ncsum, const F_INT &ncgto1, const F_INT &ncgto2,
                            const F_INT &ncgto3, const F_INT &ncgto4, const F_INT &npgto1,
                            const F_INT &npgto2, const F_INT &npgto3, const F_INT &npgto4,
                            const F_INT &shell1, const F_INT &shell2, const F_INT &shell3,
                            const F_INT &shell4, const double &x1, const double &y1, const double &z1,
                            const double &x2,const double &y2,const double &z2, const double &x3,
                            const double &y3,const double &z3,const double &x4, const double &y4, const double &z4,
                            const double *alpha, const double *cc, const F_INT *ccbeg, const F_INT *ccend,
                            const F_BOOL &spheric,  const F_BOOL &screen, F_INT *icore,
                            F_INT &nbatch, F_INT & nfirst, double *zcore );
void C_ERD__MEMORY_ERI_BATCH(const F_INT &nalpha, const F_INT &ncoeff,
                             const F_INT &ncgto1, const F_INT &ncgto2, const F_INT &ncgto3, const F_INT &ncgto4,
                             const F_INT &npgto1, const F_INT &npgto2, const F_INT &npgto3, const F_INT &npgto4,
                             const F_INT &shell1, const F_INT &shell2, const F_INT &shell3, const F_INT &shell4,
                             const double &x1, const double &y1, const double &z1, const double &x2, const double &y2,
                             const double &z2, const double &x3, const double &y3, const double &z3, const double &x4,
                             const double &y4, const double &z4, const double *alpha, const double *cc, const F_BOOL &spheric,
                             F_INT &imin, F_INT &iopt, F_INT &zmin, F_INT &zopt);
}

namespace psi{

ERDTwoElectronInt::ERDTwoElectronInt(const IntegralFactory* integral, int deriv, bool use_shell_pairs):
    TwoBodyAOInt(integral, deriv),
    d_buffer_size_(0L),
    i_buffer_size_(0L),
    buffer_offset_(0L)
{
    bs1_ = original_bs1_;
    bs2_ = original_bs2_;
    bs3_ = original_bs3_;
    bs4_ = original_bs4_;
    size_t max_cart = INT_NCART(basis1()->max_am()) * INT_NCART(basis2()->max_am()) *
                      INT_NCART(basis3()->max_am()) * INT_NCART(basis4()->max_am());

    try {
        tformbuf_ = new double[max_cart];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating tformbuf_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(tformbuf_, 0, sizeof(double)*max_cart);


    try {
        target_ = new double[max_cart];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating target_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(target_, 0, sizeof(double)*max_cart);

    size_t max_nprim = basis1()->max_nprimitive() +
                       basis2()->max_nprimitive() +
                       basis3()->max_nprimitive() +
                       basis4()->max_nprimitive();
    cc_ = new double[max_nprim];
    alpha_ = new double[max_nprim];
    new_cc_1_ = new double[original_bs1_->nprimitive()];
    new_cc_2_ = new double[original_bs2_->nprimitive()];
    new_cc_3_ = new double[original_bs3_->nprimitive()];
    new_cc_4_ = new double[original_bs4_->nprimitive()];
    cc_shell_offsets_1_ = new int[original_bs1_->nshell()];
    cc_shell_offsets_2_ = new int[original_bs2_->nshell()];
    cc_shell_offsets_3_ = new int[original_bs3_->nshell()];
    cc_shell_offsets_4_ = new int[original_bs4_->nshell()];

    // Ask ERD for the maximum amount of scratch it'll need
    compute_scratch_size();
    try {
        dscratch_ = new double[d_buffer_size_];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating dscratch_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(dscratch_, 0, sizeof(double)*d_buffer_size_);

    try {
        iscratch_ = new F_INT[i_buffer_size_];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating iscratch.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }

    normalize_basis();
}


ERDTwoElectronInt::~ERDTwoElectronInt()
{
    delete[] alpha_;
    delete[] cc_;
    delete[] cc_shell_offsets_1_;
    delete[] cc_shell_offsets_2_;
    delete[] cc_shell_offsets_3_;
    delete[] cc_shell_offsets_4_;
    delete[] new_cc_1_;
    delete[] new_cc_2_;
    delete[] new_cc_3_;
    delete[] new_cc_4_;
    delete[] tformbuf_;
    delete[] target_;
    delete[] dscratch_;
    delete[] iscratch_;
}


size_t ERDTwoElectronInt::compute_shell(const psi::AOShellCombinationsIterator& shellIter)
{
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}


void ERDTwoElectronInt::compute_scratch_size()
{
    int shell1;
    int npgto1 = 0;
    for(int shell = 0; shell < original_bs1_->nshell(); ++shell)
        if(original_bs1_->shell(shell).nprimitive() > npgto1){
            npgto1 = original_bs1_->shell(shell).nprimitive();
            shell1 = shell;
        }
    int shell2;
    int npgto2 = 0;
    for(int shell = 0; shell < original_bs2_->nshell(); ++shell)
        if(original_bs2_->shell(shell).nprimitive() > npgto2){
            npgto2 = original_bs2_->shell(shell).nprimitive();
            shell2 = shell;
        }
    int shell3;
    int npgto3 = 0;
    for(int shell = 0; shell < original_bs3_->nshell(); ++shell)
        if(original_bs3_->shell(shell).nprimitive() > npgto3){
            npgto3 = original_bs3_->shell(shell).nprimitive();
            shell3 = shell;
        }
    int shell4;
    int npgto4 = 0;
    for(int shell = 0; shell < original_bs4_->nshell(); ++shell)
        if(original_bs4_->shell(shell).nprimitive() > npgto4){
            npgto4 = original_bs4_->shell(shell).nprimitive();
            shell4 = shell;
        }
    const GaussianShell &gs1 = original_bs1_->shell(shell1);
    const GaussianShell &gs2 = original_bs2_->shell(shell2);
    const GaussianShell &gs3 = original_bs3_->shell(shell3);
    const GaussianShell &gs4 = original_bs4_->shell(shell4);
    double x1 = 1.0;
    double y1 = 1.0;
    double z1 = 1.0;
    double x2 = 2.0;
    double y2 = 2.0;
    double z2 = 2.0;
    double x3 = 3.0;
    double y3 = 3.0;
    double z3 = 3.0;
    double x4 = 4.0;
    double y4 = 4.0;
    double z4 = 4.0;
    F_INT ncgto1 = 1;
    F_INT ncgto2 = 1;
    F_INT ncgto3 = 1;
    F_INT ncgto4 = 1;
    F_INT ncgto = 4;
    F_INT npgto = npgto1 + npgto2 + npgto3 + npgto4;
    F_INT ncoeff = npgto;
    F_INT am1 = original_bs1_->max_am();
    F_INT am2 = original_bs2_->max_am();
    F_INT am3 = original_bs3_->max_am();
    F_INT am4 = original_bs4_->max_am();
    F_INT imin = 0;
    F_INT iopt = 0;
    F_INT zmin = 0;
    F_INT zopt = 0;
    F_INT last_pgto = 0;
    F_BOOL screen = 1;
    F_BOOL spheric = 0;
    for(int pgto1 = 0; pgto1 < npgto1; ++pgto1){
        cc_[last_pgto] = gs1.coef(pgto1);
        alpha_[last_pgto] = gs1.exp(pgto1);
        ++last_pgto;
    }
    for(int pgto2 = 0; pgto2 < npgto2; ++pgto2){
        cc_[last_pgto] = gs2.coef(pgto2);
        alpha_[last_pgto] = gs2.exp(pgto2);
        ++last_pgto;
    }
    for(int pgto3 = 0; pgto3 < npgto3; ++pgto3){
        cc_[last_pgto] = gs3.coef(pgto3);
        alpha_[last_pgto] = gs3.exp(pgto3);
        ++last_pgto;
    }
    for(int pgto4 = 0; pgto4 < npgto4; ++pgto4){
        cc_[last_pgto] = gs4.coef(pgto4);
        alpha_[last_pgto] = gs4.exp(pgto4);
        ++last_pgto;
    }
    long int nbatch = 0;
    // Compute the amount of memory needed for the largest quartet
    C_ERD__MEMORY_ERI_BATCH(npgto, ncoeff, ncgto1, ncgto2, ncgto3, ncgto4,
                            npgto1, npgto2, npgto3, npgto4, am1, am2, am3, am4,
                            x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, 
                            alpha_, cc_, spheric, imin, iopt, zmin, zopt);
#if DEBUG
    fprintf(outfile, "\timin %ld iopt %ld zmin %ld zopt %ld\n", imin, iopt, zmin, zopt);
#endif
    i_buffer_size_ = iopt;
    d_buffer_size_ = zopt;
}

/**
 * Applies the normalization needed for ERD, with psi4.  Note that all calls to
 * erd__normalize_cartesian must be eliminated in the ERD code for this to work.
 */
void ERDTwoElectronInt::normalize_basis()
{
    // Basis set 1
    int count = 0;
    for(int shell = 0; shell < original_bs1_->nshell(); ++shell){
        cc_shell_offsets_1_[shell] = count;
        const GaussianShell &gs = original_bs1_->shell(shell);
        int L = gs.am();
        double sum = 0.0;
        for(int j = 0; j < gs.nprimitive(); j++){
            for(int k = 0; k <= j; k++){
                double a1 = gs.exp(j);
                double a2 = gs.exp(k);
                double temp = (gs.original_coef(j) * gs.original_coef(k));
                double temp2 = ((double) L + 1.5);
                double temp3 = (2.0 * sqrt(a1 * a2) / (a1 + a2));
                temp3 = pow(temp3, temp2);
                temp = temp * temp3;
                sum = sum + temp;
                if(j != k)
                    sum = sum + temp;
            }
        }
        double prefac = 1.0;
        if(L > 1)
            prefac =pow(2.0, 2*L) / df[2*L]; 
        double norm = sqrt(prefac / sum);
        for(int j = 0; j < gs.nprimitive(); j++)
            new_cc_1_[count++] = gs.original_coef(j) * norm;
    }
    // Basis set 2
    count = 0;
    for(int shell = 0; shell < original_bs2_->nshell(); ++shell){
        cc_shell_offsets_2_[shell] = count;
        const GaussianShell &gs = original_bs2_->shell(shell);
        int L = gs.am();
        double sum = 0.0;
        for(int j = 0; j < gs.nprimitive(); j++){
            for(int k = 0; k <= j; k++){
                double a1 = gs.exp(j);
                double a2 = gs.exp(k);
                double temp = (gs.original_coef(j) * gs.original_coef(k));
                double temp2 = ((double) L + 1.5);
                double temp3 = (2.0 * sqrt(a1 * a2) / (a1 + a2));
                temp3 = pow(temp3, temp2);
                temp = temp * temp3;
                sum = sum + temp;
                if(j != k)
                    sum = sum + temp;
            }
        }
        double prefac = 1.0;
        if(L > 1)
            prefac =pow(2.0, 2*L) / df[2*L]; 
        double norm = sqrt(prefac / sum);
        for(int j = 0; j < gs.nprimitive(); j++)
            new_cc_2_[count++] = gs.original_coef(j) * norm;
    }
    // Basis set 3
    count = 0;
    for(int shell = 0; shell < original_bs3_->nshell(); ++shell){
        cc_shell_offsets_3_[shell] = count;
        const GaussianShell &gs = original_bs3_->shell(shell);
        int L = gs.am();
        double sum = 0.0;
        for(int j = 0; j < gs.nprimitive(); j++){
            for(int k = 0; k <= j; k++){
                double a1 = gs.exp(j);
                double a2 = gs.exp(k);
                double temp = (gs.original_coef(j) * gs.original_coef(k));
                double temp2 = ((double) L + 1.5);
                double temp3 = (2.0 * sqrt(a1 * a2) / (a1 + a2));
                temp3 = pow(temp3, temp2);
                temp = temp * temp3;
                sum = sum + temp;
                if(j != k)
                    sum = sum + temp;
            }
        }
        double prefac = 1.0;
        if(L > 1)
            prefac =pow(2.0, 2*L) / df[2*L]; 
        double norm = sqrt(prefac / sum);
        for(int j = 0; j < gs.nprimitive(); j++)
            new_cc_3_[count++] = gs.original_coef(j) * norm;
    }
    // Basis set 4
    count = 0;
    for(int shell = 0; shell < original_bs4_->nshell(); ++shell){
        cc_shell_offsets_4_[shell] = count;
        const GaussianShell &gs = original_bs4_->shell(shell);
        int L = gs.am();
        double sum = 0.0;
        for(int j = 0; j < gs.nprimitive(); j++){
            for(int k = 0; k <= j; k++){
                double a1 = gs.exp(j);
                double a2 = gs.exp(k);
                double temp = (gs.original_coef(j) * gs.original_coef(k));
                double temp2 = ((double) L + 1.5);
                double temp3 = (2.0 * sqrt(a1 * a2) / (a1 + a2));
                temp3 = pow(temp3, temp2);
                temp = temp * temp3;
                sum = sum + temp;
                if(j != k)
                    sum = sum + temp;
            }
        }
        double prefac = 1.0;
        if(L > 1)
            prefac =pow(2.0, 2*L) / df[2*L]; 
        double norm = sqrt(prefac / sum);
        for(int j = 0; j < gs.nprimitive(); j++)
            new_cc_4_[count++] = gs.original_coef(j) * norm;
    }
}

size_t ERDTwoElectronInt::compute_shell(int shell_i, int shell_j, int shell_k, int shell_l)
{
    const GaussianShell &gs1 = original_bs1_->shell(shell_i);
    const GaussianShell &gs2 = original_bs2_->shell(shell_j);
    const GaussianShell &gs3 = original_bs3_->shell(shell_k);
    const GaussianShell &gs4 = original_bs4_->shell(shell_l);
    Vector3 xyz1 = gs1.center();
    Vector3 xyz2 = gs2.center();
    Vector3 xyz3 = gs3.center();
    Vector3 xyz4 = gs4.center();
    double x1 = xyz1[0];
    double y1 = xyz1[1];
    double z1 = xyz1[2];
    double x2 = xyz2[0];
    double y2 = xyz2[1];
    double z2 = xyz2[2];
    double x3 = xyz3[0];
    double y3 = xyz3[1];
    double z3 = xyz3[2];
    double x4 = xyz4[0];
    double y4 = xyz4[1];
    double z4 = xyz4[2];
    F_INT ncgto1 = 1;
    F_INT ncgto2 = 1;
    F_INT ncgto3 = 1;
    F_INT ncgto4 = 1;
    F_INT ncgto = 4;
    F_INT npgto1 = gs1.nprimitive();
    F_INT npgto2 = gs2.nprimitive();
    F_INT npgto3 = gs3.nprimitive();
    F_INT npgto4 = gs4.nprimitive();
    F_INT npgto = npgto1 + npgto2 + npgto3 + npgto4;
    F_INT ncoeff = npgto;
    F_INT am1 = gs1.am();
    F_INT am2 = gs2.am();
    F_INT am3 = gs3.am();
    F_INT am4 = gs4.am();
    F_INT ccbeg[4] = {1, 1, 1, 1};
    F_INT ccend[4] = {npgto4, npgto3, npgto2, npgto1};
    F_INT nbatch;
    int offset_i = 0;
    int offset_j = offset_i + npgto4;
    int offset_k = offset_j + npgto3;
    int offset_l = offset_k + npgto2;
    F_BOOL spheric = 0;
    F_BOOL screen = 1;

    // Copy exponents and coefficients over
    ::memcpy(&(cc_[offset_i]), &(new_cc_4_[cc_shell_offsets_4_[shell_l]]), sizeof(double)*npgto4);
    ::memcpy(&(cc_[offset_j]), &(new_cc_3_[cc_shell_offsets_3_[shell_k]]), sizeof(double)*npgto3);
    ::memcpy(&(cc_[offset_k]), &(new_cc_2_[cc_shell_offsets_2_[shell_j]]), sizeof(double)*npgto2);
    ::memcpy(&(cc_[offset_l]), &(new_cc_1_[cc_shell_offsets_1_[shell_i]]), sizeof(double)*npgto1);
    int pgto = 0;
    for(int pgto4 = 0; pgto4 < npgto4; ++pgto4) alpha_[pgto++] = gs4.exp(pgto4);
    for(int pgto3 = 0; pgto3 < npgto3; ++pgto3) alpha_[pgto++] = gs3.exp(pgto3);
    for(int pgto2 = 0; pgto2 < npgto2; ++pgto2) alpha_[pgto++] = gs2.exp(pgto2);
    for(int pgto1 = 0; pgto1 < npgto1; ++pgto1) alpha_[pgto++] = gs1.exp(pgto1);
#if DEBUG
    fprintf(outfile, "\n\nShell (%2d %2d | %2d %2d) - center (%2d %2d | %2d %2d) - angular momentum (%ld %ld | %ld %ld)\n",
                    shell_i, shell_j, shell_k, shell_l,
                    gs1.ncenter(), gs2.ncenter(), gs3.ncenter(), gs4.ncenter(), 
                    am1, am2, am3, am4);
    fprintf(outfile, "Indices %ld -> %ld\n", ccbeg[0], ccend[0]);

    fprintf(outfile, "Number of primitives: %ld %ld %ld %ld\n", npgto1, npgto2, npgto3, npgto4);
    fprintf(outfile, "Coefficients: ");
    for(int n = 0; n < ncoeff; ++n)
        fprintf(outfile, "%14.10f ", cc_[n]);
    fprintf(outfile, "\n");
    fprintf(outfile, "Exponents:    ");
    for(int n = 0; n < ncoeff; ++n)
        fprintf(outfile, "%14.10f ", alpha_[n]);
    fprintf(outfile, "\n");
    fprintf(outfile, "CCBEG %3ld %3ld %3ld %3ld\n", ccbeg[0], ccbeg[1], ccbeg[2], ccbeg[3]);
    fprintf(outfile, "CCEND %3ld %3ld %3ld %3ld\n", ccend[0], ccend[1], ccend[2], ccend[3]);
#endif
    // Call ERD.  N.B. We reverse the shell ordering, because the first index is
    // the fastest running index in the buffer, which should be l for us.
    C_ERD__GENER_ERI_BATCH(i_buffer_size_, d_buffer_size_, npgto, ncoeff, ncgto,
                           ncgto4, ncgto3, ncgto2, ncgto1,
                           npgto4, npgto3, npgto2, npgto1,
                           am4, am3, am2, am1,
                           x4, y4, z4, x3, y3, z3, x2, y2, z2, x1, y1, z1,
                           alpha_, cc_, ccbeg, ccend, spheric, screen,
                           iscratch_, nbatch, buffer_offset_, dscratch_);

#if DEBUG
    fprintf(outfile, "The buffer: ");
    for(int n=0; n<d_buffer_size_;++n)
        fprintf(outfile, "%10.6f ", dscratch_[n]);
    fprintf(outfile, "\n");

    fprintf(outfile, "Buffer offset is %ld\n", buffer_offset_-1);
    fprintf(outfile, "%ld integrals were computed\n", nbatch);
#endif

    if(nbatch == 0){
        // The code should check the return value, and ignore the integrals in the buffer if we get here
        //::memset(target_, 0, sizeof(double)*gs1.nfunction() *gs2.nfunction() *gs3.nfunction() *gs4.nfunction());
        return 0;
    }

    if(original_bs1_->has_puream()){
        source_ = &(dscratch_[buffer_offset_-1]);
        pure_transform(shell_i, shell_j, shell_k, shell_l, 1);
    }else{
        ::memcpy(target_, &(dscratch_[buffer_offset_-1]), sizeof(double)*nbatch);
    }
    return nbatch;

}


size_t ERDTwoElectronInt::compute_shell_deriv1(int, int, int, int)
{
    throw PSIEXCEPTION("Derivatives for ERD are NYI!");
}


size_t ERDTwoElectronInt::compute_shell_deriv2(int, int, int, int)
{
    throw PSIEXCEPTION("Derivatives for ERD are NYI!");
}

ERDERI::ERDERI(const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : ERDTwoElectronInt(integral, deriv, use_shell_pairs)
{
}

ERDERI::~ERDERI()
{
}

} // Namespace

#endif // HAVE_ERD
