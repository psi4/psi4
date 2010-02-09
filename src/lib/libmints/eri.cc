#include <stdexcept>
#include <string>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

#include <libmints/basisset.h>
#include <libmints/gshell.h>
#include <libmints/overlap.h>
#include <libmints/eri.h>
#include <libmints/wavefunction.h>
#include <physconst.h>
#include <exception.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define EPS 1.0e-17
#define TABLESIZE 121

using namespace psi;

inline void calc_f(double *F, int n, double t)
{
    int i, m, k;
    int m2;
    double t2;
    double num;
    double sum;
    double term1, term2;
    static double K = 1.0/M_2_SQRTPI;
    double et;

    if (t>20.0){
        t2 = 2*t;
        et = exp(-t);
        t = sqrt(t);
        F[0] = K*erf(t)/t;
        for(m=0; m<=n-1; m++){
            F[m+1] = ((2*m + 1)*F[m] - et)/(t2);
        }
    }
    else {
        et = exp(-t);
        t2 = 2*t;
        m2 = 2*n;
        num = df[m2];
        i=0;
        sum = 1.0/(m2+1);
        do{
            i++;
            num = num*t2;
            term1 = num/df[m2+2*i+2];
            sum += term1;
        } while (fabs(term1) > EPS && i < MAX_FAC);
        F[n] = sum*et;
        for(m=n-1;m>=0;m--){
            F[m] = (t2*F[m+1] + et)/(2*m+1);
        }
    }
}

ERI::ERI(shared_ptr<BasisSet> bs1, shared_ptr<BasisSet>bs2, shared_ptr<BasisSet>bs3, shared_ptr<BasisSet>bs4, int deriv, double schwarz)
    : TwoBodyInt(bs1, bs2, bs3, bs4, deriv)
{
    // Initialize libint static data
    init_libint_base();
    if (deriv_)
        init_libderiv_base();

    // Figure out some information to initialize libint/libderiv with
    // 1. Maximum angular momentum
    int max_am = MAX(MAX(bs1->max_am(), bs2->max_am()), MAX(bs3->max_am(), bs4->max_am()));
    // 2. Maximum number of primitive combinations
    int max_nprim = bs1->max_nprimitive() * bs2->max_nprimitive() * bs3->max_nprimitive() * bs4->max_nprimitive();
    // 3. Maximum Cartesian class size
    max_cart_ = ioff[bs1->max_am()] * ioff[bs2->max_am()] * ioff[bs3->max_am()] * ioff[bs4->max_am()] +1;

    // Make sure libint is compiled to handle our max AM
    if (max_am >= LIBINT_MAX_AM) {
        fprintf(stderr, "ERROR: ERI - libint cannot handle angular momentum this high.\n" \
                        "       Recompile libint for higher angular momentum, then recompile this program.\n");
        throw LimitExceeded<int>("ERI - libint cannot handle angular momentum this high.\nRecompile libint for higher angular momentum, then recompile this program.", LIBINT_MAX_AM, max_am, __FILE__, __LINE__);
        //abort();
    }
    if (deriv_ == 1 && max_am >= LIBDERIV_MAX_AM1) {
        fprintf(stderr, "ERROR: ERI - libderiv cannot handle angular momentum this high.\n" \
                        "     Recompile libderiv for higher angular momentum, then recompile this program.\n");
        throw LimitExceeded<int>("ERI - libderiv cannot handle angular momentum this high.\n"
                                 "Recompile libint for higher angular momentum, then recompile this program.",
                                 LIBDERIV_MAX_AM1, max_am, __FILE__, __LINE__);
        //abort();
    }

    try {
        // Initialize libint
        init_libint(&libint_, max_am, max_nprim);
        // and libderiv, if needed
        if (deriv_)
            init_libderiv1(&libderiv_, max_am, max_nprim, max_cart_-1);
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating memory for libint/libderiv.\n");
        exit(EXIT_FAILURE);
    }
    size_t size = INT_NCART(bs1->max_am()) * INT_NCART(bs2->max_am()) *
                  INT_NCART(bs3->max_am()) * INT_NCART(bs4->max_am());

    // Used in pure_transform
    try {
        tformbuf_ = new double[size];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating tformbuf_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(tformbuf_, 0, sizeof(double)*size);

    if (deriv_ == 1)
        size *= 3*natom_;

    try {
        target_ = new double[size];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating target_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(target_, 0, sizeof(double)*size);

    try {
        source_ = new double[size];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating source_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(source_, 0, sizeof(double)*size);

    init_fjt(4*max_am + DERIV_LVL);

    screen_ = false;
    if (schwarz != 0.0)
        schwarz2_ = schwarz*schwarz;
    else
        schwarz2_ = 0.0;

    // Precompute a bunch of information
    init_shell_pairs12();
    // If basis3 and basis4 equals basis1 and basis2, then the following function will do nothing,
    // except assign pairs34_ to pairs12_
    init_shell_pairs34();
}

ERI::~ERI()
{
    delete[] tformbuf_;
    delete[] target_;
    delete[] source_;
    delete[] denom_;
    free_block(d_);
    free_libint(&libint_);
    if (deriv_)
        free_libderiv(&libderiv_);
    if (screen_)
        free(schwarz_norm_);
    free_shell_pairs12();
    free_shell_pairs34();
}

void ERI::init_shell_pairs12()
{
    ShellPair **pairs, *sp;
    Vector3 P, PA, PB, AB, A, B;
    int i, j, si, sj, np_i, np_j;
    size_t memd;
    double a1, a2, ab2, gam, c1, c2;
    double *curr_stack_ptr;

    fprintf(outfile, "  Pre-computing some values for two-electron integrals.\n");

    // Estimate memory needed by allocated space for the dynamically allocated parts of ShellPair structure
    memd = ERI::memory_to_store_shell_pairs(basis1(), basis2());

    // Allocate a stack of memory
    stack12_ = new double[memd];
    curr_stack_ptr = stack12_;

    // Allocate shell pair memory
    pairs12_ = new ShellPair*[basis1()->nshell()];
    for (i=0; i<basis1()->nshell(); ++i)
        pairs12_[i] = new ShellPair[basis2()->nshell()];

    // Loop over all shell pairs (si, sj) and create primitive pairs pairs
    for (si=0; si<basis1()->nshell(); ++si) {
        A = basis1()->shell(si)->center();

        for (sj=0; sj<basis2()->nshell(); ++sj) {
            B = basis2()->shell(sj)->center();

            AB = A - B;
            ab2 = AB.dot(AB);

            // Get the pointer for convenience
            sp = &(pairs12_[si][sj]);

            // Save some information
            sp->i = si;
            sp->j = sj;
            sp->AB[0] = AB[0]; sp->AB[1] = AB[1]; sp->AB[2] = AB[2];

            np_i = basis1()->shell(si)->nprimitive();
            np_j = basis2()->shell(sj)->nprimitive();

            // Reserve some memory for the primitives
            sp->ai = curr_stack_ptr; curr_stack_ptr += np_i;
            sp->aj = curr_stack_ptr; curr_stack_ptr += np_j;

            // Allocate and reserve memory for gammas
            sp->gamma = new double*[np_i];
            for (i=0; i<np_i; ++i) {
                sp->gamma[i] = curr_stack_ptr; curr_stack_ptr += np_j;
            }

            // Reserve space for contraction coefficients
            sp->ci = curr_stack_ptr; curr_stack_ptr += np_i;
            sp->cj = curr_stack_ptr; curr_stack_ptr += np_j;

            // Allocate and reserve space for overlaps
            sp->overlap = new double*[np_i];
            for (i=0; i<np_i; ++i) {
                sp->overlap[i] = curr_stack_ptr; curr_stack_ptr += np_j;
            }

            // Allocate and reserve space for P, PA, and PB.
            sp->P  = new double**[np_i];
            sp->PA = new double**[np_i];
            sp->PB = new double**[np_i];
            for (i=0; i<np_i; ++i) {
                sp->P[i]  = new double*[np_j];
                sp->PA[i] = new double*[np_j];
                sp->PB[i] = new double*[np_j];

                for (j=0; j<np_j; ++j) {
                    sp->P[i][j]  = curr_stack_ptr; curr_stack_ptr += 3;
                    sp->PA[i][j] = curr_stack_ptr; curr_stack_ptr += 3;
                    sp->PB[i][j] = curr_stack_ptr; curr_stack_ptr += 3;
                }
            }

            // All memory has been reserved/allocated for this shell primitive pair pair.
            // Pre-compute all data that we can:
            for (i=0; i<np_i; ++i) {
                a1 = basis1()->shell(si)->exp(i);
                c1 = basis1()->shell(si)->coef(0, i);

                // Save some information
                sp->ai[i] = a1;
                sp->ci[i] = c1;

                for (j=0; j<np_j; ++j) {
                    a2 = basis2()->shell(sj)->exp(j);
                    c2 = basis2()->shell(sj)->coef(0, j);

                    gam = a1 + a2;

                    // Compute Gaussian product and component distances
                    P = ( A * a1 + B * a2 ) / gam;
                    PA = P - A;
                    PB = P - B;

                    // Copy data into pairs array
                    sp->aj[j] = a2;
                    sp->cj[j] = c2;
                    sp->gamma[i][j] = gam;
                    sp->P[i][j][0]  = P[0];  sp->P[i][j][1]  = P[1];  sp->P[i][j][2]  = P[2];
                    sp->PA[i][j][0] = PA[0]; sp->PA[i][j][1] = PA[1]; sp->PA[i][j][2] = PA[2];
                    sp->PB[i][j][0] = PB[0]; sp->PB[i][j][1] = PB[1]; sp->PB[i][j][2] = PB[2];
                    sp->overlap[i][j] = pow(M_PI/gam, 3.0/2.0) * exp(-a1*a2*ab2/gam) * c1 * c2;
                }
            }
        }
    }
}

void ERI::init_shell_pairs34()
{
    ShellPair **pairs, *sp;
    Vector3 P, PA, PB, AB, A, B;
    int i, j, si, sj, np_i, np_j;
    size_t memd;
    double a1, a2, ab2, gam, c1, c2;
    double *curr_stack_ptr;

    // If basis1 == basis3 && basis2 == basis4, then we don't need to do anything except use the pointer
    // of pairs12_.
    if (basis1() == basis3() && basis2() == basis4()) {
        pairs34_ = pairs12_;
        stack34_ = NULL;
        return;
    }

    fprintf(outfile, "  Pre-computing additional values for two-electron integrals. [ |34) does not equal (12| ]\n");

    // Estimate memory needed by allocated space for the dynamically allocated parts of ShellPair structure
    memd = ERI::memory_to_store_shell_pairs(basis3(), basis4());

    // Allocate a stack of memory
    stack34_ = new double[memd];
    curr_stack_ptr = stack34_;

    // Allocate shell pair memory
    pairs34_ = new ShellPair*[basis3()->nshell()];
    for (i=0; i<basis3()->nshell(); ++i)
        pairs34_[i] = new ShellPair[basis4()->nshell()];

    // Loop over all shell pairs (si, sj) and create primitive pairs pairs
    for (si=0; si<basis3()->nshell(); ++si) {
        A = basis3()->shell(si)->center();

        for (sj=0; sj<basis4()->nshell(); ++sj) {
            B = basis4()->shell(sj)->center();

            AB = A - B;
            ab2 = AB.dot(AB);

            // Get the pointer for convenience
            sp = &(pairs34_[si][sj]);

            // Save some information
            sp->i = si;
            sp->j = sj;
            sp->AB[0] = AB[0]; sp->AB[1] = AB[1]; sp->AB[2] = AB[2];

            np_i = basis3()->shell(si)->nprimitive();
            np_j = basis4()->shell(sj)->nprimitive();

            // Reserve some memory for the primitives
            sp->ai = curr_stack_ptr; curr_stack_ptr += np_i;
            sp->aj = curr_stack_ptr; curr_stack_ptr += np_j;

            // Allocate and reserve memory for gammas
            sp->gamma = new double*[np_i];
            for (i=0; i<np_i; ++i) {
                sp->gamma[i] = curr_stack_ptr; curr_stack_ptr += np_j;
            }

            // Reserve space for contraction coefficients
            sp->ci = curr_stack_ptr; curr_stack_ptr += np_i;
            sp->cj = curr_stack_ptr; curr_stack_ptr += np_j;

            // Allocate and reserve space for overlaps
            sp->overlap = new double*[np_i];
            for (i=0; i<np_i; ++i) {
                sp->overlap[i] = curr_stack_ptr; curr_stack_ptr += np_j;
            }

            // Allocate and reserve space for P, PA, and PB.
            sp->P  = new double**[np_i];
            sp->PA = new double**[np_i];
            sp->PB = new double**[np_i];
            for (i=0; i<np_i; ++i) {
                sp->P[i]  = new double*[np_j];
                sp->PA[i] = new double*[np_j];
                sp->PB[i] = new double*[np_j];

                for (j=0; j<np_j; ++j) {
                    sp->P[i][j]  = curr_stack_ptr; curr_stack_ptr += 3;
                    sp->PA[i][j] = curr_stack_ptr; curr_stack_ptr += 3;
                    sp->PB[i][j] = curr_stack_ptr; curr_stack_ptr += 3;
                }
            }

            // All memory has been reserved/allocated for this shell primitive pair pair.
            // Pre-compute all data that we can:
            for (i=0; i<np_i; ++i) {
                a1 = basis3()->shell(si)->exp(i);
                c1 = basis3()->shell(si)->coef(0, i);

                // Save some information
                sp->ai[i] = a1;
                sp->ci[i] = c1;

                for (j=0; j<np_j; ++j) {
                    a2 = basis4()->shell(sj)->exp(j);
                    c2 = basis4()->shell(sj)->coef(0, j);

                    gam = a1 + a2;

                    // Compute some distances
                    P = ( A * a1 + B * a2 ) / gam;
                    PA = P - A;
                    PB = P - B;

                    // Copy data into pairs array
                    sp->aj[j] = a2;
                    sp->cj[j] = c2;
                    sp->gamma[i][j] = gam;
                    sp->P[i][j][0]  = P[0];  sp->P[i][j][1]  = P[1];  sp->P[i][j][2]  = P[2];
                    sp->PA[i][j][0] = PA[0]; sp->PA[i][j][1] = PA[1]; sp->PA[i][j][2] = PA[2];
                    sp->PB[i][j][0] = PB[0]; sp->PB[i][j][1] = PB[1]; sp->PB[i][j][2] = PB[2];
                    sp->overlap[i][j] = pow(M_PI/gam, 3.0/2.0) * exp(-a1*a2*ab2/gam) * c1 * c2;
                }
            }
        }
    }
}

void ERI::free_shell_pairs12()
{
    int i, si, sj;
    ShellPair *sp;
    int np_i;

    free(stack12_);
    for (si=0; si<basis1()->nshell(); ++si) {
        for (sj=0; sj<basis2()->nshell(); ++sj) {
            np_i = basis1()->shell(si)->nprimitive();
            sp = &(pairs12_[si][sj]);

            delete[] sp->gamma;
            delete[] sp->overlap;

            if (sp->P != NULL) {
                for (i=0; i<np_i; ++i)
                    delete[] sp->P[i];
                delete[] sp->P;
            }
            if (sp->PA != NULL) {
                for (i=0; i<np_i; ++i)
                    delete[] sp->PA[i];
                delete[] sp->PA;
            }
            if (sp->PB != NULL) {
                for (i=0; i<np_i; ++i)
                    delete[] sp->PB[i];
                delete[] sp->PB;
            }
        }
    }

    for (si=0; si<basis1()->nshell(); ++si)
        delete[] pairs12_[si];
    delete[] pairs12_;
}

void ERI::free_shell_pairs34()
{
    int i, si, sj;
    ShellPair *sp;
    int np_i;

    // If stack34_ is NULL then we only used pairs12.
    if (stack34_ == NULL)
        return;

    free(stack34_);
    for (si=0; si<basis3()->nshell(); ++si) {
        for (sj=0; sj<basis4()->nshell(); ++sj) {
            np_i = basis3()->shell(si)->nprimitive();
            sp = &(pairs34_[si][sj]);

            delete[] sp->gamma;
            delete[] sp->overlap;

            if (sp->P != NULL) {
                for (i=0; i<np_i; ++i)
                    delete[] sp->P[i];
                delete[] sp->P;
            }
            if (sp->PA != NULL) {
                for (i=0; i<np_i; ++i)
                    delete[] sp->PA[i];
                delete[] sp->PA;
            }
            if (sp->PB != NULL) {
                for (i=0; i<np_i; ++i)
                    delete[] sp->PB[i];
                delete[] sp->PB;
            }
        }
    }

    for (si=0; si<basis3()->nshell(); ++si)
        delete[] pairs34_[si];
    delete[] pairs34_;
}

size_t ERI::memory_to_store_shell_pairs(const shared_ptr<BasisSet> &bs1, const shared_ptr<BasisSet> &bs2)
{
    int i, j, np_i, np_j;
    size_t mem = 0;

    for (i=0; i<bs1->nshell(); ++i) {
        np_i = bs1->shell(i)->nprimitive();
        for (j=0; j<bs2->nshell(); ++j) {
            np_j = bs2->shell(j)->nprimitive();
            mem += (2*(np_i + np_j) + 11*np_i*np_j);
        }
    }
    return mem;
}

// Taken from CINTS
void ERI::init_fjt(int max)
{
    int i,j;
    double denom,d2jmax1,r2jmax1,wval,d2wval,sum,term,rexpw;

    int n1 = max+7;
    int n2 = TABLESIZE;
    d_ = block_matrix(n1, n2);

    /* Tabulate the gamma function for t(=wval)=0.0. */
    denom = 1.0;
    for (i=0; i<n1; i++) {
        d_[i][0] = 1.0/denom;
        denom += 2.0;
    }

    /* Tabulate the gamma function from t(=wval)=0.1, to 12.0. */
    d2jmax1 = 2.0*(n1-1) + 1.0;
    r2jmax1 = 1.0/d2jmax1;
    for (i=1; i<TABLESIZE; i++) {
        wval = 0.1 * i;
        d2wval = 2.0 * wval;
        term = r2jmax1;
        sum = term;
        denom = d2jmax1;
        for (j=2; j<=200; j++) {
            denom = denom + 2.0;
            term = term * d2wval / denom;
            sum = sum + term;
            if (term <= 1.0e-15) break;
        }
        rexpw = exp(-wval);

      /* Fill in the values for the highest j gtable entries (top row). */
        d_[n1-1][i] = rexpw * sum;

      /* Work down the table filling in the rest of the column. */
        denom = d2jmax1;
        for (j=n1 - 2; j>=0; j--) {
            denom = denom - 2.0;
            d_[j][i] = (d_[j+1][i]*d2wval + rexpw)/denom;
        }
    }

    /* Form some denominators, so divisions can be eliminated below. */
    denom_ = new double[max+1];
    denom_[0] = 0.0;
    for (i=1; i<=max; i++) {
        denom_[i] = 1.0/(2*i - 1);
    }

    wval_infinity_ = 2*max + 37.0;
    itable_infinity_ = (int) (10 * wval_infinity_);
}

void ERI::int_fjt(double *F, int J, double wval)
{
    const double sqrpih =  0.886226925452758;
    const double coef2 =  0.5000000000000000;
    const double coef3 = -0.1666666666666667;
    const double coef4 =  0.0416666666666667;
    const double coef5 = -0.0083333333333333;
    const double coef6 =  0.0013888888888889;
    const double gfac30 =  0.4999489092;
    const double gfac31 = -0.2473631686;
    const double gfac32 =  0.321180909;
    const double gfac33 = -0.3811559346;
    const double gfac20 =  0.4998436875;
    const double gfac21 = -0.24249438;
    const double gfac22 =  0.24642845;
    const double gfac10 =  0.499093162;
    const double gfac11 = -0.2152832;
    const double gfac00 = -0.490;

    double wdif, d2wal, rexpw, /* denom, */ gval, factor, rwval, term;
    int i, itable, irange;

    /* Compute an index into the table. */
    /* The test is needed to avoid floating point exceptions for
    * large values of wval. */
    if (wval > wval_infinity_) {
        itable = itable_infinity_;
    }
    else {
        itable = (int) (10.0 * wval);
    }

    /* If itable is small enough use the table to compute int_fjttable. */
    if (itable < TABLESIZE) {

        wdif = wval - 0.1 * itable;

      /* Compute fjt for J. */
        F[J] = (((((coef6 * d_[J+6][itable]*wdif
            + coef5 * d_[J+5][itable])*wdif
            + coef4 * d_[J+4][itable])*wdif
            + coef3 * d_[J+3][itable])*wdif
            + coef2 * d_[J+2][itable])*wdif
            -  d_[J+1][itable])*wdif
            +  d_[J][itable];

      /* Compute the rest of the fjt. */
        d2wal = 2.0 * wval;
        rexpw = exp(-wval);
      /* denom = 2*J + 1; */
        for (i=J; i>0; i--) {
        /* denom = denom - 2.0; */
            F[i-1] = (d2wal*F[i] + rexpw)*denom_[i];
        }
    }
    /* If wval <= 2*J + 36.0, use the following formula. */
    else if (itable <= 20*J + 360) {
        rwval = 1.0/wval;
        rexpw = exp(-wval);

      /* Subdivide wval into 6 ranges. */
        irange = itable/30 - 3;
        if (irange == 1) {
            gval = gfac30 + rwval*(gfac31 + rwval*(gfac32 + rwval*gfac33));
            F[0] = sqrpih*sqrt(rwval) - rexpw*gval*rwval;
        }
        else if (irange == 2) {
            gval = gfac20 + rwval*(gfac21 + rwval*gfac22);
            F[0] = sqrpih*sqrt(rwval) - rexpw*gval*rwval;
        }
        else if (irange == 3 || irange == 4) {
            gval = gfac10 + rwval*gfac11;
            F[0] = sqrpih*sqrt(rwval) - rexpw*gval*rwval;
        }
        else if (irange == 5 || irange == 6) {
            gval = gfac00;
            F[0] = sqrpih*sqrt(rwval) - rexpw*gval*rwval;
        }
        else {
            F[0] = sqrpih*sqrt(rwval);
        }

      /* Compute the rest of the int_fjttable from table->d[0]. */
        factor = 0.5 * rwval;
        term = factor * rexpw;
        for (i=1; i<=J; i++) {
            F[i] = factor * F[i-1] - term;
            factor = rwval + factor;
        }
    }
    /* For large values of wval use this algorithm: */
    else {
        rwval = 1.0/wval;
        F[0] = sqrpih*sqrt(rwval);
        factor = 0.5 * rwval;
        for (i=1; i<=J; i++) {
            F[i] = factor * F[i-1];
            factor = rwval + factor;
        }
    }
}

void ERI::compute_shell(int sh1, int sh2, int sh3, int sh4)
{
    //SCHWARZ SIEVE
    //Determines if shell needs to be computed.
    //If false, places zeros in target_
    //If true, standard algorithm continues
    //
    //NOTE: Schwarz Sieve only holds if all indices correspond to the same basis.
    //A calling code may only access the schwarz screening functionality by literally
    //specifying a cutoff index. This should only be done for ERIs all based on the same basis set
    if (schwarz2_ != 0.0)
    {
#ifdef MINTS_TIMER
    timer_on("sieving");
#endif

        if (screen_ == false)
            form_sieve();

        if (schwarz_norm_[ioff[((sh1>sh2)?sh1:sh2)]+((sh1>sh2)?sh2:sh1)]*schwarz_norm_[ioff[((sh3>sh4)?sh3:sh4)]+((sh3>sh4)?sh4:sh3)]<schwarz2_)
        {
            size_t nfill = bs1_->shell(sh1)->nfunction();
            nfill*=     bs2_->shell(sh2)->nfunction();
            nfill*=     bs3_->shell(sh3)->nfunction();
            nfill*=     bs4_->shell(sh4)->nfunction();
            memset(target_, 0, sizeof(double)*nfill);
            return;
        }
#ifdef MINTS_TIMER
    timer_off("sieving");
#endif
    }
    //END OF SCHWARZ SIEVE

#ifdef MINTS_TIMER
    timer_on("ERI::compute_shell");
#endif
    // Need to ensure the ordering asked by the user is valid for libint
    // compute_quartet does NOT check this. SEGFAULTS should occur if order
    // is not guaranteed.
#ifdef MINTS_TIMER
    timer_on("reorder");
#endif

    int s1, s2, s3, s4;
    int am1, am2, am3, am4, temp;
    shared_ptr<BasisSet> bs_temp;

    bool p13p24 = false, p12 = false, p34 = false;

    // AM used for ordering
    am1 = original_bs1_->shell(sh1)->am(0);
    am2 = original_bs2_->shell(sh2)->am(0);
    am3 = original_bs3_->shell(sh3)->am(0);
    am4 = original_bs4_->shell(sh4)->am(0);

    int n1 = original_bs1_->shell(sh1)->nfunction(0);
    int n2 = original_bs2_->shell(sh2)->nfunction(0);
    int n3 = original_bs3_->shell(sh3)->nfunction(0);
    int n4 = original_bs4_->shell(sh4)->nfunction(0);

    // l(a) >= l(b), l(c) >= l(d), and l(c) + l(d) >= l(a) + l(b).
    if (am1 >= am2) {
        s1 = sh1;
        s2 = sh2;

        bs1_ = original_bs1_;
        bs2_ = original_bs2_;
    } else {
        s1 = sh2;
        s2 = sh1;

        bs1_ = original_bs2_;
        bs2_ = original_bs1_;

        p12 = true;
    }

    if (am3 >= am4) {
        s3 = sh3;
        s4 = sh4;

        bs3_ = original_bs3_;
        bs4_ = original_bs4_;

    } else {
        s3 = sh4;
        s4 = sh3;

        bs3_ = original_bs4_;
        bs4_ = original_bs3_;

        p34 = true;
    }

    if ((am1 + am2) > (am3 + am4)) {
        // Swap s1 and s2 with s3 and s4
        temp = s1;
        s1 = s3;
        s3 = temp;

        temp = s2;
        s2 = s4;
        s4 = temp;

        bs_temp = bs1_;
        bs1_ = bs3_;
        bs3_ = bs_temp;

        bs_temp = bs2_;
        bs2_ = bs4_;
        bs4_ = bs_temp;

        p13p24 = true;
    }
#ifdef MINTS_TIMER
    timer_off("reorder");
#endif

    // s1, s2, s3, s4 contain the shells to do in libint order
    compute_quartet(s1, s2, s3, s4);

    // Permute integrals back, if needed
    if (p12 || p34 || p13p24) {
#ifdef MINTS_TIMER
        timer_on("permute_target");
#endif
        permute_target(source_, target_, s1, s2, s3, s4, p12, p34, p13p24);
#ifdef MINTS_TIMER
        timer_off("permute_target");
#endif
    }
    else {
#ifdef MINTS_TIMER
        timer_on("memcpy - no resort");
#endif
        // copy the integrals to the target_
        memcpy(target_, source_, n1 * n2 * n3 * n4 *sizeof(double));
#ifdef MINTS_TIMER
        timer_off("memcpy - no resort");
#endif
    }
#ifdef MINTS_TIMER
    timer_off("ERI::compute_shell");
#endif
}

void ERI::compute_quartet(int sh1, int sh2, int sh3, int sh4)
{
#ifdef MINTS_TIMER
    timer_on("setup");
#endif
    shared_ptr<GaussianShell> s1, s2, s3, s4;

    s1 = bs1_->shell(sh1);
    s2 = bs2_->shell(sh2);
    s3 = bs3_->shell(sh3);
    s4 = bs4_->shell(sh4);

    int am1 = s1->am(0);
    int am2 = s2->am(0);
    int am3 = s3->am(0);
    int am4 = s4->am(0);
    int am = am1 + am2 + am3 + am4; // total am
    int nprim1 = s1->nprimitive();
    int nprim2 = s2->nprimitive();
    int nprim3 = s3->nprimitive();
    int nprim4 = s4->nprimitive();
    size_t nprim, nprim_combination = nprim1 * nprim2 * nprim3 * nprim4;
    double A[3], B[3], C[3], D[3];

    A[0] = s1->center()[0];
    A[1] = s1->center()[1];
    A[2] = s1->center()[2];
    B[0] = s2->center()[0];
    B[1] = s2->center()[1];
    B[2] = s2->center()[2];
    C[0] = s3->center()[0];
    C[1] = s3->center()[1];
    C[2] = s3->center()[2];
    D[0] = s4->center()[0];
    D[1] = s4->center()[1];
    D[2] = s4->center()[2];

    // compute intermediates
//    double AB2 = 0.0;
//    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
//    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
//    AB2 += (A[2] - B[2]) * (A[2] - B[2]);
//    double CD2 = 0.0;
//    CD2 += (C[0] - D[0]) * (C[0] - D[0]);
//    CD2 += (C[1] - D[1]) * (C[1] - D[1]);
//    CD2 += (C[2] - D[2]) * (C[2] - D[2]);

    libint_.AB[0] = A[0] - B[0];
    libint_.AB[1] = A[1] - B[1];
    libint_.AB[2] = A[2] - B[2];
    libint_.CD[0] = C[0] - D[0];
    libint_.CD[1] = C[1] - D[1];
    libint_.CD[2] = C[2] - D[2];

#ifdef MINTS_TIMER
    timer_off("setup");
#endif

#ifdef MINTS_TIMER
    timer_on("Primitive setup");
#endif

    // Prepare all the data needed by libint
    nprim = 0;
    // New version - with ShellPair
    ShellPair *p12 = &(pairs12_[sh1][sh2]);
    ShellPair *p34 = &(pairs12_[sh3][sh4]);

    for (int p1=0; p1<nprim1; ++p1) {
        for (int p2=0; p2<nprim2; ++p2) {
            double zeta = p12->gamma[p1][p2];
            double overlap12 = p12->overlap[p1][p2];
            for (int p3=0; p3<nprim3; ++p3) {
                for (int p4=0; p4<nprim4; ++p4) {
                    double eta  = p34->gamma[p3][p4];
                    double oozn = 1.0 / (zeta+eta);
                    libint_.PrimQuartet[nprim].poz = eta * oozn;
                    double rho = zeta * libint_.PrimQuartet[nprim].poz;
//                    double rho = (zeta*eta) * oozn;
                    double coef1 = 2.0 * sqrt(rho*M_1_PI) * overlap12 * p34->overlap[p3][p4];
                    double PQ[3];
                    PQ[0] = p12->P[p1][p2][0] - p34->P[p3][p4][0];
                    PQ[1] = p12->P[p1][p2][1] - p34->P[p3][p4][1];
                    PQ[2] = p12->P[p1][p2][2] - p34->P[p3][p4][2];
                    double PQ2 = PQ[0]*PQ[0] + PQ[1]*PQ[1] + PQ[2]*PQ[2];
                    double T = rho*PQ2;

                    calc_f(libint_.PrimQuartet[nprim].F, am+1, T);

                    // Modify F to include overlap of ab and cd, eqs 14, 15, 16 of libint manual
                    for (int i=0; i<=am; ++i) {
                        libint_.PrimQuartet[nprim].F[i] *= coef1;
                    }

                    libint_.PrimQuartet[nprim].oo2zn = 0.5 * oozn;
                    libint_.PrimQuartet[nprim].pon   = zeta * oozn;
                    libint_.PrimQuartet[nprim].oo2z  = 0.5 / zeta;
                    libint_.PrimQuartet[nprim].oo2n  = 0.5 / eta;
                    double W[3];
                    W[0] = (p12->P[p1][p2][0] * zeta + p34->P[p3][p4][0] * eta) * oozn;
                    W[1] = (p12->P[p1][p2][1] * zeta + p34->P[p3][p4][1] * eta) * oozn;
                    W[2] = (p12->P[p1][p2][2] * zeta + p34->P[p3][p4][2] * eta) * oozn;

                    // PA
                    libint_.PrimQuartet[nprim].U[0][0] = p12->PA[p1][p2][0];
                    libint_.PrimQuartet[nprim].U[0][1] = p12->PA[p1][p2][1];
                    libint_.PrimQuartet[nprim].U[0][2] = p12->PA[p1][p2][2];
                    // QC
                    libint_.PrimQuartet[nprim].U[2][0] = p34->PA[p3][p4][0];
                    libint_.PrimQuartet[nprim].U[2][1] = p34->PA[p3][p4][1];
                    libint_.PrimQuartet[nprim].U[2][2] = p34->PA[p3][p4][2];
                    // WP
                    libint_.PrimQuartet[nprim].U[4][0] = W[0] - p12->P[p1][p2][0];
                    libint_.PrimQuartet[nprim].U[4][1] = W[1] - p12->P[p1][p2][1];
                    libint_.PrimQuartet[nprim].U[4][2] = W[2] - p12->P[p1][p2][2];
                    // WQ
                    libint_.PrimQuartet[nprim].U[5][0] = W[0] - p34->P[p3][p4][0];
                    libint_.PrimQuartet[nprim].U[5][1] = W[1] - p34->P[p3][p4][1];
                    libint_.PrimQuartet[nprim].U[5][2] = W[2] - p34->P[p3][p4][2];

//                    fprintf(outfile, "----- %d\n", nprim);
//                    fprintf(outfile, "poz = %lf\n", libint_.PrimQuartet[nprim].poz);
//                    fprintf(outfile, "T   = %lf\n", T);
//                    fprintf(outfile, "oo2zn = %lf\n", libint_.PrimQuartet[nprim].oo2zn);
//                    fprintf(outfile, "pon = %lf\n", libint_.PrimQuartet[nprim].pon);
//                    fprintf(outfile, "oo2z = %lf\n", libint_.PrimQuartet[nprim].oo2z);
//                    fprintf(outfile, "oo2n = %lf\n", libint_.PrimQuartet[nprim].oo2n);
//                    fprintf(outfile, "U[0] = %lf %lf %lf\n", libint_.PrimQuartet[nprim].U[0][0], libint_.PrimQuartet[nprim].U[0][1], libint_.PrimQuartet[nprim].U[0][2]);
//                    fprintf(outfile, "U[2] = %lf %lf %lf\n", libint_.PrimQuartet[nprim].U[2][0], libint_.PrimQuartet[nprim].U[2][1], libint_.PrimQuartet[nprim].U[2][2]);
//                    fprintf(outfile, "U[4] = %lf %lf %lf\n", libint_.PrimQuartet[nprim].U[4][0], libint_.PrimQuartet[nprim].U[4][1], libint_.PrimQuartet[nprim].U[4][2]);
//                    fprintf(outfile, "U[5] = %lf %lf %lf\n", libint_.PrimQuartet[nprim].U[5][0], libint_.PrimQuartet[nprim].U[5][1], libint_.PrimQuartet[nprim].U[5][2]);
//                    fprintf(outfile, "W = %lf %lf %lf\n", W[0], W[1], W[2]);
//                    fprintf(outfile, "Q = %lf %lf %lf\n", Q[0], Q[1], Q[2]);
                    nprim++;
                }
            }
        }
    }
    // Old version - without ShellPair
//    for (int p1=0; p1<nprim1; ++p1) {
//        double a1 = s1->exp(p1);
//        double c1 = s1->coef(0, p1);
//        for (int p2=0; p2<nprim2; ++p2) {
//            double a2 = s2->exp(p2);
//            double c2 = s2->coef(0, p2);
//            double zeta = a1 + a2;
//            double ooz = 1.0/zeta;
//            double oo2z = 1.0/(2.0 * zeta);
//
//            double PA[3], PB[3];
//            double P[3];
//
//            P[0] = (a1*A[0] + a2*B[0])*ooz;
//            P[1] = (a1*A[1] + a2*B[1])*ooz;
//            P[2] = (a1*A[2] + a2*B[2])*ooz;
//            PA[0] = P[0] - A[0];
//            PA[1] = P[1] - A[1];
//            PA[2] = P[2] - A[2];
//            PB[0] = P[0] - B[0];
//            PB[1] = P[1] - B[1];
//            PB[2] = P[2] - B[2];
//
//            double Sab = pow(M_PI*ooz, 3.0/2.0) * exp(-a1*a2*ooz*AB2) * c1 * c2;
//
//            for (int p3=0; p3<nprim3; ++p3) {
//                double a3 = s3->exp(p3);
//                double c3 = s3->coef(0, p3);
//                for (int p4=0; p4<nprim4; ++p4) {
//                    double a4 = s4->exp(p4);
//                    double c4 = s4->coef(0, p4);
//                    double nu = a3 + a4;
//                    double oon = 1.0/nu;
//                    double oo2n = 1.0/(2.0*nu);
//                    double oo2zn = 1.0/(2.0*(zeta+nu));
//                    double rho = (zeta*nu)/(zeta+nu);
//                    double oo2rho = 1.0 / (2.0*rho);
//
//                    double QC[3], QD[3], WP[3], WQ[3], PQ[3];
//                    double Q[3], W[3], a3C[3], a4D[3];
//
//                    a3C[0] = a3*C[0];
//                    a3C[1] = a3*C[1];
//                    a3C[2] = a3*C[2];
//
//                    a4D[0] = a4*D[0];
//                    a4D[1] = a4*D[1];
//                    a4D[2] = a4*D[2];
//
////                  Q[0] = (a3*C[0] + a4*D[0])*oon;
////                  Q[1] = (a3*C[1] + a4*D[1])*oon;
////                  Q[2] = (a3*C[2] + a4*D[2])*oon;
//
//                    Q[0] = (a3C[0] + a4D[0])*oon;
//                    Q[1] = (a3C[1] + a4D[1])*oon;
//                    Q[2] = (a3C[2] + a4D[2])*oon;
//
//                    QC[0] = Q[0] - C[0];
//                    QC[1] = Q[1] - C[1];
//                    QC[2] = Q[2] - C[2];
//                    QD[0] = Q[0] - D[0];
//                    QD[1] = Q[1] - D[1];
//                    QD[2] = Q[2] - D[2];
//                    PQ[0] = P[0] - Q[0];
//                    PQ[1] = P[1] - Q[1];
//                    PQ[2] = P[2] - Q[2];
//
//                    double PQ2 = 0.0;
//                    PQ2 += (P[0] - Q[0]) * (P[0] - Q[0]);
//                    PQ2 += (P[1] - Q[1]) * (P[1] - Q[1]);
//                    PQ2 += (P[2] - Q[2]) * (P[2] - Q[2]);
//
//                    W[0] = (zeta*P[0] + nu*Q[0]) / (zeta + nu);
//                    W[1] = (zeta*P[1] + nu*Q[1]) / (zeta + nu);
//                    W[2] = (zeta*P[2] + nu*Q[2]) / (zeta + nu);
//                    WP[0] = W[0] - P[0];
//                    WP[1] = W[1] - P[1];
//                    WP[2] = W[2] - P[2];
//                    WQ[0] = W[0] - Q[0];
//                    WQ[1] = W[1] - Q[1];
//                    WQ[2] = W[2] - Q[2];
//
//                    for (int i=0; i<3; ++i) {
//                        libint_.PrimQuartet[nprim].U[0][i] = PA[i];
//                        libint_.PrimQuartet[nprim].U[2][i] = QC[i];
//                        libint_.PrimQuartet[nprim].U[4][i] = WP[i];
//                        libint_.PrimQuartet[nprim].U[5][i] = WQ[i];
//                    }
//                    libint_.PrimQuartet[nprim].oo2z = oo2z;
//                    libint_.PrimQuartet[nprim].oo2n = oo2n;
//                    libint_.PrimQuartet[nprim].oo2zn = oo2zn;
//                    libint_.PrimQuartet[nprim].poz = rho * ooz;
//                    libint_.PrimQuartet[nprim].pon = rho * oon;
//                    libint_.PrimQuartet[nprim].oo2p = oo2rho;
//
//                    double T = rho * PQ2;
//                    calc_f(libint_.PrimQuartet[nprim].F, am+1, T);
//
//                    // Modify F to include overlap of ab and cd, eqs 14, 15, 16 of libint manual
//                    double Scd = pow(M_PI*oon, 3.0/2.0) * exp(-a3*a4*oon*CD2) * c3 * c4;
//                    double val = 2.0 * sqrt(rho * M_1_PI) * Sab * Scd;
//                    for (int i=0; i<=am; ++i) {
//                        libint_.PrimQuartet[nprim].F[i] *= val;
//                    }
//                    fprintf(outfile, "----- %d\n", nprim);
//                    fprintf(outfile, "poz = %lf\n", libint_.PrimQuartet[nprim].poz);
//                    fprintf(outfile, "T   = %lf\n", T);
//                    fprintf(outfile, "oo2zn = %lf\n", libint_.PrimQuartet[nprim].oo2zn);
//                    fprintf(outfile, "pon = %lf\n", libint_.PrimQuartet[nprim].pon);
//                    fprintf(outfile, "oo2z = %lf\n", libint_.PrimQuartet[nprim].oo2z);
//                    fprintf(outfile, "oo2n = %lf\n", libint_.PrimQuartet[nprim].oo2n);
//                    fprintf(outfile, "U[0] = %lf %lf %lf\n", libint_.PrimQuartet[nprim].U[0][0], libint_.PrimQuartet[nprim].U[0][1], libint_.PrimQuartet[nprim].U[0][2]);
//                    fprintf(outfile, "U[2] = %lf %lf %lf\n", libint_.PrimQuartet[nprim].U[2][0], libint_.PrimQuartet[nprim].U[2][1], libint_.PrimQuartet[nprim].U[2][2]);
//                    fprintf(outfile, "U[4] = %lf %lf %lf\n", libint_.PrimQuartet[nprim].U[4][0], libint_.PrimQuartet[nprim].U[4][1], libint_.PrimQuartet[nprim].U[4][2]);
//                    fprintf(outfile, "U[5] = %lf %lf %lf\n", libint_.PrimQuartet[nprim].U[5][0], libint_.PrimQuartet[nprim].U[5][1], libint_.PrimQuartet[nprim].U[5][2]);
//                    fprintf(outfile, "W = %lf %lf %lf\n", W[0], W[1], W[2]);
//                    fprintf(outfile, "Q = %lf %lf %lf\n", Q[0], Q[1], Q[2]);
//                    nprim++;
//                }
//            }
//        }
//    }

#ifdef MINTS_TIMER
    timer_off("Primitive setup");
#endif

    // How many are there?
    size_t size = INT_NCART(am1) * INT_NCART(am2) * INT_NCART(am3) * INT_NCART(am4);

#ifdef MINTS_TIMER
    timer_on("libint overhead");
#endif

    // Compute the integral
    if (am) {
        REALTYPE *target_ints;

        target_ints = build_eri[am1][am2][am3][am4](&libint_, nprim);

        memcpy(source_, target_ints, sizeof(double)*size);
    } else {
        // Handle (ss|ss)
        double temp = 0.0;
        for (size_t i=0; i<nprim_combination; ++i)
            temp += (double)libint_.PrimQuartet[i].F[0];
        source_[0] = temp;
    }

#ifdef MINTS_TIMER
    timer_off("libint overhead");
#endif

    // The following two functions time themselves.

    // Normalize the integrals for angular momentum
    normalize_am(s1, s2, s3, s4);

    // Transform the integrals to the spherical basis
    pure_transform(sh1, sh2, sh3, sh4, 1);

    // Results are in source_
}

void ERI::compute_shell_deriv1(int sh1, int sh2, int sh3, int sh4)
{
    if (deriv_ >= 1) {
        fprintf(stderr, "ERROR - ERI: ERI object not initialized to handle derivatives.\n");
        abort();
    }
    // Need to ensure the ordering asked by the user is valid for libint
    // compute_quartet does NOT check this. SEGFAULTS should occur if order
    // is not guaranteed.
    int s1, s2, s3, s4;
    int am1, am2, am3, am4, temp;
    shared_ptr<BasisSet> bs_temp;
    bool p13p24 = false, p12 = false, p34 = false;

    // AM used for ordering
    am1 = original_bs1_->shell(sh1)->am(0);
    am2 = original_bs2_->shell(sh2)->am(0);
    am3 = original_bs3_->shell(sh3)->am(0);
    am4 = original_bs4_->shell(sh4)->am(0);

    int n1 = original_bs1_->shell(sh1)->nfunction(0);
    int n2 = original_bs2_->shell(sh2)->nfunction(0);
    int n3 = original_bs3_->shell(sh3)->nfunction(0);
    int n4 = original_bs4_->shell(sh4)->nfunction(0);

    // l(a) >= l(b), l(c) >= l(d), and l(c) + l(d) >= l(a) + l(b).
    if (am1 >= am2) {
        s1 = sh1;
        s2 = sh2;

        bs1_ = original_bs1_;
        bs2_ = original_bs2_;
    } else {
        s1 = sh2;
        s2 = sh1;

        bs1_ = original_bs2_;
        bs2_ = original_bs1_;

        p12 = true;
    }

    if (am3 >= am4) {
        s3 = sh3;
        s4 = sh4;

        bs3_ = original_bs3_;
        bs4_ = original_bs4_;

    } else {
        s3 = sh4;
        s4 = sh3;

        bs3_ = original_bs4_;
        bs4_ = original_bs3_;

        p34 = true;
    }

    if ((am1 + am2) > (am3 + am4)) {
        // Swap s1 and s2 with s3 and s4
        temp = s1;
        s1 = s3;
        s3 = temp;

        temp = s2;
        s2 = s4;
        s4 = temp;

        bs_temp = bs1_;
        bs1_ = bs3_;
        bs3_ = bs_temp;

        bs_temp = bs2_;
        bs2_ = bs4_;
        bs4_ = bs_temp;

        p13p24 = true;
    }

    // s1, s2, s3, s4 contain the shells to do in libderive order
    compute_quartet_deriv1(s1, s2, s3, s4);    // compute 9 sets of integral derivatives

    size_t size = n1 * n2 * n3 * n4;
    // Permute integrals back, if needed
    if (p12 || p34 || p13p24) {
        // 3n of them
        for (int i=0; i<3*natom_; ++i)
            permute_target(source_+(i*size), target_+(i*size), s1, s2, s3, s4, p12, p34, p13p24);
    }
    else {
        // copy the integrals to the target_, 3n of them
        memcpy(target_, source_, 3 * natom_ * size *sizeof(double));
    }
}

void ERI::compute_quartet_deriv1(int sh1, int sh2, int sh3, int sh4)
{
    shared_ptr<GaussianShell> s1, s2, s3, s4;

    s1 = bs1_->shell(sh1);
    s2 = bs2_->shell(sh2);
    s3 = bs3_->shell(sh3);
    s4 = bs4_->shell(sh4);

    int am1 = s1->am(0);
    int am2 = s2->am(0);
    int am3 = s3->am(0);
    int am4 = s4->am(0);
    int am = am1 + am2 + am3 + am4; // total am
    int nprim1 = s1->nprimitive();
    int nprim2 = s2->nprimitive();
    int nprim3 = s3->nprimitive();
    int nprim4 = s4->nprimitive();
    size_t nprim;
    double A[3], B[3], C[3], D[3];
    A[0] = s1->center()[0];
    A[1] = s1->center()[1];
    A[2] = s1->center()[2];
    B[0] = s2->center()[0];
    B[1] = s2->center()[1];
    B[2] = s2->center()[2];
    C[0] = s3->center()[0];
    C[1] = s3->center()[1];
    C[2] = s3->center()[2];
    D[0] = s4->center()[0];
    D[1] = s4->center()[1];
    D[2] = s4->center()[2];

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);
    double CD2 = 0.0;
    CD2 += (C[0] - D[0]) * (C[0] - D[0]);
    CD2 += (C[1] - D[1]) * (C[1] - D[1]);
    CD2 += (C[2] - D[2]) * (C[2] - D[2]);

    libderiv_.AB[0] = A[0] - B[0];
    libderiv_.AB[1] = A[1] - B[1];
    libderiv_.AB[2] = A[2] - B[2];
    libderiv_.CD[0] = C[0] - D[0];
    libderiv_.CD[1] = C[1] - D[1];
    libderiv_.CD[2] = C[2] - D[2];

    // Prepare all the data needed by libint
    nprim = 0;
    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1->exp(p1);
        double c1 = s1->coef(0, p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2->exp(p2);
            double c2 = s2->coef(0, p2);
            double zeta = a1 + a2;
            double ooz = 1.0/zeta;
            double oo2z = 1.0/(2.0 * zeta);

            double PA[3], PB[3];
            double P[3];

            P[0] = (a1*A[0] + a2*B[0])*ooz;
            P[1] = (a1*A[1] + a2*B[1])*ooz;
            P[2] = (a1*A[2] + a2*B[2])*ooz;
            PA[0] = P[0] - A[0];
            PA[1] = P[1] - A[1];
            PA[2] = P[2] - A[2];
            PB[0] = P[0] - B[0];
            PB[1] = P[1] - B[1];
            PB[2] = P[2] - B[2];

            double Sab = pow(M_PI*ooz, 3.0/2.0) * exp(-a1*a2*ooz*AB2) * c1 * c2;

            for (int p3=0; p3<nprim3; ++p3) {
                double a3 = s3->exp(p3);
                double c3 = s3->coef(0, p3);
                for (int p4=0; p4<nprim4; ++p4) {
                    double a4 = s4->exp(p4);
                    double c4 = s4->coef(0, p4);
                    double nu = a3 + a4;
                    double oon = 1.0/nu;
                    double oo2n = 1.0/(2.0*nu);
                    double oo2zn = 1.0/(2.0*(zeta+nu));
                    double rho = (zeta*nu)/(zeta+nu);
                    double oo2rho = 1 / (2.0*rho);

                    double QC[3], QD[3], WP[3], WQ[3], PQ[3];
                    double Q[3], W[3];

                    Q[0] = (a3*C[0] + a4*D[0])*oon;
                    Q[1] = (a3*C[1] + a4*D[1])*oon;
                    Q[2] = (a3*C[2] + a4*D[2])*oon;
                    QC[0] = Q[0] - C[0];
                    QC[1] = Q[1] - C[1];
                    QC[2] = Q[2] - C[2];
                    QD[0] = Q[0] - D[0];
                    QD[1] = Q[1] - D[1];
                    QD[2] = Q[2] - D[2];
                    PQ[0] = P[0] - Q[0];
                    PQ[1] = P[1] - Q[1];
                    PQ[2] = P[2] - Q[2];

                    double PQ2 = 0.0;
                    PQ2 += (P[0] - Q[0]) * (P[0] - Q[0]);
                    PQ2 += (P[1] - Q[1]) * (P[1] - Q[1]);
                    PQ2 += (P[2] - Q[2]) * (P[2] - Q[2]);

                    W[0] = (zeta*P[0] + nu*Q[0]) / (zeta + nu);
                    W[1] = (zeta*P[1] + nu*Q[1]) / (zeta + nu);
                    W[2] = (zeta*P[2] + nu*Q[2]) / (zeta + nu);
                    WP[0] = W[0] - P[0];
                    WP[1] = W[1] - P[1];
                    WP[2] = W[2] - P[2];
                    WQ[0] = W[0] - Q[0];
                    WQ[1] = W[1] - Q[1];
                    WQ[2] = W[2] - Q[2];

                    for (int i=0; i<3; ++i) {
                        libderiv_.PrimQuartet[nprim].U[0][i] = PA[i];
                        libderiv_.PrimQuartet[nprim].U[1][i] = PB[i];
                        libderiv_.PrimQuartet[nprim].U[2][i] = QC[i];
                        libderiv_.PrimQuartet[nprim].U[3][i] = QD[i];
                        libderiv_.PrimQuartet[nprim].U[4][i] = WP[i];
                        libderiv_.PrimQuartet[nprim].U[5][i] = WQ[i];
                    }
                    libderiv_.PrimQuartet[nprim].oo2z = oo2z;
                    libderiv_.PrimQuartet[nprim].oo2n = oo2n;
                    libderiv_.PrimQuartet[nprim].oo2zn = oo2zn;
                    libderiv_.PrimQuartet[nprim].poz = rho * ooz;
                    libderiv_.PrimQuartet[nprim].pon = rho * oon;
                    // libderiv_.PrimQuartet[nprim].oo2p = oo2rho;   // NOT SET IN CINTS
                    libderiv_.PrimQuartet[nprim].twozeta_a = 2.0 * a1;
                    libderiv_.PrimQuartet[nprim].twozeta_b = 2.0 * a2;
                    libderiv_.PrimQuartet[nprim].twozeta_c = 2.0 * a3;
                    libderiv_.PrimQuartet[nprim].twozeta_d = 2.0 * a4;

                    double T = rho * PQ2;
                    int_fjt(libderiv_.PrimQuartet[nprim].F, am+DERIV_LVL, T);

                    // Modify F to include overlap of ab and cd, eqs 14, 15, 16 of libint manual
                    double Scd = pow(M_PI*oon, 3.0/2.0) * exp(-a3*a4*oon*CD2) * c3 * c4;
                    double val = 2.0 * sqrt(rho * M_1_PI) * Sab * Scd;
                    for (int i=0; i<=am+DERIV_LVL; ++i) {
                        libderiv_.PrimQuartet[nprim].F[i] *= val;
                    }

                    nprim++;
                }
            }
        }
    }

    // How many are there?
    size_t size = INT_NCART(am1) * INT_NCART(am2) * INT_NCART(am3) * INT_NCART(am4);

    // Compute the integral
    build_deriv1_eri[am1][am2][am3][am4](&libderiv_, nprim);

    // Sum results into derivs_
    int center_i = s1->ncenter()*3*size;
    int center_j = s2->ncenter()*3*size;
    int center_k = s3->ncenter()*3*size;
    int center_l = s4->ncenter()*3*size;

    // Zero out memory
    memset(source_, 0, sizeof(double) * size * 3 * natom_);

    for (size_t i=0; i < size; ++i) {
        source_[center_i+(0*size)+i] += libderiv_.ABCD[0][i];
        source_[center_i+(1*size)+i] += libderiv_.ABCD[1][i];
        source_[center_i+(2*size)+i] += libderiv_.ABCD[2][i];

        // Use translational invariance to determine center_j derivatives
        source_[center_j+(0*size)+i] -= (libderiv_.ABCD[0][i] + libderiv_.ABCD[6][i] + libderiv_.ABCD[9][i]);
        source_[center_j+(1*size)+i] -= (libderiv_.ABCD[1][i] + libderiv_.ABCD[7][i] + libderiv_.ABCD[10][i]);
        source_[center_j+(2*size)+i] -= (libderiv_.ABCD[2][i] + libderiv_.ABCD[8][i] + libderiv_.ABCD[11][i]);

        source_[center_k+(0*size)+i] += libderiv_.ABCD[6][i];
        source_[center_k+(1*size)+i] += libderiv_.ABCD[7][i];
        source_[center_k+(2*size)+i] += libderiv_.ABCD[8][i];

        source_[center_l+(0*size)+i] += libderiv_.ABCD[9][i];
        source_[center_l+(1*size)+i] += libderiv_.ABCD[10][i];
        source_[center_l+(2*size)+i] += libderiv_.ABCD[11][i];
    }

    // Normalize the 3n sets of integrals
    normalize_am(s1, s2, s3, s4, 3*natom_);

    // Transform the integrals to the spherical basis
    pure_transform(sh1, sh2, sh3, sh4, 3*natom_);

    // Results are in source_
}
int ERI::shell_is_zero(int sh1, int sh2, int sh3, int sh4)
{

        if (schwarz2_ != 0.0)
            if (screen_ == false)
                form_sieve();

        //fprintf(outfile,"\nSchwarz val is %f",schwarz_norm_[ioff[((sh1>sh2)?sh1:sh2)]+((sh1>sh2)?sh2:sh1)]*schwarz_norm_[ioff[((sh3>sh4)?sh3:sh4)]+((sh3>sh4)?sh4:sh3)]);
        if (screen_ == true)
        {
            if (schwarz_norm_[ioff[((sh1>sh2)?sh1:sh2)]+((sh1>sh2)?sh2:sh1)]*schwarz_norm_[ioff[((sh3>sh4)?sh3:sh4)]+((sh3>sh4)?sh4:sh3)]<schwarz2_)
            {
                //fprintf(outfile," Shell (%d,%d|%d,%d) is zero",sh1,sh2,sh3,sh4);
                return true;
            }
            else
            {
                //fprintf(outfile," Shell (%d,%d|%d,%d) is nonzero",sh1,sh2,sh3,sh4);
                return false;
            }
        }
        return false;
}

void ERI::form_sieve()
{
    //fprintf(outfile,"Starting Sieve"); fflush(outfile);

  int nshell = original_bs1_->nshell();
  //fprintf(outfile,"Read"); fflush(outfile);
  schwarz_norm_ = init_array(nshell*(nshell+1)/2);

  double cut = schwarz2_;
  schwarz2_ = 0.0;

  double max;
  int MU,NU,numMU,numNU,N,M, MN, ind;
  for (MU = 0, MN = 0; MU < nshell; MU++)
    for (NU = 0; NU <= MU; NU++, MN++)
    {
        compute_shell(MU,NU,MU,NU);
        numMU = original_bs1_->shell(MU)->nfunction();
        numNU = original_bs1_->shell(NU)->nfunction();
        max = 0.0;

        for (M = 0, ind = 0; M<numMU*numMU; M++)
            for (N = 0; N<numNU*numNU; N++, ind++)
                if (fabs(target_[ind])>max)
                    max = fabs(target_[ind]);

        schwarz_norm_[ioff[MU]+NU] = max;
    }
  schwarz2_ = cut;
  screen_ = true;

  //fprintf(outfile,"Norm:\n");
  //for (int ij = 0; ij<nshell*(nshell+1)/2; ij++)
    //fprintf(outfile,"%20.10f\n",schwarz_norm_[ij]);

  //fprintf(outfile,"Iterations done\n"); fflush(outfile);
}

