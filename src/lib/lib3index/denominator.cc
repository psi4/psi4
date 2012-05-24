#include <boost/regex.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/algorithm/string.hpp>

#include "3index.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace boost;
using namespace std;
using namespace psi;

namespace psi {

// the third parameter of from_string() should be
// one of std::hex, std::dec or std::oct
template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}
Denominator::Denominator(boost::shared_ptr<Vector> eps_occ, boost::shared_ptr<Vector> eps_vir, double delta) :
    eps_occ_(eps_occ), eps_vir_(eps_vir), delta_(delta)
{
}
Denominator::~Denominator()
{
}
boost::shared_ptr<Denominator> Denominator::buildDenominator(const std::string& algorithm,
    boost::shared_ptr<Vector> eps_occ, boost::shared_ptr<Vector> eps_vir, double delta)
{
    Denominator* d;
    if (algorithm == "LAPLACE") {
        d = new LaplaceDenominator(eps_occ, eps_vir, delta);
    } else if (algorithm == "CHOLESKY") {
        d = new CholeskyDenominator(eps_occ, eps_vir, delta);
    } else {
        throw PSIEXCEPTION("Denominator: algorithm is not LAPLACE or CHOLESKY");
    }

    return boost::shared_ptr<Denominator>(d);
}

LaplaceDenominator::LaplaceDenominator(boost::shared_ptr<Vector> eps_occ, boost::shared_ptr<Vector> eps_vir, double delta) :
    Denominator(eps_occ, eps_vir, delta)
{
    decompose();
}
void Denominator::debug()
{
    int nocc = eps_occ_->dimpi()[0];
    int nvir = eps_vir_->dimpi()[0];

    double* e_o = eps_occ_->pointer();
    double* e_v = eps_vir_->pointer();
    double** denp = denominator_->pointer();

    SharedMatrix true_denom(new Matrix("Exact Delta Tensor", nocc*nvir, nocc*nvir));
    SharedMatrix app_denom(new Matrix("Approximate Delta Tensor", nocc*nvir, nocc*nvir));
    SharedMatrix err_denom(new Matrix("Error in Delta Tensor", nocc*nvir, nocc*nvir));

    double** tp = true_denom->pointer();
    double** ap = app_denom->pointer();
    double** ep = err_denom->pointer();

    for (int i = 0; i < nocc; i++)
    for (int a = 0; a < nvir; a++)
    for (int j = 0; j < nocc; j++)
    for (int b = 0; b < nvir; b++)
        tp[i*nvir + a][j*nvir + b] =  1.0 / (e_v[a] + e_v[b] - e_o[i] - e_o[j]);

    for (int alpha = 0; alpha < nvector_; alpha++)
    for (int i = 0; i < nocc; i++)
    for (int a = 0; a < nvir; a++)
    for (int j = 0; j < nocc; j++)
    for (int b = 0; b < nvir; b++)
        ap[i*nvir + a][j*nvir + b] +=  denp[alpha][i*nvir + a]*denp[alpha][j*nvir + b];

    C_DCOPY(nocc*nvir*nocc*nvir, ap[0], 1, ep[0], 1);
    C_DAXPY(nocc*nvir*nocc*nvir, -1.0, tp[0], 1, ep[0], 1);

    true_denom->print();
    app_denom->print();
    err_denom->print();
}
LaplaceDenominator::~LaplaceDenominator()
{
}
void LaplaceDenominator::decompose()
{
    int nocc = eps_occ_->dimpi()[0];
    int nvir = eps_vir_->dimpi()[0];

    double E_LOMO = eps_occ_->get(0, 0);
    double E_HOMO = eps_occ_->get(0, nocc - 1);
    double E_LUMO = eps_vir_->get(0, 0);
    double E_HUMO = eps_vir_->get(0, nvir - 1);

    double A = 2.0*(E_LUMO - E_HOMO);
    double B = 2.0*(E_HUMO - E_LOMO);
    double R = B / A;

    // Pick appropriate quadrature file and read contents
    std::string PSIDATADIR = Process::environment("PSIDATADIR");
    std::string err_table_filename = PSIDATADIR + "/quadratures/1_x/error.bin";
    std::string R_filename = PSIDATADIR + "/quadratures/1_x/R_avail.bin";

    ifstream err_table_file(err_table_filename.c_str(), ios::in | ios::binary);
    ifstream R_avail_file(R_filename.c_str(), ios::in | ios::binary);

    if (!err_table_file)
        throw PSIEXCEPTION("LaplaceQuadrature: Cannot locate error property file for quadrature rules (should be PSIDATADIR/quadratures/1_x/error.bin)");
    if (!R_avail_file)
        throw PSIEXCEPTION("LaplaceQuadrature: Cannot locate R property file for quadrature rules (should be PSIDATADIR/quadratures/1_x/R_avail.bin)");

    int nk = 53;
    int nR = 99;

    // Read in the R available
    double* R_availp = new double[nR];
    R_avail_file.read((char*) R_availp, nR*sizeof(double));

    SharedMatrix err_table(new Matrix("Error Table (nR x nk)", nR, nk));
    double** err_tablep = err_table->pointer();
    err_table_file.read((char*) err_tablep[0], nR*nk*sizeof(double));

    R_avail_file.close();
    err_table_file.close();

    //for (int r2 = 0; r2 < nR; r2++)
    //    fprintf(outfile, "  R[%4d] = %20.14E\n", r2+1, R_availp[r2]);
    //err_table->print();

    int indR;
    for (indR = 0; indR < nR; indR++) {
        if (R < R_availp[indR])
            break;
    }
    if (indR == nR) {
        // TODO: Relax this
        throw PSIEXCEPTION("Laplace Quadrature requested for (E_HUMO - E_LOMO)/(E_LUMO-E_HOMO) > 7.0 * 10^12, quadratures are not designed for this range.");
    }

    double accuracy;
    int k, r;
    bool found = false;
    for (k = 0; k < nk; k++) {
        for (r = indR; r < nR; r++) {
            double err = err_tablep[r][k];
            if (err != 0.0 && err < delta_) {
                accuracy = err;
                found = true;
                break;
            }
        }
        if (found)
            break;
    }


    if (!found) {
        throw PSIEXCEPTION("Laplace Quadrature rule could not be found with specified accuracy for this system");
    }

    nvector_ = k + 1;

    // A bit hacky, but OK
    int exponent = (int) floor(log(R_availp[r])/log(10.0));
    int mantissa = (int) round(R_availp[r]/pow(10.0,exponent));
    if (mantissa == 10) {
        exponent++;
        mantissa = 1;
    }

    std::stringstream st;
    st << setfill('0');
    st << "1_xk" <<  setw(2) << nvector_;
    st << "_" << mantissa;
    st << "E" << exponent;

    std::string quadfile = PSIDATADIR + "/quadratures/1_x/" + st.str().c_str();

    fprintf(outfile, "\n  ==> Laplace Denominator <==\n\n");
    fprintf(outfile, "  This system has an intrinsic R = (E_HUMO - E_LOMO)/(E_LUMO - E_HOMO) of %7.4E.\n", R);
    fprintf(outfile, "  A %d point minimax quadrature with R of %1.0E will be used for the denominator.\n", nvector_, R_availp[r]);
    fprintf(outfile, "  The worst-case Chebyshev norm for this quadrature rule is %7.4E.\n", accuracy);
    fprintf(outfile, "  Quadrature rule read from file %s.\n\n", quadfile.c_str());

    // The quadrature is defined as \omega_v exp(-\alpha_v x) = 1/x
    double* alpha = new double[nvector_];
    double* omega = new double[nvector_];

    std::vector<std::string> lines;
    std::string text;
    ifstream infile(quadfile.c_str());
    if (!infile)
        throw PSIEXCEPTION("LaplaceDenominator: Unable to open quadrature rule file: " + quadfile);
    while (infile.good()) {
        getline(infile, text);
        lines.push_back(text);
    }

#define NUMBER "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"
    regex numberline("^\\s*(" NUMBER ").*");
    smatch what;

    // We'll be rigorous, the files are extremely well defined
    int lineno = 0;
    for (int index = 0; index < nvector_; index++) {
        std::string line  = lines[lineno++];
        if (!regex_match(line, what, numberline))
            throw PSIEXCEPTION("LaplaceDenominator: Unable to read grid file line: \n" + line);
        if (!from_string<double>(omega[index], what[1], std::dec))
            throw PSIEXCEPTION("LaplaceDenominator: Unable to convert grid file line: \n" + line);
    }
    for (int index = 0; index < nvector_; index++) {
        std::string line  = lines[lineno++];
        if (!regex_match(line, what, numberline))
            throw PSIEXCEPTION("LaplaceDenominator: Unable to read grid file line: \n" + line);
        if (!from_string<double>(alpha[index], what[1], std::dec))
            throw PSIEXCEPTION("LaplaceDenominator: Unable to convert grid file line: \n" + line);
    }

    //for (int k = 0; k < nvector_; k++)
    //    printf("  %24.16E, %24.16E\n", omega[k], alpha[k]);

    // Cast weights back to problem size
    for (int k = 0; k < nvector_; k++) {
        alpha[k] /= A;
        omega[k] /= A;
    }

    denominator_occ_ = SharedMatrix(new Matrix("Occupied Laplace Delta Tensor", nvector_, nocc));
    denominator_vir_ = SharedMatrix(new Matrix("Virtual Laplace Delta Tensor", nvector_, nvir));
    denominator_ = SharedMatrix(new Matrix("OV Laplace Delta Tensor", nvector_, nocc*nvir));

    double** dop = denominator_occ_->pointer();
    double** dvp = denominator_vir_->pointer();
    double** dovp = denominator_->pointer();

    double* e_o = eps_occ_->pointer();
    double* e_v = eps_vir_->pointer();

    for (int k = 0; k < nvector_; k++) {
        for (int i = 0; i < nocc; i++) {
            dop[k][i] = pow(omega[k],0.25)*exp(alpha[k]*e_o[i]);
        }
        for (int a = 0; a < nvir; a++) {
            dvp[k][a] = pow(omega[k],0.25)*exp(-alpha[k]*e_v[a]);
        }
        for (int i = 0; i < nocc; i++) {
            for (int a = 0; a < nvir; a++) {
                dovp[k][i*nvir + a] = pow(omega[k],0.5)*exp(-alpha[k]*(e_v[a] - e_o[i]));
            }
        }
    }

    delete[] alpha;
    delete[] omega;
    delete[] R_availp;
}
void LaplaceDenominator::debug()
{
    fprintf(outfile, "\n  DEBUG: Laplace Denominator. Compound results: \n");

    Denominator::debug();

    fprintf(outfile, "\n  DEBUG: Laplace Denominator. Compound results: \n");

    int nocc = eps_occ_->dimpi()[0];
    int nvir = eps_vir_->dimpi()[0];

    double* e_o = eps_occ_->pointer();
    double* e_v = eps_vir_->pointer();
    double** denop = denominator_occ_->pointer();
    double** denvp = denominator_vir_->pointer();

    SharedMatrix true_denom(new Matrix("Exact Delta Tensor", nocc*nvir, nocc*nvir));
    SharedMatrix app_denom(new Matrix("Approximate Delta Tensor (Fully Separated)", nocc*nvir, nocc*nvir));
    SharedMatrix err_denom(new Matrix("Error in Delta Tensor", nocc*nvir, nocc*nvir));

    double** tp = true_denom->pointer();
    double** ap = app_denom->pointer();
    double** ep = err_denom->pointer();

    for (int i = 0; i < nocc; i++)
    for (int a = 0; a < nvir; a++)
    for (int j = 0; j < nocc; j++)
    for (int b = 0; b < nvir; b++)
        tp[i*nvir + a][j*nvir + b] =  1.0 / (e_v[a] + e_v[b] - e_o[i] - e_o[j]);

    for (int alpha = 0; alpha < nvector_; alpha++)
    for (int i = 0; i < nocc; i++)
    for (int a = 0; a < nvir; a++)
    for (int j = 0; j < nocc; j++)
    for (int b = 0; b < nvir; b++)
        ap[i*nvir + a][j*nvir + b] +=  denop[alpha][i]*denop[alpha][j]*denvp[alpha][a]*denvp[alpha][b];

    C_DCOPY(nocc*nvir*nocc*nvir, ap[0], 1, ep[0], 1);
    C_DAXPY(nocc*nvir*nocc*nvir, -1.0, tp[0], 1, ep[0], 1);

    true_denom->print();
    app_denom->print();
    err_denom->print();

}
CholeskyDenominator::CholeskyDenominator(boost::shared_ptr<Vector> eps_occ, boost::shared_ptr<Vector> eps_vir, double delta) :
    Denominator(eps_occ, eps_vir, delta)
{
    decompose();
}
CholeskyDenominator::~CholeskyDenominator()
{
}
void CholeskyDenominator::decompose()
{
    double* eps_occp = eps_occ_->pointer();
    double* eps_virp = eps_vir_->pointer();

    int nocc = eps_occ_->dimpi()[0];
    int nvir = eps_vir_->dimpi()[0];
    int nspan = nocc * nvir;

    double* diagonal = new double[nspan];

    for (int i = 0; i < nocc; i++) {
        for (int a = 0; a < nvir; a++) {
            diagonal[i*nvir + a] = 1.0 / (2.0 * (eps_virp[a] - eps_occp[i]));
        }
    }

    std::vector<double*> L;
    std::vector<int> order;

    int Q = -1;
    nvector_ = 0;
    double max_err;
    while (nvector_ < nspan) {

        max_err = diagonal[0];
        int max_index = 0;

        for (int ia = 0; ia < nspan; ia++) {
            if (max_err <= diagonal[ia]) {
                max_index = ia;
                max_err = diagonal[ia];
            }
        }

        if (fabs(max_err) < delta_)
            break;

        int j = max_index / nvir;
        int b = max_index % nvir;

        Q++;
        nvector_++;
        L.push_back(new double[nspan]);

        ::memset(static_cast<void*>(L[Q]), '\0', nspan*sizeof(double));

        // Find the diagonal
        double LL = sqrt(max_err);

        // Update the vector
        for (int i = 0; i < nocc; i++) {
            for (int a = 0; a < nvir; a++) {
                L[Q][i*nvir + a] = 1.0 / (eps_virp[a] + eps_virp[b] -
                    eps_occp[i] - eps_occp[j]);
            }
        }

        for (int P = 0; P < Q; P++) {
            C_DAXPY(nspan,-L[P][max_index],L[P],1,L[Q],1);
        }

        C_DSCAL(nspan,1.0 / LL, L[Q], 1);

        // Explicitly zero out elements of the vector
        // Which are psychologically upper triangular
        for (int i = 0; i < order.size(); i++)
            L[Q][order[i]] = 0.0;

        // Place the diagonal
        L[Q][max_index] = LL;

        // Update the Schur complement
        for (int i = 0; i < nspan; i++)
            diagonal[i] -= L[Q][i] * L[Q][i];

        // Explicitly zero out elements of the Schur complement
        // Which are already exact, and do not really belong
        // This prevents false selection due to roundoff
        diagonal[max_index] = 0.0;

        // Add the diagonal index to the list of completed indices
        order.push_back(max_index);
    }

    fprintf(outfile, "\n  ==> Cholesky Denominator <==\n\n");
    fprintf(outfile, "  A %d point partial Cholesky decomposition will be used for the denominator.\n", nvector_);
    fprintf(outfile, "  The worst-case Chebyshev norm for this quadrature rule is %7.4E.\n\n", max_err);

    denominator_ = SharedMatrix(new Matrix("Cholesky Delta Tensor", nvector_, nspan));
    double** Lar = denominator_->pointer();

    for (int d = 0; d < nvector_; d++) {
        C_DCOPY(nspan, L[d], 1, Lar[d], 1);
        delete[] L[d];
    }
}
void CholeskyDenominator::debug()
{
    fprintf(outfile, "\n  DEBUG: Cholesky Denominator. Compound results: \n");

    Denominator::debug();
}

SAPTDenominator::SAPTDenominator(boost::shared_ptr<Vector> eps_occA,
  boost::shared_ptr<Vector> eps_virA, boost::shared_ptr<Vector> eps_occB,
  boost::shared_ptr<Vector> eps_virB, double delta, bool debug) :
    eps_occA_(eps_occA), eps_virA_(eps_virA), eps_occB_(eps_occB),
    eps_virB_(eps_virB), delta_(delta), debug_(debug)
{
}
SAPTDenominator::~SAPTDenominator()
{
}
boost::shared_ptr<SAPTDenominator> SAPTDenominator::buildDenominator(const std::string& algorithm,
    boost::shared_ptr<Vector> eps_occA, boost::shared_ptr<Vector> eps_virA,
    boost::shared_ptr<Vector> eps_occB, boost::shared_ptr<Vector> eps_virB,
    double delta, bool debug)
{
    SAPTDenominator* d;
    if (algorithm == "LAPLACE") {
        d = new SAPTLaplaceDenominator(eps_occA, eps_virA, eps_occB, eps_virB, delta, debug);
    } else if (algorithm == "CHOLESKY") {
        d = new SAPTCholeskyDenominator(eps_occA, eps_virA, eps_occB, eps_virB, delta, debug);
    } else {
        throw PSIEXCEPTION("Denominator: algorithm is not LAPLACE or CHOLESKY");
    }

    return boost::shared_ptr<SAPTDenominator>(d);
}
void SAPTDenominator::debug()
{
    fprintf(outfile, "\n  ==> Debug Monomer A Denominator <==\n\n");
    check_denom(eps_occA_,eps_virA_,denominatorA_);
    fprintf(outfile, "\n  ==> Debug Monomer B Denominator <==\n\n");
    check_denom(eps_occB_,eps_virB_,denominatorB_);
}
void SAPTDenominator::check_denom(boost::shared_ptr<Vector> eps_occ,
  boost::shared_ptr<Vector> eps_vir, SharedMatrix denominator)
{
    int nocc = eps_occ->dimpi()[0];
    int nvir = eps_vir->dimpi()[0];

    double* e_o = eps_occ->pointer();
    double* e_v = eps_vir->pointer();
    double** denp = denominator->pointer();

    SharedMatrix true_denom(new Matrix("Exact Delta Tensor", nocc*nvir, nocc*nvir));
    SharedMatrix app_denom(new Matrix("Approximate Delta Tensor", nocc*nvir, nocc*nvir));
    SharedMatrix err_denom(new Matrix("Error in Delta Tensor", nocc*nvir, nocc*nvir));

    double** tp = true_denom->pointer();
    double** ap = app_denom->pointer();
    double** ep = err_denom->pointer();

    for (int i = 0; i < nocc; i++)
    for (int a = 0; a < nvir; a++)
    for (int j = 0; j < nocc; j++)
    for (int b = 0; b < nvir; b++)
        tp[i*nvir + a][j*nvir + b] =  1.0 / (e_v[a] + e_v[b] - e_o[i] - e_o[j]);

    for (int alpha = 0; alpha < nvector_; alpha++)
    for (int i = 0; i < nocc; i++)
    for (int a = 0; a < nvir; a++)
    for (int j = 0; j < nocc; j++)
    for (int b = 0; b < nvir; b++)
        ap[i*nvir + a][j*nvir + b] +=  denp[alpha][i*nvir + a]*denp[alpha][j*nvir + b];

    C_DCOPY(nocc*nvir*nocc*nvir, ap[0], 1, ep[0], 1);
    C_DAXPY(nocc*nvir*nocc*nvir, -1.0, tp[0], 1, ep[0], 1);

    true_denom->print();
    app_denom->print();
    err_denom->print();
}
SAPTLaplaceDenominator::SAPTLaplaceDenominator(boost::shared_ptr<Vector> eps_occA,
  boost::shared_ptr<Vector> eps_virA, boost::shared_ptr<Vector> eps_occB,
  boost::shared_ptr<Vector> eps_virB, double delta, bool debug) :
    SAPTDenominator(eps_occA, eps_virA, eps_occB, eps_virB, delta, debug)
{
    decompose();
}
SAPTLaplaceDenominator::~SAPTLaplaceDenominator()
{
}
void SAPTLaplaceDenominator::decompose()
{
    int noccA = eps_occA_->dimpi()[0];
    int nvirA = eps_virA_->dimpi()[0];
    int noccB = eps_occB_->dimpi()[0];
    int nvirB = eps_virB_->dimpi()[0];

    double E_LOMOA = eps_occA_->get(0, 0);
    double E_HOMOA = eps_occA_->get(0, noccA - 1);
    double E_LUMOA = eps_virA_->get(0, 0);
    double E_HUMOA = eps_virA_->get(0, nvirA - 1);
    double E_LOMOB = eps_occB_->get(0, 0);
    double E_HOMOB = eps_occB_->get(0, noccB - 1);
    double E_LUMOB = eps_virB_->get(0, 0);
    double E_HUMOB = eps_virB_->get(0, nvirB - 1);

    double A = (E_LUMOA - E_HOMOA) + (E_LUMOB - E_HOMOB);
    double B = (E_HUMOA - E_LOMOA) + (E_HUMOB - E_LOMOB);
    double R = B / A;

    // Pick appropriate quadrature file and read contents
    std::string PSIDATADIR = Process::environment("PSIDATADIR");
    std::string err_table_filename = PSIDATADIR + "/quadratures/1_x/error.bin";
    std::string R_filename = PSIDATADIR + "/quadratures/1_x/R_avail.bin";

    ifstream err_table_file(err_table_filename.c_str(), ios::in | ios::binary);
    ifstream R_avail_file(R_filename.c_str(), ios::in | ios::binary);

    if (!err_table_file)
        throw PSIEXCEPTION("LaplaceQuadrature: Cannot locate error property file for quadrature rules (should be PSIDATADIR/quadratures/1_x/error.bin)");
    if (!R_avail_file)
        throw PSIEXCEPTION("LapalceQuadrature: Cannot locate R property file for quadrature rules (should be PSIDATADIR/quadratures/1_x/R_avail.bin)");

    int nk = 53;
    int nR = 99;

    // Read in the R available
    double* R_availp = new double[nR];
    R_avail_file.read((char*) R_availp, nR*sizeof(double));

    SharedMatrix err_table(new Matrix("Error Table (nR x nk)", nR, nk));
    double** err_tablep = err_table->pointer();
    err_table_file.read((char*) err_tablep[0], nR*nk*sizeof(double));

    R_avail_file.close();
    err_table_file.close();

    //for (int r2 = 0; r2 < nR; r2++)
    //    fprintf(outfile, "  R[%4d] = %20.14E\n", r2+1, R_availp[r2]);
    //err_table->print();

    int indR;
    for (indR = 0; indR < nR; indR++) {
        if (R < R_availp[indR])
            break;
    }
    if (indR == nR) {
        // TODO: Relax this
        throw PSIEXCEPTION("Laplace Quadrature requested for (E_HUMO - E_LOMO)/(E_LUMO-E_HOMO) > 7.0 * 10^12, quadratures are not designed for this range.");
    }

    double accuracy;
    int k, r;
    bool found = false;
    for (k = 0; k < nk; k++) {
        for (r = indR; r < nR; r++) {
            double err = err_tablep[r][k];
            if (err != 0.0 && err < delta_) {
                accuracy = err;
                found = true;
                break;
            }
        }
        if (found)
            break;
    }


    if (!found) {
        throw PSIEXCEPTION("Laplace Quadrature rule could not be found with specified accuracy for this system");
    }

    nvector_ = k + 1;

                // A bit hacky, but OK
    int exponent = (int) floor(log(R_availp[r])/log(10.0));
    int mantissa = (int) round(R_availp[r]/pow(10.0,exponent));
                if (mantissa == 10) {
                          exponent++;
        mantissa = 1;
                }

    std::stringstream st;
    st << setfill('0');
    st << "1_xk" <<  setw(2) << nvector_;
    st << "_" << mantissa;
    st << "E" << exponent;

    std::string quadfile = PSIDATADIR + "/quadratures/1_x/" + st.str().c_str();

    if (debug_) {
        fprintf(outfile, "\n  ==> Laplace Denominator <==\n\n");
        fprintf(outfile, "  This system has an intrinsic R = (E_HUMO - E_LOMO)/(E_LUMO - E_HOMO) of %7.4E.\n", R);
        fprintf(outfile, "  A %d point minimax quadrature with R of %1.0E will be used for the denominator.\n", nvector_, R_availp[r]);
        fprintf(outfile, "  The worst-case Chebyshev norm for this quadrature rule is %7.4E.\n", accuracy);
        fprintf(outfile, "  Quadrature rule read from file %s.\n", quadfile.c_str());
    }

    // The quadrature is defined as \omega_v exp(-\alpha_v x) = 1/x
    double* alpha = new double[nvector_];
    double* omega = new double[nvector_];

    std::vector<std::string> lines;
    std::string text;
    ifstream infile(quadfile.c_str());
    if (!infile)
        throw PSIEXCEPTION("LaplaceDenominator: Unable to open quadrature rule file: " + quadfile);
    while (infile.good()) {
        getline(infile, text);
        lines.push_back(text);
    }

#define NUMBER "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"
    regex numberline("^\\s*(" NUMBER ").*");
    smatch what;

    // We'll be rigorous, the files are extremely well defined
    int lineno = 0;
    for (int index = 0; index < nvector_; index++) {
        std::string line  = lines[lineno++];
        if (!regex_match(line, what, numberline))
            throw PSIEXCEPTION("LaplaceDenominator: Unable to read grid file line: \n" + line);
        if (!from_string<double>(omega[index], what[1], std::dec))
            throw PSIEXCEPTION("LaplaceDenominator: Unable to convert grid file line: \n" + line);
    }
    for (int index = 0; index < nvector_; index++) {
        std::string line  = lines[lineno++];
        if (!regex_match(line, what, numberline))
            throw PSIEXCEPTION("LaplaceDenominator: Unable to read grid file line: \n" + line);
        if (!from_string<double>(alpha[index], what[1], std::dec))
            throw PSIEXCEPTION("LaplaceDenominator: Unable to convert grid file line: \n" + line);
    }

    //for (int k = 0; k < nvector_; k++)
    //    printf("  %24.16E, %24.16E\n", omega[k], alpha[k]);

    // Cast weights back to problem size
    for (int k = 0; k < nvector_; k++) {
        alpha[k] /= A;
        omega[k] /= A;
    }

    denominator_occA_ = SharedMatrix(new Matrix("Occupied Laplace Delta Tensor (A)", nvector_, noccA));
    denominator_virA_ = SharedMatrix(new Matrix("Virtual Laplace Delta Tensor (A)", nvector_, nvirA));
    denominatorA_ = SharedMatrix(new Matrix("OV Laplace Delta Tensor (A)", nvector_, noccA*nvirA));

    denominator_occB_ = SharedMatrix(new Matrix("Occupied Laplace Delta Tensor (B)", nvector_, noccB));
    denominator_virB_ = SharedMatrix(new Matrix("Virtual Laplace Delta Tensor (B)", nvector_, nvirB));
    denominatorB_ = SharedMatrix(new Matrix("OV Laplace Delta Tensor (B)", nvector_, noccB*nvirB));

    double** doA = denominator_occA_->pointer();
    double** dvA = denominator_virA_->pointer();
    double** dovA = denominatorA_->pointer();

    double** doB = denominator_occB_->pointer();
    double** dvB = denominator_virB_->pointer();
    double** dovB = denominatorB_->pointer();

    double* e_oA = eps_occA_->pointer();
    double* e_vA = eps_virA_->pointer();

    double* e_oB = eps_occB_->pointer();
    double* e_vB = eps_virB_->pointer();

    for (int k = 0; k < nvector_; k++) {
        for (int i = 0; i < noccA; i++) {
            doA[k][i] = pow(omega[k],0.25)*exp(alpha[k]*e_oA[i]);
        }
        for (int a = 0; a < nvirA; a++) {
            dvA[k][a] = pow(omega[k],0.25)*exp(-alpha[k]*e_vA[a]);
        }
        for (int i = 0; i < noccA; i++) {
            for (int a = 0; a < nvirA; a++) {
                dovA[k][i*nvirA + a] = pow(omega[k],0.5)*exp(-alpha[k]*(e_vA[a] - e_oA[i]));
            }
        }

        for (int i = 0; i < noccB; i++) {
            doB[k][i] = pow(omega[k],0.25)*exp(alpha[k]*e_oB[i]);
        }
        for (int a = 0; a < nvirB; a++) {
            dvB[k][a] = pow(omega[k],0.25)*exp(-alpha[k]*e_vB[a]);
        }
        for (int i = 0; i < noccB; i++) {
            for (int a = 0; a < nvirB; a++) {
                dovB[k][i*nvirB + a] = pow(omega[k],0.5)*exp(-alpha[k]*(e_vB[a] - e_oB[i]));
            }
        }
    }


    delete[] alpha;
    delete[] omega;
    delete[] R_availp;
}
void SAPTLaplaceDenominator::debug()
{
    SAPTDenominator::debug();
    fprintf(outfile, "\n  ==> Debug Monomer A Split Denominator <==\n\n");
    check_split(eps_occA_,eps_virA_,denominator_occA_,denominator_virA_);
    fprintf(outfile, "\n  ==> Debug Monomer B Split Denominator <==\n\n");
    check_split(eps_occB_,eps_virB_,denominator_occB_,denominator_virB_);
}
void SAPTLaplaceDenominator::check_split(boost::shared_ptr<Vector> eps_occ,
  boost::shared_ptr<Vector> eps_vir, SharedMatrix denominator_occ,
  SharedMatrix denominator_vir)
{
    int nocc = eps_occ->dimpi()[0];
    int nvir = eps_vir->dimpi()[0];

    double* e_o = eps_occ->pointer();
    double* e_v = eps_vir->pointer();
    double** denop = denominator_occ->pointer();
    double** denvp = denominator_vir->pointer();

    SharedMatrix true_denom(new Matrix("Exact Delta Tensor", nocc*nvir, nocc*nvir));
    SharedMatrix app_denom(new Matrix("Approximate Delta Tensor (Fully Separated)", nocc*nvir, nocc*nvir));
    SharedMatrix err_denom(new Matrix("Error in Delta Tensor", nocc*nvir, nocc*nvir));

    double** tp = true_denom->pointer();
    double** ap = app_denom->pointer();
    double** ep = err_denom->pointer();

    for (int i = 0; i < nocc; i++)
    for (int a = 0; a < nvir; a++)
    for (int j = 0; j < nocc; j++)
    for (int b = 0; b < nvir; b++)
        tp[i*nvir + a][j*nvir + b] =  1.0 / (e_v[a] + e_v[b] - e_o[i] - e_o[j]);

    for (int alpha = 0; alpha < nvector_; alpha++)
    for (int i = 0; i < nocc; i++)
    for (int a = 0; a < nvir; a++)
    for (int j = 0; j < nocc; j++)
    for (int b = 0; b < nvir; b++)
        ap[i*nvir + a][j*nvir + b] +=  denop[alpha][i]*denop[alpha][j]*denvp[alpha][a]*denvp[alpha][b];

    C_DCOPY(nocc*nvir*nocc*nvir, ap[0], 1, ep[0], 1);
    C_DAXPY(nocc*nvir*nocc*nvir, -1.0, tp[0], 1, ep[0], 1);

    true_denom->print();
    app_denom->print();
    err_denom->print();
}
SAPTCholeskyDenominator::SAPTCholeskyDenominator(boost::shared_ptr<Vector> eps_occA,
  boost::shared_ptr<Vector> eps_virA, boost::shared_ptr<Vector> eps_occB,
  boost::shared_ptr<Vector> eps_virB, double delta, bool debug) :
    SAPTDenominator(eps_occA, eps_virA, eps_occB, eps_virB, delta, debug)
{
    decompose();
}
SAPTCholeskyDenominator::~SAPTCholeskyDenominator()
{
}
void SAPTCholeskyDenominator::decompose()
{
    int noccA = eps_occA_->dimpi()[0];
    int nvirA = eps_virA_->dimpi()[0];
    int noccB = eps_occB_->dimpi()[0];
    int nvirB = eps_virB_->dimpi()[0];

    double* eps_occAp = eps_occA_->pointer();
    double* eps_occBp = eps_occB_->pointer();
    double* eps_virAp = eps_virA_->pointer();
    double* eps_virBp = eps_virB_->pointer();

    // Build the schur complement
    boost::shared_ptr<Vector> schurA(new Vector("Diagonal Complement A", noccA * nvirA));
    boost::shared_ptr<Vector> schurB(new Vector("Diagonal Complement B", noccB * nvirB));
    double* schurAp = schurA->pointer();
    double* schurBp = schurB->pointer();

    for (int i = 0; i < noccA; i++) {
        for (int a = 0; a < nvirA; a++) {
            schurAp[i * nvirA + a] = 1.0 / (2.0 * (eps_virAp[a] - eps_occAp[i]));
        }
    }

    for (int j = 0; j < noccB; j++) {
        for (int b = 0; b < nvirB; b++) {
            schurBp[j * nvirB + b] = 1.0 / (2.0 * (eps_virBp[b] - eps_occBp[j]));
        }
    }

    std::vector<double*> denA;
    std::vector<double*> denB;

    std::vector<std::pair<bool,int> > w_order;

    nvector_ = 0;
    int Q = -1;

    double max_err = 0.0;

    while (nvector_ < noccA * nvirA + noccB * nvirB) {
        Q++;
        nvector_++;

        // Locate pivot
        double maxA = 0.0;
        int indA = 0;
        for (int ia = 0; ia < noccA * nvirA; ia++) {
            if (fabs(maxA) < fabs(schurAp[ia])) {
                     maxA  = fabs(schurAp[ia]);
                     indA  = ia;
            }
        }
        double maxB = 0.0;
        int indB = 0;
        for (int ia = 0; ia < noccB * nvirB; ia++) {
            if (fabs(maxA) < fabs(schurBp[ia])) {
                     maxB  = fabs(schurBp[ia]);
                     indB  = ia;
            }
        }

        if (maxB > maxA) {
            w_order.push_back(make_pair(true,indB));
        } else {
            w_order.push_back(make_pair(false,indA));
        }

        bool onB = w_order[Q].first;
        int Q_global = w_order[Q].second;
        int i_global = Q_global / (onB ? nvirB : nvirA);
        int a_global = Q_global % (onB ? nvirB : nvirA);

        // Build row
        denA.push_back(new double[noccA * nvirA]);
        denB.push_back(new double[noccB * nvirB]);
        std::vector<double> L_PQ;
        double L_QQ;
        if (!onB) {
            L_QQ = sqrt(schurAp[Q_global]);

            for (int i = 0; i < noccA; i++) {
                for (int a = 0; a < nvirA; a++) {
                    denA[Q][i * nvirA + a] = 1.0 / (eps_virAp[a_global] - eps_occAp[i_global] + eps_virAp[a] - eps_occAp[i]);
                }
            }

            for (int i = 0; i < noccB; i++) {
                for (int a = 0; a < nvirB; a++) {
                    denB[Q][i * nvirB + a] = 1.0 / (eps_virAp[a_global] - eps_occAp[i_global] + eps_virBp[a] - eps_occBp[i]);
                }
            }

            for (int P = 0; P < Q; P++) {
                L_PQ.push_back(denA[P][Q_global]);
            }

        } else {
            L_QQ = sqrt(schurBp[Q_global]);

            for (int i = 0; i < noccA; i++) {
                for (int a = 0; a < nvirA; a++) {
                    denA[Q][i * nvirA + a] = 1.0 / (eps_virBp[a_global] - eps_occBp[i_global] + eps_virAp[a] - eps_occAp[i]);
                }
            }

            for (int i = 0; i < noccB; i++) {
                for (int a = 0; a < nvirB; a++) {
                    denB[Q][i * nvirB + a] = 1.0 / (eps_virBp[a_global] - eps_occBp[i_global] + eps_virBp[a] - eps_occBp[i]);
                }
            }

            for (int P = 0; P < Q; P++) {
                L_PQ.push_back(denB[P][Q_global]);
            }
        }

        for (int P = 0; P < Q; P++) {
            C_DAXPY(noccA * nvirA, -L_PQ[P], denA[P], 1, denA[Q], 1);
            C_DAXPY(noccB * nvirB, -L_PQ[P], denB[P], 1, denB[Q], 1);
        }

        C_DSCAL(noccA * nvirA, 1.0 / L_QQ, denA[Q], 1);
        C_DSCAL(noccB * nvirB, 1.0 / L_QQ, denB[Q], 1);

        for (int P = 0; P < Q; P++) {
            bool PonB = w_order[P].first;
            int  Pind = w_order[P].second;
            if (!PonB) {
                denA[Q][Pind] = 0.0;
            } else {
                denB[Q][Pind] = 0.0;
            }
        }

        if (!onB) {
            denA[Q][Q_global] = L_QQ;
        } else {
            denB[Q][Q_global] = L_QQ;
        }

        // Update schur complement diagonal
        for (int ia = 0; ia < noccA * nvirA; ia++) {
            schurAp[ia] -= denA[Q][ia] * denA[Q][ia];
        }

        for (int ia = 0; ia < noccB * nvirB; ia++) {
            schurBp[ia] -= denB[Q][ia] * denB[Q][ia];
        }

        for (int P = 0; P <= Q; P++) {
            bool PonB = w_order[P].first;
            int  Pind = w_order[P].second;
            if (!PonB) {
                schurAp[Pind] = 0.0;
            } else {
                schurBp[Pind] = 0.0;
            }
        }

        // Termination condition (schwarz)
        maxA = 0.0;
        for (int ia = 0; ia < noccA * nvirA; ia++) {
            if (fabs(maxA) < fabs(schurAp[ia]))
                     maxA  = fabs(schurAp[ia]);
        }
        maxB = 0.0;
        for (int ia = 0; ia < noccB * nvirB; ia++) {
            if (fabs(maxB) < fabs(schurBp[ia]))
                     maxB  = fabs(schurBp[ia]);
        }

        max_err = sqrt(maxA * maxB);

        if (max_err < delta_) {
            break;
        }
    }

    // Copy Cholesky vectors into permanent matrices
    denominatorA_ = SharedMatrix(new Matrix("Denominator A", nvector_, noccA * nvirA));
    denominatorB_ = SharedMatrix(new Matrix("Denominator B", nvector_, noccB * nvirB));
    double** denAp = denominatorA_->pointer();
    double** denBp = denominatorB_->pointer();

    for (int P = 0; P < nvector_; P++) {
        ::memcpy(static_cast<void*>(denAp[P]), static_cast<void*>(denA[P]), noccA * nvirA * sizeof(double));
        ::memcpy(static_cast<void*>(denBp[P]), static_cast<void*>(denB[P]), noccB * nvirB * sizeof(double));

        delete[] denA[P];
        delete[] denB[P];
    }

    if (debug_) {
        fprintf(outfile, "\n  ==> Cholesky Denominator <==\n\n");
        fprintf(outfile, "  A %d point partial Cholesky decomposition will be used for the denominator.\n", nvector_);
        fprintf(outfile, "  The worst-case Chebyshev norm for this quadrature rule is %7.4E.\n\n", max_err);
    }
}

TLaplaceDenominator::TLaplaceDenominator(boost::shared_ptr<Vector> eps_occ, boost::shared_ptr<Vector> eps_vir, double delta) :
    eps_occ_(eps_occ), eps_vir_(eps_vir), delta_(delta)
{
    decompose();
}
TLaplaceDenominator::~TLaplaceDenominator()
{
}
void TLaplaceDenominator::decompose()
{
    int nocc = eps_occ_->dimpi()[0];
    int nvir = eps_vir_->dimpi()[0];

    double E_LOMO = eps_occ_->get(0, 0);
    double E_HOMO = eps_occ_->get(0, nocc - 1);
    double E_LUMO = eps_vir_->get(0, 0);
    double E_HUMO = eps_vir_->get(0, nvir - 1);

    double A = 3.0*(E_LUMO - E_HOMO);
    double B = 3.0*(E_HUMO - E_LOMO);
    double R = B / A;

    // Pick appropriate quadrature file and read contents
    std::string PSIDATADIR = Process::environment("PSIDATADIR");
    std::string err_table_filename = PSIDATADIR + "/quadratures/1_x/error.bin";
    std::string R_filename = PSIDATADIR + "/quadratures/1_x/R_avail.bin";

    ifstream err_table_file(err_table_filename.c_str(), ios::in | ios::binary);
    ifstream R_avail_file(R_filename.c_str(), ios::in | ios::binary);

    if (!err_table_file)
        throw PSIEXCEPTION("LaplaceQuadrature: Cannot locate error property file for quadrature rules (should be PSIDATADIR/quadratures/1_x/error.bin)");
    if (!R_avail_file)
        throw PSIEXCEPTION("LaplaceQuadrature: Cannot locate R property file for quadrature rules (should be PSIDATADIR/quadratures/1_x/R_avail.bin)");

    int nk = 53;
    int nR = 99;

    // Read in the R available
    double* R_availp = new double[nR];
    R_avail_file.read((char*) R_availp, nR*sizeof(double));

    SharedMatrix err_table(new Matrix("Error Table (nR x nk)", nR, nk));
    double** err_tablep = err_table->pointer();
    err_table_file.read((char*) err_tablep[0], nR*nk*sizeof(double));

    R_avail_file.close();
    err_table_file.close();

    //for (int r2 = 0; r2 < nR; r2++)
    //    fprintf(outfile, "  R[%4d] = %20.14E\n", r2+1, R_availp[r2]);
    //err_table->print();

    int indR;
    for (indR = 0; indR < nR; indR++) {
        if (R < R_availp[indR])
            break;
    }
    if (indR == nR) {
        // TODO: Relax this
        throw PSIEXCEPTION("Laplace Quadrature requested for (E_HUMO - E_LOMO)/(E_LUMO-E_HOMO) > 7.0 * 10^12, quadratures are not designed for this range.");
    }

    double accuracy;
    int k, r;
    bool found = false;
    for (k = 0; k < nk; k++) {
        for (r = indR; r < nR; r++) {
            double err = err_tablep[r][k];
            if (err != 0.0 && err < delta_) {
                accuracy = err;
                found = true;
                break;
            }
        }
        if (found)
            break;
    }


    if (!found) {
        throw PSIEXCEPTION("Laplace Quadrature rule could not be found with specified accuracy for this system");
    }

    nvector_ = k + 1;

    // A bit hacky, but OK
    int exponent = (int) floor(log(R_availp[r])/log(10.0));
    int mantissa = (int) round(R_availp[r]/pow(10.0,exponent));
    if (mantissa == 10) {
        exponent++;
        mantissa = 1;
    }

    std::stringstream st;
    st << setfill('0');
    st << "1_xk" <<  setw(2) << nvector_;
    st << "_" << mantissa;
    st << "E" << exponent;

    std::string quadfile = PSIDATADIR + "/quadratures/1_x/" + st.str().c_str();

    fprintf(outfile, "\n  ==> (T) Laplace Denominator <==\n\n");
    fprintf(outfile, "  This system has an intrinsic R = (E_HUMO - E_LOMO)/(E_LUMO - E_HOMO) of %7.4E.\n", R);
    fprintf(outfile, "  A %d point minimax quadrature with R of %1.0E will be used for the denominator.\n", nvector_, R_availp[r]);
    fprintf(outfile, "  The worst-case Chebyshev norm for this quadrature rule is %7.4E.\n", accuracy);
    fprintf(outfile, "  Quadrature rule read from file %s.\n\n", quadfile.c_str());

    // The quadrature is defined as \omega_v exp(-\alpha_v x) = 1/x
    double* alpha = new double[nvector_];
    double* omega = new double[nvector_];

    std::vector<std::string> lines;
    std::string text;
    ifstream infile(quadfile.c_str());
    if (!infile)
        throw PSIEXCEPTION("LaplaceDenominator: Unable to open quadrature rule file: " + quadfile);
    while (infile.good()) {
        getline(infile, text);
        lines.push_back(text);
    }

#define NUMBER "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"
    regex numberline("^\\s*(" NUMBER ").*");
    smatch what;

    // We'll be rigorous, the files are extremely well defined
    int lineno = 0;
    for (int index = 0; index < nvector_; index++) {
        std::string line  = lines[lineno++];
        if (!regex_match(line, what, numberline))
            throw PSIEXCEPTION("LaplaceDenominator: Unable to read grid file line: \n" + line);
        if (!from_string<double>(omega[index], what[1], std::dec))
            throw PSIEXCEPTION("LaplaceDenominator: Unable to convert grid file line: \n" + line);
    }
    for (int index = 0; index < nvector_; index++) {
        std::string line  = lines[lineno++];
        if (!regex_match(line, what, numberline))
            throw PSIEXCEPTION("LaplaceDenominator: Unable to read grid file line: \n" + line);
        if (!from_string<double>(alpha[index], what[1], std::dec))
            throw PSIEXCEPTION("LaplaceDenominator: Unable to convert grid file line: \n" + line);
    }

    //for (int k = 0; k < nvector_; k++)
    //    printf("  %24.16E, %24.16E\n", omega[k], alpha[k]);

    // Cast weights back to problem size
    for (int k = 0; k < nvector_; k++) {
        alpha[k] /= A;
        omega[k] /= A;
    }

    denominator_occ_ = SharedMatrix(new Matrix("Occupied Laplace Delta Tensor", nvector_, nocc));
    denominator_vir_ = SharedMatrix(new Matrix("Virtual Laplace Delta Tensor", nvector_, nvir));

    double** dop = denominator_occ_->pointer();
    double** dvp = denominator_vir_->pointer();

    double* e_o = eps_occ_->pointer();
    double* e_v = eps_vir_->pointer();

    for (int k = 0; k < nvector_; k++) {
        for (int i = 0; i < nocc; i++) {
            dop[k][i] = pow(omega[k],1.0/6.0)*exp(alpha[k]*e_o[i]);
        }
        for (int a = 0; a < nvir; a++) {
            dvp[k][a] = pow(omega[k],1.0/6.0)*exp(-alpha[k]*e_v[a]);
        }
    }

    delete[] alpha;
    delete[] omega;
    delete[] R_availp;
}
void TLaplaceDenominator::debug()
{
    int nocc = eps_occ_->dimpi()[0];
    int nvir = eps_vir_->dimpi()[0];

    double* e_o = eps_occ_->pointer();
    double* e_v = eps_vir_->pointer();

    double** d_o = denominator_occ_->pointer();
    double** d_v = denominator_vir_->pointer();

    SharedMatrix true_denom(new Matrix("Exact Delta Tensor",      nocc*nocc*nocc,nvir*nvir*nvir));
    SharedMatrix app_denom(new Matrix("Approximate Delta Tensor", nocc*nocc*nocc,nvir*nvir*nvir));
    SharedMatrix err_denom(new Matrix("Error in Delta Tensor",    nocc*nocc*nocc,nvir*nvir*nvir));

    double** tp = true_denom->pointer();
    double** ap = app_denom->pointer();
    double** ep = err_denom->pointer();

    for (int i = 0; i < nocc; i++)
    for (int j = 0; j < nocc; j++)
    for (int k = 0; k < nocc; k++)
    for (int a = 0; a < nvir; a++)
    for (int b = 0; b < nvir; b++)
    for (int c = 0; c < nvir; c++)
        tp[i*nocc*nocc + j*nocc + k][a*nvir*nvir + b*nvir + c] =  1.0 / (e_v[a] + e_v[b] + e_v[c] - e_o[i] - e_o[j] - e_o[k]);

    for (int alpha = 0; alpha < nvector_; alpha++)
    for (int i = 0; i < nocc; i++)
    for (int j = 0; j < nocc; j++)
    for (int k = 0; k < nocc; k++)
    for (int a = 0; a < nvir; a++)
    for (int b = 0; b < nvir; b++)
    for (int c = 0; c < nvir; c++)
        ap[i*nocc*nocc + j*nocc + k][a*nvir*nvir + b*nvir + c] += d_o[alpha][i] * d_o[alpha][j] * d_o[alpha][k] *
                                                                  d_v[alpha][a] * d_v[alpha][b] * d_v[alpha][c]; 

    err_denom->copy(app_denom);
    err_denom->subtract(true_denom);

    denominator_occ_->print();
    denominator_vir_->print();

    true_denom->print();
    app_denom->print();
    err_denom->print();
}

} // Namespace psi
