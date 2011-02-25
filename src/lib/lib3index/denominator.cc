#include "3index.h"

#include <boost/regex.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/algorithm/string.hpp>

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
Denominator::Denominator(shared_ptr<Vector> eps_occ, shared_ptr<Vector> eps_vir, double delta) : 
    eps_occ_(eps_occ), eps_vir_(eps_vir), delta_(delta)
{
}
Denominator::~Denominator()
{
}
LaplaceDenominator::LaplaceDenominator(shared_ptr<Vector> eps_occ, shared_ptr<Vector> eps_vir, double delta) : 
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

    shared_ptr<Matrix> true_denom(new Matrix("Exact Delta Tensor", nocc*nvir, nocc*nvir));
    shared_ptr<Matrix> app_denom(new Matrix("Approximate Delta Tensor", nocc*nvir, nocc*nvir));
    shared_ptr<Matrix> err_denom(new Matrix("Error in Delta Tensor", nocc*nvir, nocc*nvir));

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
    std::string err_table_filename = PSIDATADIR + "quadratures/1_x/error.bin";
    std::string R_filename = PSIDATADIR + "quadratures/1_x/R_avail.bin";

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
 
    shared_ptr<Matrix> err_table(new Matrix("Error Table (nR x nk)", nR, nk));
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

    int exponent = (int) floor(log(R_availp[r])/log(10.0));
    int mantissa = (int) round(R_availp[r]/pow(10.0,exponent));   

    std::stringstream st;
    st << setfill('0');
    st << "1_xk" <<  setw(2) << nvector_;
    st << "_" << mantissa; 
    st << "E" << exponent;  
 
    std::string quadfile = PSIDATADIR + "quadratures/1_x/" + st.str().c_str();

    fprintf(outfile, "\n  ==> Laplace Denominator <==\n\n"); 
    fprintf(outfile, "  This system has an intrinsic R = (E_HUMO - E_LOMO)/(E_LUMO - E_HOMO) of %7.4E.\n", R);
    fprintf(outfile, "  A %d point minimax quadrature with R of %1.0E will be used for the denominator.\n", nvector_, R_availp[r]);
    fprintf(outfile, "  The worst-case Chebyshev norm for this quadrature rule is %7.4E.\n", accuracy); 
    fprintf(outfile, "  Quadrature rule read from file %s.\n", quadfile.c_str());    

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

    denominator_occ_ = shared_ptr<Matrix>(new Matrix("Occupied Laplace Delta Tensor", nvector_, nocc)); 
    denominator_vir_ = shared_ptr<Matrix>(new Matrix("Virtual Laplace Delta Tensor", nvector_, nvir)); 
    denominator_ = shared_ptr<Matrix>(new Matrix("OV Laplace Delta Tensor", nvector_, nocc*nvir)); 

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

    shared_ptr<Matrix> true_denom(new Matrix("Exact Delta Tensor", nocc*nvir, nocc*nvir));
    shared_ptr<Matrix> app_denom(new Matrix("Approximate Delta Tensor (Fully Separated)", nocc*nvir, nocc*nvir));
    shared_ptr<Matrix> err_denom(new Matrix("Error in Delta Tensor", nocc*nvir, nocc*nvir));

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
CholeskyDenominator::CholeskyDenominator(shared_ptr<Vector> eps_occ, shared_ptr<Vector> eps_vir, double delta, double dm) : 
    Denominator(eps_occ, eps_vir, delta), degeneracy_multiplier_(dm)
{
    decompose();
}
CholeskyDenominator::~CholeskyDenominator()
{
}
bool CholeskyDenominator::criteria(std::pair<int, double> a, std::pair<int, double> b)
{
    return a.second < b.second;
}
void CholeskyDenominator::decompose()
{
    int noccA = eps_occ_->dimpi()[0];    
    int nvirA = eps_vir_->dimpi()[0];    
    int nspan = noccA*nvirA; 

    std::vector<std::pair<int, double> > eps(nspan);    
    std::vector<int> degenerate_index(nspan);

    for (int i = 0; i < noccA; i++) {
        for (int a = 0; a < nvirA; a++) {
            eps[i*nvirA + a] = make_pair(i*nvirA + a, eps_vir_->get(0,a) - eps_occ_->get(0,i));
        }
    }

    std::sort(eps.begin(), eps.end(), &CholeskyDenominator::criteria);

    double* w_ia = new double[nspan];
    w_ia[0] = eps[0].second;
    degenerate_index[0] = 0;
    
    int N = 1;
    for (int k = 1; k < nspan; k++) {
        if (fabs(eps[k].second - w_ia[N-1]) > delta_ / degeneracy_multiplier_) {
            w_ia[N] = eps[k].second;
            N++;
        }
        degenerate_index[k] = N - 1;
    }
    
    int Ndelta; double max_error = 0.0;
    for (Ndelta = 1; Ndelta <= N; Ndelta++) {
        max_error = 0.0;
        bool OK = true;
        for (int p = Ndelta; p < N; p++) {
            double D = 1.0 / (2.0 * w_ia[p]);
            double Q = 0.0;
            for (int m = 0; m < Ndelta - 1; m++) {
                Q = (w_ia[p] - w_ia[m]) / (w_ia[p] + w_ia[m]);
                D *= Q * Q;
            }
            if (fabs(D) > max_error)
                max_error = fabs(D);           
 
            //printf("  D(%d)^%d = %14.10E\n", p, Ndelta, D);
            if (fabs(D) > delta_) {
                OK = false;
                break;
            }
        }
        if (OK) 
            break;
    }

    nvector_ = Ndelta;

    shared_ptr<Matrix> L(new Matrix("Cholesky L, Energy Denominator", N, Ndelta));
    double** Lp = L->pointer(0);

    for (int n = 0; n < Ndelta; n++) {
        for (int p = n; p < N; p++) {
            double wn = w_ia[n];
            double wp = w_ia[p];
            double wm = 0.0;
            double Q = 0.0;
            double Lpp = sqrt(2.0*wn) / (wp + wn);
            for (int m = 0; m < n; m++) {
                wm = w_ia[m];
                Q = (wp - wm) / (wp + wm);
                Lpp *= Q;
            }
            Lp[p][n] = Lpp;
        }
    }

    delete[] w_ia;

    fprintf(outfile, "\n  ==> Cholesky Denominator <==\n\n"); 
    fprintf(outfile, "  A %d point partial Cholesky decomposition will be used for the denominator.\n", nvector_);
    fprintf(outfile, "  The worst-case Chebyshev norm for this quadrature rule is %7.4E.\n", max_error); 

    denominator_ = shared_ptr<Matrix>(new Matrix("Cholesky Delta Tensor", Ndelta, nspan));
    double** Lar = denominator_->pointer();   
 
    for (int k = 0; k < nspan; k++) {
        int ind = eps[k].first;
        C_DCOPY(Ndelta, Lp[degenerate_index[k]], 1, &Lar[0][ind], noccA*nvirA);
    }

    //fprintf(outfile, " Lar\n");
    //print_mat(Lar, Ndelta, noccA*nvirA, outfile);
    //fprintf(outfile, " Lbs\n");
    //print_mat(Lbs, Ndelta, noccA*nvirA, outfile);

    //shared_ptr<Matrix> Delta(new Matrix(noccA*nvirA, noccB*nvirB));
    //double** Delt = Delta->get_pointer();

    //C_DGEMM('T','N', noccA*nvirA, noccB*nvirB, Ndelta, 1.0, Lar[0], noccA*nvirA, Lbs[0], noccB*nvirB, 0.0, Delt[0], noccB*nvirB); 
    //Delta->print();
 
}
void CholeskyDenominator::debug()
{
    fprintf(outfile, "\n  DEBUG: Cholesky Denominator. Compound results: \n");
    
    Denominator::debug();
}

}

