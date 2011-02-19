#include "3index.h"

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

Denominator::Denominator(shared_ptr<Vector> eps_occ, shared_ptr<Vector> eps_vir) : 
    eps_occ_(eps_occ), eps_vir_(eps_vir)
{
    gauge_ = 0.5 * (eps_occ_->get(0, eps_occ_->dimpi()[0] - 1) + eps_vir_->get(0, 0));
}
Denominator::~Denominator()
{
}
LaplaceDenominator::LaplaceDenominator(shared_ptr<Vector> eps_occ, shared_ptr<Vector> eps_vir, int nvector) : 
    Denominator(eps_occ, eps_vir), nvector_(nvector)
{
    decompose();
}
LaplaceDenominator::~LaplaceDenominator()
{
}
void LaplaceDenominator::decompose()
{
    //fprintf(outfile, "  Energy Denominator Laplace Decomposition\n");
    //fprintf(outfile, "   Number of vectors required = %d [vectors]\n", npoint);
    //fprintf(outfile, "   Memory required = %ld [bytes]\n\n", npoint*nspan*8L);
    //fprintf(outfile, "\n");
   /** 
    double** Lar = block_matrix(npoint, noccA*nvirA); 
    double** Lbs = block_matrix(npoint, noccB*nvirB); 

    shared_ptr<ChebyshevIIQuadrature> quad(new ChebyshevIIQuadrature(npoint, center));
    quad->print();  

    int vec = 0;
    for (quad->reset(); !quad->isDone(); quad->nextPoint()) {
        double t = quad->getPoint();
        double w = sqrt(quad->getWeight()); 
        for (int i = 0; i < noccA; i++) {
            for (int a = 0; a < nvirA; a++) {
                Lar[vec][i * nvirA + a] = w * exp((calc_info_.evalsA[i] - calc_info_.evalsA[a + noccA]) * t);
            }
        } 
        for (int i = 0; i < noccB; i++) {
            for (int a = 0; a < nvirB; a++) {
                Lbs[vec][i * nvirB + a] = w * exp((calc_info_.evalsB[i] - calc_info_.evalsB[a + noccB]) * t);
            }
        } 
        vec++ ;
    } 

    calc_info_.Lar = Lar;
    calc_info_.Lbs = Lbs;
    calc_info_.ndelta = npoint;
    **/
}
CholeskyDenominator::CholeskyDenominator(shared_ptr<Vector> eps_occ, shared_ptr<Vector> eps_vir, double delta, double dm) : 
    Denominator(eps_occ, eps_vir), delta_(delta), degeneracy_multiplier_(dm)
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
    
    int Ndelta;
    for (Ndelta = 1; Ndelta <= N; Ndelta++) {
        bool OK = true;
        for (int p = Ndelta; p < N; p++) {
            double D = 1.0 / (2.0 * w_ia[p]);
            double Q = 0.0;
            for (int m = 0; m < Ndelta - 1; m++) {
                Q = (w_ia[p] - w_ia[m]) / (w_ia[p] + w_ia[m]);
                D *= Q * Q;
            }
            //printf("  D(%d)^%d = %14.10E\n", p, Ndelta, D);
            if (fabs(D) > delta_) {
                OK = false;
                break;
            }
        }
        if (OK) 
            break;
    }


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

    //fprintf(outfile, "  Energy Denominator Cholesky Decomposition\n");
    //fprintf(outfile, "   Delta = %.3E [H^(-1/2)]\n", delta);
    //fprintf(outfile, "   Number of vectors required = %d [vectors]\n", Ndelta);
    //fprintf(outfile, "   Memory required = %ld [bytes]\n\n", Ndelta*nspan*8L);
    //L->print();

    denominator_ = shared_ptr<Matrix>(new Matrix("Cholesky Delta Tensor", nspan, Ndelta));
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

}

