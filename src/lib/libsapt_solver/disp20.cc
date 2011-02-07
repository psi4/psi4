/* This function calculates the Disp20 energy */

#define EXTERN

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <cstring>
#include <iostream>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libutil/quad.h>
#include "structs.h"
#include "sapt0.h"
#include "sapt_dft.h"
#include "sapt2p.h"

#include <vector>
#include <algorithm>
#include <utility>

using namespace std;

namespace psi { namespace sapt {

void SAPT0::disp20()
{
  double energy=0.0;
  long int avail_mem = params_.memory - sizeof(double)*calc_info_.noccB*
    calc_info_.nvirB*(long int) calc_info_.nrio;

  if (params_.print)
    fprintf(outfile,"Begining Disp20 Calculation\n");

  long int temp_size = avail_mem / (sizeof(double) * (2*calc_info_.noccB*
    calc_info_.nvirB + (long int) calc_info_.nrio));

  if (temp_size < 1) {
    fprintf(outfile,"Not enough memory in Disp20\n\n");
    exit(0);
  }

  if (temp_size > calc_info_.noccA*calc_info_.nvirA)
    temp_size = calc_info_.noccA*calc_info_.nvirA;

  double **B_p_AR = block_matrix(temp_size,calc_info_.nrio);
  double **B_p_BS = get_BS_ints(0);
  double **ARBS = block_matrix(temp_size,calc_info_.noccB*calc_info_.nvirB);
  double **tARBS = block_matrix(temp_size,calc_info_.noccB*calc_info_.nvirB);

  int blocks = (calc_info_.noccA*calc_info_.nvirA)/temp_size;
  if ((calc_info_.noccA*calc_info_.nvirA)%temp_size) blocks++;

  if (params_.print) {
    fprintf(outfile,"T2 ARBS Amplitudes formed in %d chunks\n\n",blocks);
    fflush(outfile);
  }

  psio_address next_PSIF_DF_AR = PSIO_ZERO;
  psio_address next_PSIF_tARBS = PSIO_ZERO;

  for (int t_ar=0; t_ar<blocks; t_ar++) {
    int ar_start = temp_size*t_ar;
    int ar_stop = temp_size*(t_ar+1);
    if (ar_stop > calc_info_.noccA*calc_info_.nvirA)
      ar_stop = calc_info_.noccA*calc_info_.nvirA;

    psio_->read(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *) &(B_p_AR[0][0]),
      sizeof(double)*(ar_stop-ar_start)*(ULI) calc_info_.nrio,next_PSIF_DF_AR,
      &next_PSIF_DF_AR);

    C_DGEMM('N','T',ar_stop-ar_start,calc_info_.noccB*calc_info_.nvirB,
      calc_info_.nrio,1.0,&(B_p_AR[0][0]),calc_info_.nrio,&(B_p_BS[0][0]),
      calc_info_.nrio,0.0,&(ARBS[0][0]),calc_info_.noccB*calc_info_.nvirB);

    //shared_ptr<Matrix> Delta(new Matrix(calc_info_.noccA*calc_info_.nvirA, calc_info_.noccB*calc_info_.nvirB));
    //double** Delt = Delta->get_pointer();
    
    #pragma omp for schedule(static)
      for (int ar=ar_start; ar<ar_stop; ar++) {
        int a = ar/calc_info_.nvirA;
        int r = ar%calc_info_.nvirA + calc_info_.noccA;
        for (int b=0,bs=0; b<calc_info_.noccB; b++) {
          for (int s=calc_info_.noccB; s<calc_info_.nmo; s++,bs++) {
            double denom = calc_info_.evalsA[a]+calc_info_.evalsB[b]-
                           calc_info_.evalsA[r]-calc_info_.evalsB[s];
            tARBS[ar-ar_start][bs] = ARBS[ar-ar_start][bs]/denom;
            //Delt[a*calc_info_.nvirA + r - calc_info_.noccA][b*calc_info_.nvirB + s - calc_info_.noccB] = 1.0 / denom;
        }}
      }

    //Delta->print();

    energy += C_DDOT((ar_stop-ar_start)*calc_info_.noccB*calc_info_.nvirB,
              &(ARBS[0][0]),1,&(tARBS[0][0]),1);

    psio_->write(PSIF_SAPT_AMPS,"T ARBS Amplitudes",(char *) &(tARBS[0][0]),
      sizeof(double)*(ar_stop-ar_start)*calc_info_.noccB*
      (ULI) calc_info_.nvirB,next_PSIF_tARBS,&next_PSIF_tARBS);
  }

  free_block(ARBS);
  free_block(tARBS);
  free_block(B_p_AR);
  free_block(B_p_BS);

  results_.disp20 = energy*4.0;

  if (params_.print) {
    fprintf(outfile,"Dispersion Energy     = %18.12lf mH\n\n",energy*4000.0);
    fflush(outfile);
  }
}

void SAPT_DFT::disp20()
{
  double disp=0.0;
  double exdisp=0.0;

  double **B_p_AR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    calc_info_.noccA*calc_info_.nvirA);
  double **B_p_BS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    calc_info_.noccB*calc_info_.nvirB);
  double **ARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccB*calc_info_.nvirB);
  double **tARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccB*calc_info_.nvirB);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(B_p_BS[0][0]),calc_info_.nrio,0.0,&(ARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  for (int a=0, ar=0; a < calc_info_.noccA; a++) {
  for (int r=0; r < calc_info_.nvirA; r++, ar++) {
    for (int b=0, bs=0; b < calc_info_.noccB; b++) {
    for (int s=0; s < calc_info_.nvirB; s++, bs++) {
      double denom = calc_info_.evalsA[a]+calc_info_.evalsB[b]-
        calc_info_.evalsA[r+calc_info_.noccA]-
        calc_info_.evalsB[s+calc_info_.noccB];
      tARBS[ar][bs] = ARBS[ar][bs]/denom;
    }}
  }}
  disp = 4.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.noccB*
    calc_info_.nvirB,&(tARBS[0][0]),1,&(ARBS[0][0]),1);

  psio_->read_entry(PSIF_SAPT_AMPS,"Exch-Disp V_ARBS",(char *)
    &(ARBS[0][0]),calc_info_.noccA*calc_info_.nvirA*calc_info_.noccB*
    calc_info_.nvirB*(ULI) sizeof(double));

  exdisp = C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.noccB*
    calc_info_.nvirB,&(tARBS[0][0]),1,&(ARBS[0][0]),1);

  free_block(ARBS);
  free_block(tARBS);

  if (params_.print) {
    fprintf(outfile,"Disp20                = %18.12lf H\n\n",disp);
    fprintf(outfile,"Exch-Disp20           = %18.12lf H\n\n",exdisp);
    fflush(outfile);
  }

  results_.disp20 = disp;
  results_.exch_disp20 = exdisp;
}

void SAPT2p::disp20()
{
  double energy=0.0;

  if (params_.print)
    fprintf(outfile,"Begining Disp20 Calculation\n\n");

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"T(BS) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);
  double **B_p_AR = get_AR_ints(1);

  energy = C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nrio,
    B_p_AR[0],1,T_p_AR[0],1);

  free_block(B_p_AR);
  free_block(T_p_AR);

  results_.disp20 = energy*4.0;

  if (params_.print) {
    fprintf(outfile,"Dispersion Energy     = %18.12lf mH\n\n",energy*4000.0);
    fflush(outfile);
  }
}
bool SAPT0::denominator_criteria(std::pair<int, double> a, std::pair<int, double> b)
{
    return a.second < b.second;
}
void SAPT0::cholesky_denominator(double delta)
{
    int noccA = calc_info_.noccA;    
    int noccB = calc_info_.noccB;    
    int nvirA = calc_info_.nvirA;    
    int nvirB = calc_info_.nvirB;
    int nspan = noccA*nvirA + noccB*nvirB; 

    std::vector<std::pair<int, double> > eps(nspan);    
    std::vector<int> degenerate_index(nspan);

    for (int i = 0; i < noccA; i++) {
        for (int a = 0; a < nvirA; a++) {
            // minus means monomer A
            // and we're using FORTRAN indexing. Woot!
            eps[i*nvirA + a] = make_pair(-(i*nvirA + a + 1), calc_info_.evalsA[a + noccA] - calc_info_.evalsA[i]);
        }
    }

    for (int j = 0; j < noccB; j++) {
        for (int b = 0; b < nvirB; b++) {
            // Remember, we're using FORTRAN indexing
            eps[j*nvirB + b + noccA*nvirA] = make_pair(j*nvirB + b + 1, calc_info_.evalsB[b + noccB] - calc_info_.evalsB[j]);
        }
    }

    std::sort(eps.begin(), eps.end(), &SAPT0::denominator_criteria);

    double* w_ia = new double[nspan];
    w_ia[0] = eps[0].second;
    degenerate_index[0] = 0;
    
    int N = 1;
    for (int k = 1; k < nspan; k++) {
        if (fabs(eps[k].second - w_ia[N-1]) > 1.0E-12) {
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
            if (fabs(D) > delta) {
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

    fprintf(outfile, "  Energy Denominator Cholesky Decomposition\n");
    fprintf(outfile, "   Delta = %.3E [H^(-1/2)]\n", delta);
    fprintf(outfile, "   Number of vectors required = %d [vectors]\n", Ndelta);
    fprintf(outfile, "   Memory required = %ld [bytes]\n\n", Ndelta*nspan*8L);
    //L->print();

    double** Lar = block_matrix(Ndelta, noccA*nvirA); 
    double** Lbs = block_matrix(Ndelta, noccB*nvirB); 

    for (int k = 0; k < nspan; k++) {
        int ind = eps[k].first;
        if (ind < 0) {
            ind = -ind - 1;
            C_DCOPY(Ndelta, Lp[degenerate_index[k]], 1, &Lar[0][ind], noccA*nvirA);
        } else {
            ind = ind - 1;
            C_DCOPY(Ndelta, Lp[degenerate_index[k]], 1, &Lbs[0][ind], noccB*nvirB);
        }
    }

    //fprintf(outfile, " Lar\n");
    //print_mat(Lar, Ndelta, noccA*nvirA, outfile);
    //fprintf(outfile, " Lbs\n");
    //print_mat(Lbs, Ndelta, noccA*nvirA, outfile);

    //shared_ptr<Matrix> Delta(new Matrix(noccA*nvirA, noccB*nvirB));
    //double** Delt = Delta->get_pointer();

    //C_DGEMM('T','N', noccA*nvirA, noccB*nvirB, Ndelta, 1.0, Lar[0], noccA*nvirA, Lbs[0], noccB*nvirB, 0.0, Delt[0], noccB*nvirB); 
    //Delta->print();
   
    calc_info_.Lar = Lar;
    calc_info_.Lbs = Lbs;
    calc_info_.ndelta = Ndelta;
 
}
void SAPT0::laplace_denominator(int npoint, double center)
{
    int noccA = calc_info_.noccA;    
    int noccB = calc_info_.noccB;    
    int nvirA = calc_info_.nvirA;    
    int nvirB = calc_info_.nvirB;
    int nspan = noccA*nvirA + noccB*nvirB; 

    fprintf(outfile, "  Energy Denominator Laplace Decomposition\n");
    fprintf(outfile, "   Number of vectors required = %d [vectors]\n", npoint);
    fprintf(outfile, "   Memory required = %ld [bytes]\n\n", npoint*nspan*8L);
    fprintf(outfile, "\n");
    
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
    
}

}}
