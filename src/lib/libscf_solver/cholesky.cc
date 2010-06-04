#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include "hf.h"
#include "rhf.h"
#include "uhf.h"
#include "rohf.h"

#include <libmints/mints.h>

using namespace std;
using namespace psi;

namespace psi { namespace scf {

void HF::form_CD()
{
    fprintf(outfile, "\n  Computing Integrals using Cholesky Decomposition\n");
    //TODO: Add support for molecular symmetry
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1 for now.\n"); fflush(outfile);
        abort();
    } 

  double tol = options_.get_double("CHOLESKY_CUTOFF");
  double schwartz_tol = options_.get_double("SCHWARZ_CUTOFF");

  // Create integral factory
  IntegralFactory cdfactory(basisset_, basisset_, basisset_, basisset_);
  TwoBodyInt* eri = cdfactory.eri();

  int ntri = basisset_->nbf()*(basisset_->nbf()+1)/2;
  int shelltri = basisset_->nshell()*(basisset_->nshell()+1)/2;

  int *ij2i = init_int_array(ntri);
  int *ij2j = init_int_array(ntri);

  for (int i=0,ij=0;i<basisset_->nbf();i++) {
    for (int j=0;j<=i;j++,ij++) {
      ij2i[ij] = i;
      ij2j[ij] = j;
    }
  }

  timer_on("Form Schwartz");

  int screened=0;
  double X_max=0.0;
  const double *buffer = eri->buffer();
  double *Schwartz = init_array(shelltri);
  double *diag = init_array(ntri);
  int maxPshell=0;
  int *reorder = init_int_array(shelltri);
  int *snuc = init_int_array(basisset_->nshell());

  snuc = chkpt_->rd_snuc();

  for(int P=0,PQ=0;P<basisset_->nshell();P++) { // Get diagonal 2e terms
    int numw = basisset_->shell(P)->nfunction();

    if (numw > maxPshell) maxPshell = numw;

    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      reorder[PQ] = PQ;

      eri->compute_shell(P, Q, P, Q);

      for (int w=0;w<numw;w++) {
        int absw = basisset_->shell(P)->function_index()+w;
        for(int x=0;x<numx;x++) { 
          int absx = basisset_->shell(Q)->function_index()+x;
          int index = ( ( (w*numx + x) * numw + w) * numx + x);
          double tei = buffer[index];
          diag[INDEX2(absw,absx)] = tei;
          if (tei > X_max) X_max = tei;
          if(fabs(tei) > Schwartz[PQ]) {
            Schwartz[PQ] = fabs(tei);
          }
        }
      }
    }
  }

  sort_cholesky(Schwartz, reorder, shelltri);

  timer_off("Form Schwartz");
  timer_on("Full CD");

  int numL = 0;

  for (int i=0; i<basisset_->nbf()*(basisset_->nbf()+1)/2; i++) {
    if (sqrt(diag[i]*X_max) > tol) { 
      numL++;                  // How many Cholesky vectors are neccesary
    }
  }

  double **L = (double **) malloc(numL * sizeof(double**));

  for (int i=0; i<numL; i++) {
    L[i] = NULL;
  }

  double **temp = block_matrix(maxPshell*maxPshell,ntri);

  int pos=0;

  // Loop over all possible Cholesky vectors
  for (int PQshell=0; PQshell<shelltri; PQshell++) {
    int P = ij2i[reorder[PQshell]];
    int Q = ij2j[reorder[PQshell]];
    // Determine absolute start and end to temp
    int Pstart =  basisset_->shell(P)->function_index();
    int Qstart =  basisset_->shell(Q)->function_index();
    int numw = basisset_->shell(P)->nfunction();
    int numx = basisset_->shell(Q)->nfunction();

    // Only continue if the shells are centered on the same nuclei
    if (scf_type_ == "1C_CD" &&  snuc[P] != snuc[Q]) continue;

    // Only continue if some diagonal element is greater than tol
    if (tol > Schwartz[PQshell]) break;

    double dmax = 0.0;

    for(int w=0;w<numw;w++) {
      for(int x=0;x<numx;x++) {
        int ij = INDEX2(Pstart+w,Qstart+x);
        int k = ij2i[ij];
        int l = ij2j[ij];
        int kl = INDEX2(k,l);
        double tval = diag[kl];
        if (Pstart+w>=Qstart+x && tval > tol) {

        for(int i=0; i<pos; i++) {
          tval -= L[i][kl] * L[i][kl];
        }

        if (dmax < tval) dmax = tval;
    }}}

    if (tol > dmax) continue;

    bzero(&(temp[0][0]), maxPshell*maxPshell*ntri*sizeof(double));

timer_on("Compute Integrals");

      // Compute all possible integrals for that PQ shell
      for (int RSshell=0; RSshell<shelltri; RSshell++) {

        if (schwartz_tol > sqrt(Schwartz[RSshell]*Schwartz[PQshell])) {
          screened++;
          continue;
          }

        int R = ij2i[reorder[RSshell]];
        int S = ij2j[reorder[RSshell]];
        int numy = basisset_->shell(R)->nfunction();
        int numz = basisset_->shell(S)->nfunction();

        eri->compute_shell(P, Q, R, S);

        for(int w=0,wx=0,index=0;w<numw;w++) {
          for(int x=0;x<numx;x++,wx++) {
            for(int y=0;y<numy;y++) {
              for(int z=0;z<numz;z++,++index) {
                int k = basisset_->shell(R)->function_index()+y;
                int l = basisset_->shell(S)->function_index()+z;
                double tei = buffer[index];
                temp[wx][INDEX2(k,l)] = tei;
              }
            }
          }
        }
      }

timer_off("Compute Integrals");
timer_on("Cholesky Decomp");

      // Loop over unique elements of temp
      for(int w=0,wx=0;w<numw;w++) {
        for(int x=0;x<numx;x++,wx++) {
          int ij = INDEX2(Pstart+w,Qstart+x);
          int k = ij2i[ij];
          int l = ij2j[ij];
          int kl = INDEX2(k,l);
          if (Pstart+w>=Qstart+x && temp[wx][kl] > tol) {

            if (L[pos] == NULL) 
              L[pos] = init_array(ntri);

            double tval = temp[wx][kl];
            for(int i=0; i<pos; i++) {
              tval -= L[i][kl] * L[i][kl];
            }

            if (tval > 0.0) {
              L[pos][kl] = sqrt(tval);
            }
            else {
              L[pos][kl] = 0.0;
            }

            if (L[pos][kl]*L[pos][kl] > tol) {

              for(int n=0; n<pos; n++) {
                C_DAXPY(ntri,-L[n][kl],&(L[n][0]),1,&(temp[wx][0]),1);
              }

              tval = L[pos][kl];
              C_DSCAL(ntri,1.0/tval,&(temp[wx][0]),1);
              C_DCOPY(ntri,&(temp[wx][0]),1,&(L[pos][0]),1);
              L[pos][kl] = tval;
              pos++;
              }
            }
          }
        }
timer_off("Cholesky Decomp");
      }

  fprintf(outfile,"  %d shell triplets screened via Schwartz inequality\n\n",
    screened);
  fflush(outfile);

  timer_off("Full CD");

  free_block(temp);

  ntri_naive_ = ntri;
  ntri_ = ntri;
  ri_pair_mu_ = init_int_array(ntri_);  
  ri_pair_nu_ = init_int_array(ntri_);  
  for (int i = 0, ij = 0; i < basisset_->nbf(); i++)
    for (int j = 0; j<=i; j++, ij++) { 
        ri_pair_mu_[ij] = i;
        ri_pair_nu_[ij] = j;
    }
  
  ri_nbf_ = pos;
  df_storage_ = full;
	
  psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
  psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
  
  for (int Q = 0; Q < ri_nbf_; Q++) {
    psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(L[Q][0]),sizeof(double)*ntri_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ); 
    free(L[Q]);  
  }
  free(L);

  psio_->close(PSIF_DFSCF_BJ,1); 

  B_ia_P_ = block_matrix(ri_nbf_,ntri_);

  if (print_> 3) 
    fprintf(outfile,"  Cholesky Tensor: ntri_ = %d, ri_nbf_ = %d:\n",ntri_,ri_nbf_);
    

  if (print_ > 7) {
    print_mat(B_ia_P_,ri_nbf_,ntri_,outfile);
  }

  psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
  next_PSIF_DFSCF_BJ = PSIO_ZERO;
  
  for (int Q = 0; Q < ri_nbf_; Q++) {
    psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(B_ia_P_[Q][0]),sizeof(double)*ntri_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);   
  }
  
  psio_->close(PSIF_DFSCF_BJ,0); 

}
void HF::sort_cholesky(double *vec, int *reorder, int length)
{
  int i, j;
  double dtemp;
  int itemp;
  for (i=0; i<length; i++)
  {
    for(j = (i+1); j < length; j++)
    {
      if (vec[i] < vec[j])
      {
        dtemp = vec[i];
        vec[i] = vec[j];
        vec[j] = dtemp;

        itemp = reorder[i];
        reorder[i] = reorder[j];
        reorder[j] = itemp;
      }
    }
  }
  return;
}
}}
