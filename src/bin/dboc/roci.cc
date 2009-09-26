/*! \file
 \ingroup DBOC
 \brief Enter brief description of file here 
 */
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libqt/slaterdset.h>
#include <psifiles.h>
#include "moinfo.h"
#include "float.h"
#include "linalg.h"
#include "mo_overlap.h"
#include "hfwfn.h"
#include <psi4-dec.h>

// Set to 1 to try using unit MO overlap between + and - displacements
#define TRY_UNIT_OVERLAP 0

using namespace std;
namespace psi {
  namespace DBOC {
    
    // Wrap a,b indices into one composite index assuming S2 symmetry
#define INDEX2(a,b) ((a) > (b)) ? ( (((a)*(a+1)) >> 1) + (b) ) : ( (((b)*(b+1)) >> 1) + (a) )
    // Wrap a>=b indices into one composite index assuming S2 symmetry
#define INDEX2_ORD(a,b) ( (((a)*(a+1))/2) + (b) )
    // Wrap a>=b>=c indices into one composite index assuming S3 symmetry
#define INDEX3_ORD(a,b,c) ( ((a)*(((a)+4)*((a)-1)+6)/6) + (((b)*(b+1))/2) + (c) )
    
    extern char *CI_Vector_Labels[MAX_NUM_DISP];
    extern HFWavefunction* HFVectors[MAX_NUM_DISP];
    extern void done(const char *);
    extern void mo_maps(short int**, short int**);

    double eval_roci_derwfn_overlap(DisplacementIndex LDisp, DisplacementIndex RDisp)
    {
      // Read in CI vectors
      SlaterDetVector *vecm, *vecp;
      slaterdetvector_read(PSIF_CIVECT,CI_Vector_Labels[RDisp],&vecm);
      slaterdetvector_read(PSIF_CIVECT,CI_Vector_Labels[LDisp],&vecp);

      int nfzc = vecm->sdset->alphastrings->nfzc;
#if USE_MOINFO
      extern MOInfo_t MOInfo;
      int nalpha = MOInfo.nalpha;
      int nbeta = MOInfo.nbeta;
#else
      int nalpha = HFVectors[LDisp]->nalpha();
      int nbeta = HFVectors[LDisp]->nbeta();
#endif
      int nact_a = nalpha - nfzc;
      int nact_b = nbeta - nfzc;

#if !TRY_UNIT_OVERLAP
      FLOAT **CSC_full = eval_S_alpha(LDisp,RDisp);
#else
      const int num_mo = HFVectors[LDisp]->num_mo();
      FLOAT **CSC_full = create_matrix(num_mo, num_mo);
      for(int i=0; i<num_mo; ++i)
        CSC_full[i][i] = 1.0;
#endif
      FLOAT **CSC_a = create_matrix(nalpha,nalpha);
      FLOAT **CSC_b = create_matrix(nbeta,nbeta);
      int *tmpintvec = new int[nalpha];

      // Compute overlap between strings for alpha spin case
      StringSet *ssetm;
      ssetm = vecm->sdset->alphastrings;
      short int* fzc_occ = ssetm->fzc_occ;
      int nstr_a = ssetm->size;
      FLOAT **S_a = create_matrix(nstr_a,nstr_a);
      // Assume the order of strings is the same for - and + displacements
      for(int jp=0; jp<nstr_a; jp++) {
        String *str_j = &ssetm->strings[jp];
        for(int im=0; im<nstr_a; im++) {
          String *str_i = &ssetm->strings[im];

          for(int j=0;j<nact_a;j++)
            for(int i=0;i<nact_a;i++)
              CSC_a[j+nfzc][i+nfzc] = CSC_full[str_j->occ[j]][str_i->occ[i]];

          // frozen orbitals need to be mapped to pitzer order manually
          for(int i=0; i<nfzc; i++) {
            int ii = fzc_occ[i];
            for(int j=0; j<nact_a; j++)
              CSC_a[j+nfzc][i] = CSC_full[str_j->occ[j]][ii];
          }

          for(int j=0; j<nfzc; j++) {
            int jj = fzc_occ[j];
            for(int i=0; i<nact_a; i++)
              CSC_a[j][i+nfzc] = CSC_full[jj][str_i->occ[i]];
          }

          for(int i=0;i<nfzc;i++) {
            int ii = fzc_occ[i];
            for(int j=0;j<nfzc;j++)
              CSC_a[i][j] = CSC_full[ii][fzc_occ[j]];
          }

          FLOAT sign;
          lu_decom(CSC_a, nalpha, tmpintvec, &sign);
          FLOAT deter1 = 1.0;
          for(int i=0;i<nalpha;i++)
          deter1 *= CSC_a[i][i];

          S_a[jp][im] = sign*deter1;
        }
      }

      // Compute overlap between strings for beta spin case
      ssetm = vecm->sdset->betastrings;
      int nstr_b = ssetm->size;
      FLOAT **S_b = create_matrix(nstr_b,nstr_b);
      // Assume the order of strings is the same for - and + displacements
      for(int jp=0; jp<nstr_b; jp++) {
        String *str_j = &ssetm->strings[jp];
        for(int im=0; im<nstr_b; im++) {
          String *str_i = &ssetm->strings[im];

          for(int j=0;j<nact_b;j++)
            for(int i=0;i<nact_b;i++)
              CSC_b[j+nfzc][i+nfzc] = CSC_full[str_j->occ[j]][str_i->occ[i]];

          // frozen orbitals need to be mapped to pitzer order manually
          for(int i=0; i<nfzc; i++) {
            int ii = fzc_occ[i];
            for(int j=0; j<nact_b; j++)
              CSC_b[j+nfzc][i] = CSC_full[str_j->occ[j]][ii];
          }

          for(int j=0; j<nfzc; j++) {
            int jj = fzc_occ[j];
            for(int i=0; i<nact_b; i++)
              CSC_b[j][i+nfzc] = CSC_full[jj][str_i->occ[i]];
          }

          for(int i=0;i<nfzc;i++) {
            int ii = fzc_occ[i];
            for(int j=0;j<nfzc;j++)
              CSC_b[i][j] = CSC_full[ii][fzc_occ[j]];
          }

          FLOAT sign;
          lu_decom(CSC_b, nbeta, tmpintvec, &sign);
          FLOAT deter1 = 1.0;
          for(int i=0;i<nbeta;i++)
          deter1 *= CSC_b[i][i];

          S_b[jp][im] = sign*deter1;
        }
      }

      // Evaluate total overlap in the highest available precision
      int ndets = vecm->size;
      FLOAT S_tot = 0.0;
      for(int I=0; I<ndets; I++) {
        SlaterDet *detI = vecm->sdset->dets + I;
        int Istra = detI->alphastring;
        int Istrb = detI->betastring;
        FLOAT cI = vecm->coeffs[I];

        for(int J=0; J<ndets; J++) {
          SlaterDet *detJ = vecp->sdset->dets + J;
          int Jstra = detJ->alphastring;
          int Jstrb = detJ->betastring;
          FLOAT cJ = vecp->coeffs[J];

          FLOAT S = S_a[Jstra][Istra] * S_b[Jstrb][Istrb];

          FLOAT contrib = cI * S * cJ;
          S_tot += cI * S * cJ;
          /* fprintf(outfile,"  %3d %3d %+15.10Le", I, J, cI);
           fprintf(outfile," %+15.10Le", cJ);
           fprintf(outfile," %+25.15Le", S);
           fprintf(outfile," %+25.15Le", contrib);
           fprintf(outfile," %+25.15Le\n", S_tot); */

        }
      }

      // Cleanup
      slaterdetvector_delete_full(vecm);
      slaterdetvector_delete_full(vecp);
      delete[] tmpintvec;
      delete_matrix(CSC_a);
      delete_matrix(CSC_full);
      delete_matrix(CSC_b);
      delete_matrix(S_a,nstr_a,nstr_a);
      delete_matrix(S_b,nstr_b,nstr_b);
      double S_tot_double = (double) S_tot;
      return fabs(S_tot_double);
    }

  }} // namespace psi::DBOC
