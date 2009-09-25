/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <ccfiles.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

/* pair_energies(): For RHF references, compute pair energies. Spin-adapt
** pair energies if SPINADAPT_ENERGIES is set to true.
**
** E(IJ) = T2(IJ,AB) * (<ij|ab> - <ij|ba>)
** E(Ij) = T2(Ij,Ab) * <ij|ab>
**
*/

void pair_energies(double** epair_aa, double** epair_ab)
{
  dpdbuf4 tau, D, E;

  if(params.ref == 0) { /** RHF **/

    int i, j, ij;
    int irrep;
    int nocc_act = 0;
    int naa, nab;
    for(irrep=0; irrep<moinfo.nirreps; irrep++)
      nocc_act += moinfo.clsdpi[irrep];
    naa = nocc_act * (nocc_act-1)/2;
    nab = nocc_act * nocc_act;

    /* Compute alpha-alpha pair energies */
    if (naa) {
      double* eaa = init_array(naa);
      
      dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 0, 5, 1, "D <ij|ab>");
      dpd_buf4_init(&tau, CC_TAMPS, 0, 2, 5, 0, 5, 1, "tauIjAb");
      dpd_buf4_init(&E, CC_TMP0, 0, 2, 2, 2, 2, 0, "E <ij|kl>");
      dpd_contract444(&D, &tau, &E, 0, 0, 1.0, 0.0);

      //dpd_buf4_print(&E, outfile, 1);

      /* Extract diagonal elements (i.e. pair energies) and print them out nicely */
      for(irrep=0; irrep<moinfo.nirreps; irrep++) {
        double **block;
        dpdparams4 *Params = E.params;
        int p;
        int np = Params->rowtot[irrep];

        dpd_buf4_mat_irrep_init(&E, irrep);
        dpd_buf4_mat_irrep_rd(&E, irrep);
        block = E.matrix[irrep];

        for(p=0; p<np; p++) {
          int i, j, ij;
          
          i = Params->roworb[irrep][p][0];
          j = Params->roworb[irrep][p][1];
          
          ij = (i > j) ? i*(i-1)/2 + j : j*(j-1)/2 + i;
          eaa[ij] = block[p][p];
        }
	dpd_buf4_mat_irrep_close(&E, irrep);
      }

      *epair_aa = eaa;

      dpd_buf4_close(&tau);
      dpd_buf4_close(&D);
      dpd_buf4_close(&E);
    }

    /* Compute alpha-beta pair energies */
    if (nab) {
      double* eab = init_array(nab);

      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      dpd_buf4_init(&tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
      dpd_buf4_init(&E, CC_TMP0, 0, 0, 0, 0, 0, 0, "E <ij|kl>");
      dpd_contract444(&D, &tau, &E, 0, 0, 1.0, 0.0);

      //dpd_buf4_print(&E, outfile, 1);

      /* Extract diagonal elements (i.e. pair energies) and print them out nicely */
      for(irrep=0; irrep<moinfo.nirreps; irrep++) {
        double **block;
        dpdparams4 *Params = E.params;
        int p;
        int np = Params->rowtot[irrep];

        dpd_buf4_mat_irrep_init(&E, irrep);
        dpd_buf4_mat_irrep_rd(&E, irrep);
        block = E.matrix[irrep];

        for(p=0; p<np; p++) {
          int i, j, ij;

          i = Params->roworb[irrep][p][0];
          j = Params->roworb[irrep][p][1];

          ij = i*nocc_act + j;
          eab[ij] = block[p][p];
        }
	dpd_buf4_mat_irrep_close(&E, irrep);
      }

      *epair_ab = eab;

      dpd_buf4_close(&tau);
      dpd_buf4_close(&D);
      dpd_buf4_close(&E);
    }
  }
  
}

void print_pair_energies(double* emp2_aa, double* emp2_ab, double* ecc_aa, double* ecc_ab)
{
  if(params.ref == 0) { /** RHF **/

    int i, j, ij;
    int irrep;
    int nocc_act = 0;
    int naa, nab;
    for(irrep=0; irrep<moinfo.nirreps; irrep++)
      nocc_act += moinfo.clsdpi[irrep];
    naa = nocc_act * (nocc_act-1)/2;
    nab = nocc_act * nocc_act;

    if (!params.spinadapt_energies) {
      double emp2_aa_tot = 0.0;
      double emp2_ab_tot = 0.0;
      double ecc_aa_tot = 0.0;  
      double ecc_ab_tot = 0.0;  

      fprintf(outfile, "\tAlpha-alpha pair energies\n");
      fprintf(outfile, "\t    i       j         MP2             %s\n",params.wfn);
      fprintf(outfile, "\t  -----   -----   ------------   ------------\n");
      if (naa) {
        ij = 0;
        for(i=0; i<nocc_act; i++)
          for(j=0; j<i; j++,ij++) {
            fprintf(outfile, "\t  %3d     %3d     %12.9lf   %12.9lf\n", i+1, j+1, emp2_aa[ij], ecc_aa[ij]);
            emp2_aa_tot += emp2_aa[ij];
            ecc_aa_tot += ecc_aa[ij];
          }
      }
      fprintf(outfile, "\t  -------------   ------------   ------------\n");
      fprintf(outfile, "\t      Total       %12.9lf   %12.9lf\n\n", emp2_aa_tot, ecc_aa_tot);

      fprintf(outfile, "\tAlpha-beta pair energies\n");
      fprintf(outfile, "\t    i       j         MP2             %s\n",params.wfn);
      fprintf(outfile, "\t  -----   -----   ------------   ------------\n");
      if (nab) {
        ij = 0;
        for(i=0; i<nocc_act; i++)
          for(j=0; j<nocc_act; j++,ij++) {
            fprintf(outfile, "\t  %3d     %3d     %12.9lf   %12.9lf\n", i+1, j+1, emp2_ab[ij], ecc_ab[ij]);
            emp2_ab_tot += emp2_ab[ij];
            ecc_ab_tot += ecc_ab[ij];
          }
      }
      fprintf(outfile, "\t  -------------   ------------   ------------\n");
      fprintf(outfile, "\t      Total       %12.9lf   %12.9lf\n", emp2_ab_tot, ecc_ab_tot);

    }
    else {
      /* Spin-adapt pair energies */
      double emp2_s_tot = 0.0;
      double emp2_t_tot = 0.0;
      double ecc_s_tot = 0.0;
      double ecc_t_tot = 0.0;

      fprintf(outfile, "\tSinglet pair energies\n");
      fprintf(outfile, "\t    i       j         MP2             %s\n",params.wfn);
      fprintf(outfile, "\t  -----   -----   ------------   ------------\n");
      ij = 0;
      for(i=0; i<nocc_act; i++)
        for(j=0; j<=i; j++,ij++) {
          double emp2_s, ecc_s;
          int ij_ab = i*nocc_act + j;
          int ij_aa = i*(i-1)/2 + j;
          double eab, eaa;
          
          /* MP2 */
          eab = emp2_ab[ij_ab];
          if (i != j)
            eaa = emp2_aa[ij_aa];
          else
            eaa = 0.0;
          emp2_s = (i != j ? 2.0 : 1.0) * eab - 0.5 * eaa;

          /* CC */
          eab = ecc_ab[ij_ab];
          if (i != j)
            eaa = ecc_aa[ij_aa];
          else
            eaa = 0.0;
          ecc_s = (i != j ? 2.0 : 1.0) * eab - 0.5 * eaa;
          
          fprintf(outfile, "\t  %3d     %3d     %12.9lf   %12.9lf\n", i+1, j+1, emp2_s, ecc_s);
          emp2_s_tot += emp2_s;
          ecc_s_tot += ecc_s;
        }
      fprintf(outfile, "\t  -------------   ------------   ------------\n");
      fprintf(outfile, "\t      Total       %12.9lf   %12.9lf\n\n", emp2_s_tot, ecc_s_tot);

      fprintf(outfile, "\tTriplet pair energies\n");
      fprintf(outfile, "\t    i       j         MP2             %s\n",params.wfn);
      fprintf(outfile, "\t  -----   -----   ------------   ------------\n");
      if (naa) {
        ij = 0;
        for(i=0; i<nocc_act; i++)
          for(j=0; j<i; j++,ij++) {
            fprintf(outfile, "\t  %3d     %3d     %12.9lf   %12.9lf\n", i+1, j+1, 1.5*emp2_aa[ij], 1.5*ecc_aa[ij]);
            emp2_t_tot += 1.5*emp2_aa[ij];
            ecc_t_tot += 1.5*ecc_aa[ij];
          }
      }
      fprintf(outfile, "\t  -------------   ------------   ------------\n");
      fprintf(outfile, "\t      Total       %12.9lf   %12.9lf\n", emp2_t_tot, ecc_t_tot);

    }      
    
    fprintf(outfile, "\n");
  }
}
}} // namespace psi::ccenergy
