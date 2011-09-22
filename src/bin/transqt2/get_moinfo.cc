/*! \file
    \ingroup TRANSQT2
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

#include <libmints/wavefunction.h>

/* get_moinfo(): Routine to obtain basic orbital information from
** chkpt and compute the associated lookup arrays.
**
** Notes on some of the quantities computed here:
**
** pitz2corr_one/pitz2corr_two: These are orbital-index translation
** arrays for the one-electron integral/frozen-core operator and
** two-electron integral transformations, respectively.  Two different
** arrays are necessary because the one-electron integral/frozen-core
** operator transformation is carried out in the *full* MO space
** (because of the default behavior of the iwl_rdone code is to
** "filter" frozen orbitals) while the two-electron integral
** transformation may be carried out only in the "active" space (if
** integrals involving core orbitals are not needed, for example).
** Furthermore, the definitions of these arrays vary depending on the
** code making use of the transformed integrals.  The CI/MCSCF codes
** require a RAS-type ordering, the CC codes require QT ordering (for
** now), and the SCF codes (for MO-basis CPHF and orbital stability
** calculations) require the original Pitzer order (i.e., no
** re-ordering).
**
** C_offset: This array indicates the number of columns (MOs) of the
** transformation matrix that should be skipped ("offset") in the
** two-electron integrals transformation for the given method.
** Although the full transformation matrix (including frozen and
** restricted orbitals) is used for all transforms, we may want to
** skip frozen-core orbitals, for example.  The C_offset array (which
** must match up in definition with the actpi/actsym arrays) is used
** in the DGEMM calls to skip the unwanted columns/MOs.
**
** core: The is the set of occupied orbitals used in the construction
** of the frozen-core operator in main().  For CI/MCSCF calculations,
** this should include both the frozen-core and the "restricted core"
** defined by the ras_set2() function.  For CC calculations, this is
** only the frozen-core.
**
** TDC, 2/2008
*/

namespace psi {
  namespace transqt2 {

    void semicanonical_fock(void);

    void get_moinfo(Options& options)
    {
      int i, j, k, h, p, q, count;
      double escf;
      int *offset, this_offset;
      int *rstr_docc, *rstr_uocc, **ras_opi;
      double ***C, ***C_a, ***C_b;
      double **C_full, **C_full_a, **C_full_b;
      int *reorder, *reorder_A, *reorder_B;
      double **TMP;

      chkpt_init(PSIO_OPEN_OLD);
      moinfo.nirreps = chkpt_rd_nirreps();
      moinfo.nmo = chkpt_rd_nmo();
      moinfo.nso = chkpt_rd_nso();
      moinfo.nao = chkpt_rd_nao();
      moinfo.labels = chkpt_rd_irr_labs();
      moinfo.enuc = chkpt_rd_enuc();
      escf = chkpt_rd_escf();
      moinfo.sopi = chkpt_rd_sopi();
      moinfo.mopi = chkpt_rd_orbspi();
      moinfo.clsdpi = chkpt_rd_clsdpi();
      moinfo.openpi = chkpt_rd_openpi();
      moinfo.frdocc = Process::environment.reference_wavefunction()->frzcpi();
      moinfo.fruocc = Process::environment.reference_wavefunction()->frzvpi();
      chkpt_close();

      moinfo.nfzc = moinfo.nfzv = 0;
      for(i=0; i < moinfo.nirreps; i++) {
    moinfo.nfzc += moinfo.frdocc[i];
    moinfo.nfzv += moinfo.fruocc[i];
      }

      moinfo.uoccpi = init_int_array(moinfo.nirreps);
      for(i=0; i < moinfo.nirreps; i++)
    moinfo.uoccpi[i] = moinfo.mopi[i] - moinfo.clsdpi[i] - moinfo.openpi[i];

      if(params.semicanonical) semicanonical_fock();

      /* SO symmetry array */
      moinfo.sosym = init_int_array(moinfo.nso);
      for(h=0,count=0; h < moinfo.nirreps; h++)
    for(i=0; i < moinfo.sopi[h]; i++,count++)
      moinfo.sosym[count] = h;

      /* AO/MO transformation matrix */
      if(params.ref == 0 || params.ref == 1) { /* RHF/ROHF */
    C = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    chkpt_init(PSIO_OPEN_OLD);
    for(h=0; h < moinfo.nirreps; h++) {
      C[h] = chkpt_rd_scf_irrep(h);
      if(params.print_lvl > 2) {
        fprintf(outfile, "\n\tMOs for irrep %d:\n",h);
        mat_print(C[h], moinfo.sopi[h], moinfo.mopi[h], outfile);
      }
    }
    C_full = chkpt_rd_scf();
    chkpt_close();
    moinfo.C = C;
    moinfo.C_full = C_full;
      }
      else if(params.ref == 2) { /* UHF */
    C_a = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    C_b = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    chkpt_init(PSIO_OPEN_OLD);
    for(h=0; h < moinfo.nirreps; h++) {
      C_a[h] = chkpt_rd_alpha_scf_irrep(h);
      if(params.print_lvl > 2) {
        fprintf(outfile, "\n\tAlpha MOs for irrep %d:\n",h);
        mat_print(C_a[h], moinfo.sopi[h], moinfo.mopi[h], outfile);
      }
      C_b[h] = chkpt_rd_beta_scf_irrep(h);
      if(params.print_lvl > 2) {
        fprintf(outfile, "\n\tBeta MOs for irrep %d:\n",h);
        mat_print(C_b[h], moinfo.sopi[h], moinfo.mopi[h], outfile);
      }
    }
    C_full_a = chkpt_rd_alpha_scf();
    C_full_b = chkpt_rd_beta_scf();
    chkpt_close();
    moinfo.C_a = C_a;
    moinfo.C_b = C_b;
    moinfo.C_full_a = C_full_a;
    moinfo.C_full_b = C_full_b;
      }

      if(ci_wfn(params.wfn)) {

    /** Compute CI spatial-orbial reordering array(s) for the one-electron transformation **/
    rstr_docc = init_int_array(moinfo.nirreps);
    rstr_uocc = init_int_array(moinfo.nirreps);
    ras_opi = init_int_matrix(4,moinfo.nirreps);

    // below, we need frdocc and fruocc as they appear in input, i.e., they
    // should not be zeroed out even if the frozen orbitals are to be
    // transformed
    moinfo.pitz2corr_one = init_int_array(moinfo.nmo);
    if (!ras_set2(moinfo.nirreps, moinfo.nmo, 1, 1,
              moinfo.mopi, moinfo.clsdpi, moinfo.openpi,
              moinfo.frdocc, moinfo.fruocc,
              rstr_docc, rstr_uocc,
              ras_opi, moinfo.pitz2corr_one, 1, 0, options))
      {
        throw PsiException("Error in ras_set(). Aborting.", __FILE__, __LINE__);
      }

    /* "core" array needed for frozen-core operator */
    /* core for CI wfns is frozen-docc plus restricted-docc */
    moinfo.core = init_int_array(moinfo.nirreps);
    moinfo.ncore = 0;
    for(h=0; h < moinfo.nirreps; h++) {
      moinfo.core[h] = moinfo.frdocc[h] + rstr_docc[h];
      moinfo.ncore += moinfo.core[h];
    }

    /* Definition of nactive, actpi, and actsym varies depending
       on the type of CI calculation */
    if((params.wfn == "MCSCF") || (params.wfn == "CASSCF") ||
       (params.wfn == "RASSCF") || (params.wfn == "DETCAS") ||
       params.dertype == 1) {

      /* Include all orbtials; partial transforms to come later */
      moinfo.nactive = moinfo.nmo;

      moinfo.actpi = init_int_array(moinfo.nirreps);
      for(h=0; h < moinfo.nirreps; h++)
        moinfo.actpi[h] = moinfo.mopi[h];

      moinfo.actsym = init_int_array(moinfo.nactive);
      for(h=0,count=0; h < moinfo.nirreps; h++)
        for(i=0; i < moinfo.actpi[h]; i++,count++)
          moinfo.actsym[count] = h;

      /** Compute CI spatial-orbial reordering array(s) for the two-electron transformation **/

      // Now we need to translate the full Pitzer -> correlated mapping
      // array to one that involves only the active orbitals
      moinfo.pitz2corr_two = init_int_array(moinfo.nactive);
      for (h=0,j=0,k=0; h<moinfo.nirreps; h++) {
        for (i=0; i<moinfo.mopi[h]; i++,j++) {
          // j is now an absolute Pitzer MO index
          // k is an index for Pitzer orbitals NOT including frozen ones
          moinfo.pitz2corr_two[k++] = moinfo.pitz2corr_one[j];
        }
      }

      /* zero offsets for MCSCF wfns and derivs b/c all orbitals needed */
      moinfo.C_offset = init_int_array(moinfo.nirreps);
    }
    else {
      /* Otherwise leave out fzc/fzv */
      moinfo.nactive = moinfo.nmo - moinfo.nfzc - moinfo.nfzv;

      moinfo.actpi = init_int_array(moinfo.nirreps);
      for(h=0; h < moinfo.nirreps; h++)
        moinfo.actpi[h] = moinfo.mopi[h] - moinfo.frdocc[h] - moinfo.fruocc[h];

      moinfo.actsym = init_int_array(moinfo.nactive);
      for(h=0,count=0; h < moinfo.nirreps; h++)
        for(i=0; i < moinfo.actpi[h]; i++,count++)
          moinfo.actsym[count] = h;

      /** Compute CI spatial-orbial reordering array(s) for the two-electron transformation **/

      // Now we need to translate the full Pitzer -> correlated mapping
      // array to one that involves only the active orbitals
      moinfo.pitz2corr_two = init_int_array(moinfo.nactive);
      for (h=0,j=0,k=0; h<moinfo.nirreps; h++) {
        for (i=0; i<moinfo.mopi[h]; i++,j++) {
          // j is now an absolute Pitzer MO index
          // k is an index for Pitzer orbitals NOT including frozen ones
          if (i < moinfo.frdocc[h] || i >= (moinfo.mopi[h]-moinfo.fruocc[h]))
        continue;
          moinfo.pitz2corr_two[k++] = moinfo.pitz2corr_one[j] - moinfo.nfzc;
        }
      }

      /* frozen-core offsets for truncated-CI energies */
      moinfo.C_offset = init_int_array(moinfo.nirreps);
      for(h=0; h < moinfo.nirreps; h++)
        moinfo.C_offset[h] = moinfo.frdocc[h];
    }

    free_int_matrix(ras_opi);
    free(rstr_docc);
    free(rstr_uocc);
      }
      else if(cc_wfn(params.wfn) || (params.wfn == "MP2")) {

    /* Leave out fzc/fzv */
    moinfo.nactive = moinfo.nmo - moinfo.nfzc - moinfo.nfzv;

    moinfo.actpi = init_int_array(moinfo.nirreps);
    for(h=0; h < moinfo.nirreps; h++)
      moinfo.actpi[h] = moinfo.mopi[h] - moinfo.frdocc[h] - moinfo.fruocc[h];

    moinfo.actsym = init_int_array(moinfo.nactive);
    for(h=0,count=0; h < moinfo.nirreps; h++)
      for(i=0; i < moinfo.actpi[h]; i++,count++)
        moinfo.actsym[count] = h;

    /* "core" orbitals needed for frozen-core operator */
    /* core for CC wfns is just frozen-docc */
    moinfo.core = init_int_array(moinfo.nirreps);
    moinfo.ncore = 0;
    for(h=0; h < moinfo.nirreps; h++) {
      moinfo.core[h] = moinfo.frdocc[h];
      moinfo.ncore += moinfo.core[h];
    }

    /* frozen-core offsets used for CC energies */
    moinfo.C_offset = init_int_array(moinfo.nirreps);
    for(h=0; h < moinfo.nirreps; h++)
      moinfo.C_offset[h] = moinfo.frdocc[h];

    /** Compute CC spatial-orbial reordering array(s) for the one-electron transformation **/

    if(params.ref == 0 || params.ref == 1) {
      moinfo.pitz2corr_one = init_int_array(moinfo.nmo);
      reorder_qt(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
             moinfo.pitz2corr_one, moinfo.mopi, moinfo.nirreps);
    }
    else if(params.ref == 2) {
      moinfo.pitz2corr_one_A = init_int_array(moinfo.nmo);
      moinfo.pitz2corr_one_B = init_int_array(moinfo.nmo);
      reorder_qt_uhf(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
             moinfo.pitz2corr_one_A, moinfo.pitz2corr_one_B, moinfo.mopi, moinfo.nirreps);
    }

    /** Compute CC spatial-orbial reordering array(s) for the two-electron transformation **/

    offset = init_int_array(moinfo.nirreps);
    for(h=1; h < moinfo.nirreps; h++)
      offset[h] = offset[h-1] + moinfo.actpi[h-1];

    if(params.ref == 0 || params.ref == 1) {

      moinfo.pitz2corr_two = init_int_array(moinfo.nactive);

      count = 0;
      for(h=0; h < moinfo.nirreps; h++) {
        this_offset = offset[h];
        for(p=0; p < moinfo.clsdpi[h] - moinfo.frdocc[h]; p++)
          moinfo.pitz2corr_two[this_offset+p] = count++;
      }
      for(h=0; h < moinfo.nirreps; h++) {
        this_offset = offset[h] + moinfo.clsdpi[h] - moinfo.frdocc[h];
        for(p=0; p < moinfo.openpi[h]; p++)
          moinfo.pitz2corr_two[this_offset+p] = count++;
      }
      for(h=0; h < moinfo.nirreps; h++) {
        this_offset = offset[h] + moinfo.clsdpi[h] - moinfo.frdocc[h] + moinfo.openpi[h];
        for(p=0; p < moinfo.uoccpi[h] - moinfo.fruocc[h]; p++)
          moinfo.pitz2corr_two[this_offset+p] = count++;
      }
    }
    else if(params.ref == 2) {
      moinfo.pitz2corr_two_A = init_int_array(moinfo.nactive);
      moinfo.pitz2corr_two_B = init_int_array(moinfo.nactive);

      count = 0;
      for(h=0; h < moinfo.nirreps; h++) {
        this_offset = offset[h];
        for(p=0; p < moinfo.clsdpi[h] - moinfo.frdocc[h] + moinfo.openpi[h]; p++)
          moinfo.pitz2corr_two_A[this_offset+p] = count++;
      }
      for(h=0; h < moinfo.nirreps; h++) {
        this_offset = offset[h] + moinfo.clsdpi[h] - moinfo.frdocc[h] + moinfo.openpi[h];
        for(p=0; p < moinfo.uoccpi[h] - moinfo.fruocc[h]; p++)
          moinfo.pitz2corr_two_A[this_offset+p] = count++;
      }

      count = 0;
      for(h=0; h < moinfo.nirreps; h++) {
        this_offset = offset[h];
        for(p=0; p < moinfo.clsdpi[h] - moinfo.frdocc[h]; p++)
          moinfo.pitz2corr_two_B[this_offset+p] = count++;
      }
      for(h=0; h < moinfo.nirreps; h++) {
        this_offset = offset[h] + moinfo.clsdpi[h] - moinfo.frdocc[h];
        for(p=0; p < moinfo.openpi[h] + moinfo.uoccpi[h] - moinfo.fruocc[h]; p++)
          moinfo.pitz2corr_two_B[this_offset+p] = count++;
      }
    }
    free(offset);

      }
      else if((params.wfn == "SCF") || (params.wfn == "SCF_MVD")
           || (params.wfn == "MRCCSD")) {  /* psimrcc requires no freezing nor reordering */

    /* Note that no frozen orbitals are allowed in this case */

    moinfo.nactive = moinfo.nmo;
    moinfo.actpi = init_int_array(moinfo.nirreps);
    for(h=0; h < moinfo.nirreps; h++) moinfo.actpi[h] = moinfo.mopi[h];

    moinfo.actsym = init_int_array(moinfo.nactive);
    for(h=0,count=0; h < moinfo.nirreps; h++)
      for(i=0; i < moinfo.actpi[h]; i++,count++)
        moinfo.actsym[count] = h;

    /* "core" orbitals needed for frozen-core operator */
    /* core is a zero-vector for SCF wfns */
    moinfo.core = init_int_array(moinfo.nirreps);
    moinfo.ncore = 0;

    /* zero offsets for SCF freqs and props */
    moinfo.C_offset = init_int_array(moinfo.nirreps);

    /** Compute SCF spatial-orbial reordering array(s) for the
        one- and two-electron transformations **/
    if(params.ref == 0 || params.ref == 1) {
      moinfo.pitz2corr_one = init_int_array(moinfo.nmo);
      moinfo.pitz2corr_two = init_int_array(moinfo.nmo);
      for(i=0; i < moinfo.nmo; i++) {
        moinfo.pitz2corr_one[i] = moinfo.pitz2corr_two[i] = i;
      }
    }
    else if(params.ref == 2) {
      moinfo.pitz2corr_one_A = init_int_array(moinfo.nmo);
      moinfo.pitz2corr_one_B = init_int_array(moinfo.nmo);
      moinfo.pitz2corr_two_A = init_int_array(moinfo.nmo);
      moinfo.pitz2corr_two_B = init_int_array(moinfo.nmo);
      for(i=0; i < moinfo.nmo; i++) {
        moinfo.pitz2corr_one_A[i] = i;
        moinfo.pitz2corr_one_B[i] = i;
        moinfo.pitz2corr_two_A[i] = i;
        moinfo.pitz2corr_two_B[i] = i;
      }
    }
      }
      else {
        char junk[32];
    sprintf(junk, "WFN %s not yet supported by transqt2.", params.wfn.c_str());
        throw PsiException(junk, __FILE__, __LINE__);
      }

      if(params.print_lvl) {
    fprintf(outfile,"\tChkpt Parameters:\n");
    fprintf(outfile,"\t--------------------\n");
    fprintf(outfile,"\tNumber of irreps     = %d\n",moinfo.nirreps);
    fprintf(outfile,"\tNumber of SOs        = %d\n",moinfo.nso);
    fprintf(outfile,"\tNumber of MOs        = %d\n",moinfo.nmo);
    fprintf(outfile,"\tNumber of active MOs = %d\n\n",moinfo.nactive);
    fprintf(outfile,
        "\tLabel\t# SOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
    fprintf(outfile,
        "\t-----\t-----\t------\t------\t------\t------\t------\n");
    for(i=0; i < moinfo.nirreps; i++) {
      fprintf(outfile,
          "\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n",
          moinfo.labels[i],moinfo.sopi[i],moinfo.frdocc[i],
          moinfo.clsdpi[i]-moinfo.frdocc[i],moinfo.openpi[i],moinfo.uoccpi[i]-moinfo.fruocc[i],
          moinfo.fruocc[i]);
    }
    fprintf(outfile,"\n\tNuclear Rep. energy (chkpt) =  %20.14f\n", moinfo.enuc);
    fprintf(outfile,  "\tSCF energy          (chkpt) =  %20.14f\n", escf);
      }
    }

    void cleanup(void)
    {
      int h;

      for(h=0; h < moinfo.nirreps; h++)
    free(moinfo.labels[h]);
      free(moinfo.labels);

      free(moinfo.C_offset);

      if(ci_wfn(params.wfn)) {
    free(moinfo.pitz2corr_one);
    free(moinfo.pitz2corr_two);
    for(h=0; h < moinfo.nirreps; h++)
      free_block(moinfo.C[h]);
    free(moinfo.C);
      }
      else if(cc_wfn(params.wfn) || (params.wfn == "SCF") ||
          (params.wfn == "SCF_MVD")) {
    if(params.ref == 0 || params.ref == 1) {
      free(moinfo.pitz2corr_one);
      free(moinfo.pitz2corr_two);
      for(h=0; h < moinfo.nirreps; h++)
        if(moinfo.sopi[h] && moinfo.mopi[h]) free_block(moinfo.C[h]);
      free(moinfo.C);
    }
    else if(params.ref == 2) {
      free(moinfo.pitz2corr_one_A);
      free(moinfo.pitz2corr_one_B);
      free(moinfo.pitz2corr_two_A);
      free(moinfo.pitz2corr_two_B);
      for(h=0; h < moinfo.nirreps; h++) {
            if(moinfo.sopi[h] && moinfo.mopi[h]) {
          free_block(moinfo.C_a[h]);
          free_block(moinfo.C_b[h]);
            }
      }
      free(moinfo.C_a); free(moinfo.C_b);
    }
      }

      free(moinfo.sopi);
      free(moinfo.sosym);
      free(moinfo.mopi);
      free(moinfo.mosym);
      free(moinfo.actpi);
      free(moinfo.actsym);
      free(moinfo.clsdpi);
      free(moinfo.openpi);
      free(moinfo.uoccpi);
      // Wavefunction owns these two arrays
//      free(moinfo.frdocc);
//      free(moinfo.fruocc);
      free(moinfo.core);
    }

  } // namespace transqt2
} // namespace psi
