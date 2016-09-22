/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup TRANSQT2
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#define EXTERN
#include "globals.h"

#include "psi4/libmints/wavefunction.h"

/* get_moinfo(): Routine to obtain basic orbital information from
** wfn and compute the associated lookup arrays.
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

void get_moinfo(SharedWavefunction wfn, Options& options)
{
    int i, j, k, h, p, q, count;
    double escf;
    int *offset, this_offset;
    int *rstr_docc, *rstr_uocc, **ras_opi;
    double ***C, ***C_a, ***C_b;
    double **C_full, **C_full_a, **C_full_b;
    int *reorder, *reorder_A, *reorder_B;
    double **TMP;

    moinfo.nirreps = wfn->nirrep();
    moinfo.nmo = wfn->nmo();
    moinfo.nso = wfn->nso();
    moinfo.nao = wfn->basisset()->nao();
    moinfo.labels = wfn->molecule()->irrep_labels();
    moinfo.enuc = wfn->molecule()->nuclear_repulsion_energy();
    escf = wfn->reference_energy(); // Is this the one I want? -TDC, 11/15/15

    // Would like to replace these with Dimension objects, TDC 11/2015
    moinfo.frdocc = init_int_array(moinfo.nirreps);
    moinfo.fruocc = init_int_array(moinfo.nirreps);
    moinfo.sopi = init_int_array(moinfo.nirreps);
    moinfo.mopi = init_int_array(moinfo.nirreps);
    moinfo.clsdpi = init_int_array(moinfo.nirreps);
    moinfo.openpi = init_int_array(moinfo.nirreps);
   for (int h=0; h<moinfo.nirreps; h++) {
     moinfo.frdocc[h] = wfn->frzcpi()[h];
     moinfo.fruocc[h] = wfn->frzvpi()[h];
     moinfo.sopi[h] = wfn->nsopi()[h];
     moinfo.mopi[h] = wfn->nmopi()[h];
     moinfo.clsdpi[h] = wfn->doccpi()[h];
     moinfo.openpi[h] = wfn->soccpi()[h];
   }

    if(options["FROZEN_DOCC"].has_changed()){
        if(options["FROZEN_DOCC"].size() != moinfo.nirreps)
            throw PSIEXCEPTION("FROZEN_DOCC array should be the same size as the number of irreps.");
        for(int h = 0; h < moinfo.nirreps; ++h)
            moinfo.frdocc[h] = options["FROZEN_DOCC"][h].to_integer();
    }
    if(options["FROZEN_UOCC"].has_changed()){
        if(options["FROZEN_UOCC"].size() != moinfo.nirreps)
            throw PSIEXCEPTION("FROZEN_UOCC array should be the same size as the number of irreps.");
        for(int h = 0; h < moinfo.nirreps; ++h)
            moinfo.fruocc[h] = options["FROZEN_UOCC"][h].to_integer();
    }

    moinfo.nfzc = moinfo.nfzv = 0;
    for(i=0; i < moinfo.nirreps; i++) {
        moinfo.nfzc += moinfo.frdocc[i];
        moinfo.nfzv += moinfo.fruocc[i];
    }

    moinfo.uoccpi = init_int_array(moinfo.nirreps);
    for(i=0; i < moinfo.nirreps; i++)
        moinfo.uoccpi[i] = moinfo.mopi[i] - moinfo.clsdpi[i] - moinfo.openpi[i];

    if(params.semicanonical)
      throw PSIEXCEPTION("SemiCanonical transform does not work at the moment");
      //Process::environment.legacy_wavefunction()->semicanonicalize();

    /* SO symmetry array */
    moinfo.sosym = init_int_array(moinfo.nso);
    for(h=0,count=0; h < moinfo.nirreps; h++)
        for(i=0; i < moinfo.sopi[h]; i++,count++)
            moinfo.sosym[count] = h;

    /* AO/MO transformation matrix */
    if(params.ref == 0 || params.ref == 1) { /* RHF/ROHF */
        C = (double ***) malloc(moinfo.nirreps * sizeof(double **));
        for(h=0; h < moinfo.nirreps; h++) {
            C[h] = wfn->Ca()->pointer(h);
            if(params.print_lvl > 2) {
                outfile->Printf( "\n\tMOs for irrep %d:\n",h);
                mat_print(C[h], moinfo.sopi[h], moinfo.mopi[h], "outfile");
            }
        }
        C_full = wfn->Ca()->to_block_matrix();
        moinfo.C = C;
        moinfo.C_full = C_full;
    }

    else if(params.ref == 2) { /* UHF */
        C_a = (double ***) malloc(moinfo.nirreps * sizeof(double **));
        C_b = (double ***) malloc(moinfo.nirreps * sizeof(double **));
        for(h=0; h < moinfo.nirreps; h++) {
            C_a[h] = wfn->Ca()->pointer(h);
            if(params.print_lvl > 2) {
                outfile->Printf( "\n\tAlpha MOs for irrep %d:\n",h);
                mat_print(C_a[h], moinfo.sopi[h], moinfo.mopi[h], "outfile");
            }
            C_b[h] = wfn->Cb()->pointer(h);
            if(params.print_lvl > 2) {
                outfile->Printf( "\n\tBeta MOs for irrep %d:\n",h);
                mat_print(C_b[h], moinfo.sopi[h], moinfo.mopi[h], "outfile");
            }
        }
        C_full_a = wfn->Ca()->to_block_matrix();
        C_full_b = wfn->Cb()->to_block_matrix();
        moinfo.C_a = C_a;
        moinfo.C_b = C_b;
        moinfo.C_full_a = C_full_a;
        moinfo.C_full_b = C_full_b;
    }

    if(cc_wfn(params.wfn) || (params.wfn == "MP2")) {

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
        outfile->Printf("\tChkpt Parameters:\n");
        outfile->Printf("\t--------------------\n");
        outfile->Printf("\tNumber of irreps     = %d\n",moinfo.nirreps);
        outfile->Printf("\tNumber of SOs        = %d\n",moinfo.nso);
        outfile->Printf("\tNumber of MOs        = %d\n",moinfo.nmo);
        outfile->Printf("\tNumber of active MOs = %d\n\n",moinfo.nactive);
        outfile->Printf(
                "\tLabel\t# SOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
        outfile->Printf(
                "\t-----\t-----\t------\t------\t------\t------\t------\n");
        for(i=0; i < moinfo.nirreps; i++) {
            outfile->Printf(
                    "\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n",
                    moinfo.labels[i],moinfo.sopi[i],moinfo.frdocc[i],
                    moinfo.clsdpi[i]-moinfo.frdocc[i],moinfo.openpi[i],moinfo.uoccpi[i]-moinfo.fruocc[i],
                    moinfo.fruocc[i]);
        }
        if (ci_wfn(params.wfn)) {
            outfile->Printf("\n\tDOCC         = ");
            for (i=0; i<moinfo.nirreps; i++)
                outfile->Printf("%2d ", moinfo.clsdpi[i]);
            outfile->Printf("\n\tSOCC         = ");
            for (i=0; i<moinfo.nirreps; i++)
                outfile->Printf("%2d ", moinfo.openpi[i]);
            outfile->Printf("\n");
            outfile->Printf("\n\tFROZEN DOCC  = ");
            for (i=0; i<moinfo.nirreps; i++)
                outfile->Printf("%2d ", moinfo.frdocc[i]);
            outfile->Printf("\n\tRESTR DOCC   = ");
            for (i=0; i<moinfo.nirreps; i++)
                outfile->Printf("%2d ", rstr_docc[i]);
            for (i=0; i<4; i++) {
                outfile->Printf("\n\tRAS %d        = ",i+1);
                for (int j=0; j<moinfo.nirreps; j++)
                    outfile->Printf("%2d ", ras_opi[i][j]);
            }
            outfile->Printf("\n\tRESTR UOCC   = ");
            for (i=0; i<moinfo.nirreps; i++)
                outfile->Printf("%2d ", rstr_uocc[i]);
            outfile->Printf("\n\tFROZEN UOCC  = ");
            for (i=0; i<moinfo.nirreps; i++)
                outfile->Printf("%2d ", moinfo.fruocc[i]);
            outfile->Printf("\n");

            free_int_matrix(ras_opi);
            free(rstr_docc);
            free(rstr_uocc);
        }
        outfile->Printf("\n\tNuclear Rep. energy (wfn)   =  %20.14f\n", moinfo.enuc);
        outfile->Printf(  "\tSCF energy          (wfn)   =  %20.14f\n", escf);
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
//        for(h=0; h < moinfo.nirreps; h++)
//            free_block(moinfo.C[h]);
        free(moinfo.C);
    }
    else if(cc_wfn(params.wfn) || (params.wfn == "SCF") ||
            (params.wfn == "SCF_MVD")) {
        if(params.ref == 0 || params.ref == 1) {
            free(moinfo.pitz2corr_one);
            free(moinfo.pitz2corr_two);
//            for(h=0; h < moinfo.nirreps; h++)
//                if(moinfo.sopi[h] && moinfo.mopi[h]) free_block(moinfo.C[h]);
            free(moinfo.C);
        }
        else if(params.ref == 2) {
            free(moinfo.pitz2corr_one_A);
            free(moinfo.pitz2corr_one_B);
            free(moinfo.pitz2corr_two_A);
            free(moinfo.pitz2corr_two_B);
//            for(h=0; h < moinfo.nirreps; h++) {
//                if(moinfo.sopi[h] && moinfo.mopi[h]) {
 //                   free_block(moinfo.C_a[h]);
//                    free_block(moinfo.C_b[h]);
//                }
//            }
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
    // Wavefunction did own the following 2 arrays, but not in current test
    free(moinfo.frdocc);
    free(moinfo.fruocc);
    free(moinfo.core);
}

} // namespace transqt2
} // namespace psi
