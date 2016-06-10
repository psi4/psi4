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

/*! \defgroup TRANSQT2 transqt2: Integral Transformation Program */


/*!
** \file
** \ingroup TRANSQT2
** \brief Integral Transformation Program
**
** A program to transform one- and two-electron integrals from the
** symmetry-orbital basis to the molecular-orbital basis.
**
** This code replaces the original transqt code developed initially
** in 1995 by TDC, CDS, and JTF.  This version is designed to take
** advantage of libdpd's ability to handle easily four-index
** quantities, including symmetry.  This version requires
** significantly less disk space (ca. 1/2) than the original code,
** and is often much faster because of its reduced I/O requirements.
**
** This version of the code can do RHF, ROHF, and UHF transformations
** that are compatible with all the coupled cluster codes, including
** frozen orbitals.
**
** Remaining tasks to achieve full replacement of transqt:
**   (1) Add reordering arrays needed for DETCI and SCF DERTYPE=2. (DONE)
**   (2) Add partial transforms for MP2 and MP2-R12. (Still needed?)
**   (3) Replace the backtransformation.  (I want to do this with
**       symmetry, though, so there's no hurry here.)
**
** TDC, 7/06 (last updated 2/08)
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include <psi4-dec.h>
#include "globals.h"

#include <libmints/wavefunction.h>
#include <libtrans/mospace.h>
#include <libmints/matrix.h>
#include <libpsio/psio.hpp>

namespace psi {
namespace transqt2 {

void init_io();
void title(void);
void get_params(Options& options);
void get_moinfo(SharedWavefunction ref_wfn, Options& options);
void cleanup(void);
void exit_io(void);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_rhf(int **);
int file_build_presort(dpdfile4 *,int,double,long int,int,int,double *,
                       double *,double *,double *,int);
void transtwo_rhf(void);
void transtwo_uhf(void);
void transone(int,int,double *,double *,double **,int,int *);

PsiReturnType transqt2(SharedWavefunction ref_wfn, Options & options)
{
    int nso, nmo, ntri_so, ntri_mo, nirreps;
    int **cachelist, *cachefiles;
    dpdfile4 I;
    int h, pq, p, q, i;
    double *H, *D, *F, *oei;
    double *H_a, *H_b, *D_a, *D_b, *F_a, *F_b;
    double **C, **C_a, **C_b;
    int stat;
    int *so_offset, *mo_offset;
    double efzc;

    init_io();
    title();

    get_params(options);
    get_moinfo(ref_wfn, options);

    nso = moinfo.nso;
    nmo = moinfo.nmo;
    ntri_so = nso*(nso+1)/2;
    ntri_mo = nmo*(nmo+1)/2;
    nirreps = moinfo.nirreps;

    cachefiles = init_int_array(PSIO_MAXUNIT);
    cachelist = cacheprep_rhf(params.cachelev, cachefiles); /* really just a placeholder */

    std::vector<int*> spaces;
    spaces.push_back(moinfo.sopi);
    spaces.push_back(moinfo.sosym);
    spaces.push_back(moinfo.actpi);
    spaces.push_back(moinfo.actsym);
    dpd_init(0, nirreps, params.memory, 0, cachefiles, cachelist, NULL, 2, spaces);

    /*** Starting one-electron transforms and presort ***/

    if(params.ref == 0 || params.ref == 1) C = moinfo.C_full;
    else {
        C_a = moinfo.C_full_a;
        C_b = moinfo.C_full_b;
    }

    /* build the frozen-core density (RHF) */
    so_offset = init_int_array(nirreps);
    mo_offset = init_int_array(nirreps);
    for(h=1; h < nirreps; h++) {
        so_offset[h] = so_offset[h-1] + moinfo.sopi[h-1];
        mo_offset[h] = mo_offset[h-1] + moinfo.mopi[h-1];
    }

    if(params.ref == 0 || params.ref == 1) { /* RHF/ROHF */
        D = init_array(ntri_so);
        for(h=0; h < nirreps; h++)
            for(p=so_offset[h]; p < so_offset[h]+moinfo.sopi[h]; p++)
                for(q=so_offset[h]; q <=p; q++) {
                    pq = INDEX(p,q);
                    for(i=mo_offset[h]; i < mo_offset[h] + moinfo.core[h]; i++)
                        D[pq] += C[p][i] * C[q][i];
                }
        if(params.print_lvl > 2) {
            outfile->Printf( "\n\tFrozen-core density (SO):\n");
            print_array(D, nso, "outfile");
        }
    }
    else { /* UHF */
        D_a = init_array(ntri_so);
        D_b = init_array(ntri_so);
        for(h=0; h < nirreps; h++)
            for(p=so_offset[h]; p < so_offset[h]+moinfo.sopi[h]; p++)
                for(q=so_offset[h]; q <=p; q++) {
                    pq = INDEX(p,q);
                    for(i=mo_offset[h]; i < mo_offset[h] + moinfo.core[h]; i++) {
                        D_a[pq] += C_a[p][i] * C_a[q][i];
                        D_b[pq] += C_b[p][i] * C_b[q][i];
                    }
                }
        if(params.print_lvl > 2) {
            outfile->Printf( "\n\tAlpha Frozen-core density (SO):\n");
            print_array(D_a, nso, "outfile");
            outfile->Printf( "\n\tBeta Frozen-core density (SO):\n");
            print_array(D_b, nso, "outfile");
        }
    }

    free(so_offset);
    free(mo_offset);

    /* pre-sort the SO-basis two-electron integrals and generate the fzc operator(s) */
    if(params.ref == 0 || params.ref == 1)
        F = init_array(ntri_so);
    else {
        F_a = init_array(ntri_so);
        F_b = init_array(ntri_so);
    }

    timer_on("presort");
    if(params.print_lvl) {
        outfile->Printf( "\n\tPresorting SO-basis two-electron integrals.\n");
        
    }
    psio_open(PSIF_SO_PRESORT, 0);
    global_dpd_->file4_init(&I, PSIF_SO_PRESORT, 0, 3, 3, "SO Ints (pq,rs)");
    if(params.ref == 0 || params.ref == 1)
        file_build_presort(&I, PSIF_SO_TEI, params.tolerance, params.memory,
                           !params.delete_tei, moinfo.ncore, D, NULL, F, NULL, params.ref);
    else
        file_build_presort(&I, PSIF_SO_TEI, params.tolerance, params.memory,
                           !params.delete_tei, moinfo.ncore, D_a, D_b, F_a, F_b, params.ref);
    global_dpd_->file4_close(&I);
    psio_close(PSIF_SO_PRESORT, 1);
    timer_off("presort");

    /* read the bare one-electron integrals */

    boost::shared_ptr<PSIO> psio_;
    psio_ = ref_wfn->psio();

    SharedMatrix H_so = ref_wfn->H()->clone();
    H_so->set_name(PSIF_SO_H);
    H_so->save(psio_, PSIF_OEI);

    oei = init_array(ntri_so);
    H = init_array(ntri_so);
    stat = iwl_rdone(PSIF_OEI, PSIF_SO_H, H, ntri_so, 0, 0, "outfile");

    /* add the remaining one-electron terms to the fzc operator(s) */
    if(params.ref == 0 || params.ref == 1) {
        for(pq=0; pq < ntri_so; pq++)
            F[pq] += H[pq];
    }
    else {
        for(pq=0; pq < ntri_so; pq++) {
            F_a[pq] += H[pq];
            F_b[pq] += H[pq];
        }
    }

    /* compute the frozen-core energy and write it to the wfn */
    efzc = 0.0;
    if(params.ref == 0 || params.ref == 1) { /* RHF/ROHF */
        for(p=0; p < nso; p++) {
            pq = INDEX(p,p);
            efzc += D[pq] * (H[pq] + F[pq]);
            for(q=0; q < p; q++) {
                pq = INDEX(p,q);
                efzc += 2.0 * D[pq] * (H[pq] + F[pq]);
            }
        }
    }
    else { /* UHF */
        for(p=0; p < nso; p++) {
            pq = INDEX(p,p);
            efzc += 0.5 * D_a[pq] * (H[pq] + F_a[pq]);
            efzc += 0.5 * D_b[pq] * (H[pq] + F_b[pq]);
            for(q=0; q < p; q++) {
                pq = INDEX(p,q);
                efzc += D_a[pq] * (H[pq] + F_a[pq]);
                efzc += D_b[pq] * (H[pq] + F_b[pq]);
            }
        }
    }
    if(params.print_lvl) {
        outfile->Printf( "\tFrozen-core energy = %20.15f\n", efzc);
        
    }
    // Add frozen-core energy to wfn
    Process::environment.legacy_wavefunction()->set_efzc(efzc);

    /*** One-electron forward transforms.  Note that all orbitals are
       transformed, including those in the inactive space. ***/

    /* transform the bare one-electron integrals */
    if(params.ref == 0 || params.ref == 1) {
        transone(nso, nmo, H, oei, C, nmo, moinfo.pitz2corr_one);
        if(params.print_lvl > 2) {
            outfile->Printf( "\n\tOne-electron integrals (MO basis):\n");
            print_array(oei, nmo, "outfile");
        }
        iwl_wrtone(PSIF_OEI, PSIF_MO_OEI, ntri_mo, oei);
    }
    else { /* UHF */
        /* alpha */
        transone(nso, nmo, H, oei, C_a, nmo, moinfo.pitz2corr_one_A);
        if(params.print_lvl > 2) {
            outfile->Printf( "\n\tAlpha one-electron integrals (MO basis):\n");
            print_array(oei, nmo, "outfile");
        }
        iwl_wrtone(PSIF_OEI, PSIF_MO_A_OEI, ntri_mo, oei);

        /* beta */
        transone(nso, nmo, H, oei, C_b, nmo, moinfo.pitz2corr_one_B);
        if(params.print_lvl > 2) {
            outfile->Printf( "\n\tBeta one-electron integrals (MO basis):\n");
            print_array(oei, nmo, "outfile");
        }
        iwl_wrtone(PSIF_OEI, PSIF_MO_B_OEI, ntri_mo, oei);
    }

    /* transform the frozen-core operator */
    if(params.ref == 0 || params.ref == 1) { /* RHF/ROHF */
        transone(nso, nmo, F, oei, C, nmo, moinfo.pitz2corr_one);
        if(params.print_lvl > 2) {
            outfile->Printf( "\n\tFrozen-core operator (MO basis):\n");
            print_array(oei, nmo, "outfile");
        }
        iwl_wrtone(PSIF_OEI, PSIF_MO_FZC, ntri_mo, oei);
    }
    else { /* UHF */

        /* alpha */
        transone(nso, nmo, F_a, oei, C_a, nmo, moinfo.pitz2corr_one_A);
        if(params.print_lvl > 2) {
            outfile->Printf( "\n\tAlpha frozen-core operator (MO basis):\n");
            print_array(oei, nmo, "outfile");
        }
        iwl_wrtone(PSIF_OEI, PSIF_MO_A_FZC, ntri_mo, oei);

        /* beta */
        transone(nso, nmo, F_b, oei, C_b, nmo, moinfo.pitz2corr_one_B);
        if(params.print_lvl > 2) {
            outfile->Printf( "\n\tBeta frozen-core operator (MO basis):\n");
            print_array(oei, nmo, "outfile");
        }
        iwl_wrtone(PSIF_OEI, PSIF_MO_B_FZC, ntri_mo, oei);
    }

    free(oei);
    free(H);
    if(params.ref == 0 || params.ref == 1) {
        free(F);
        free(D);
    }
    else {
        free(F_a);
        free(F_b);
        free(D_a);
        free(D_b);
    }

    /*** One-electron transforms complete ***/

    /*** Starting two-electron transforms ***/

    if(!options.get_bool("NO_TEI")){
        if(params.ref == 0 || params.ref == 1)
            transtwo_rhf();
        else transtwo_uhf();
    }

    /*** Two-electron transforms complete ***/

    dpd_close(0);

    cachedone_rhf(cachelist);
    free(cachefiles);

    cleanup();
    exit_io();
    return(Success);
}

void init_io()
{
    //  int i;
    //  char *progid;
    //  int num_extra_args = 0;
    //  char **extra_args;
    //  extra_args = (char **) malloc(argc*sizeof(char *));

    params.print_lvl = 1;
    params.backtr = 0;
    //  for (i=1; i<argc; i++) {
    //    if (!strcmp(argv[i], "--quiet"))
    //      params.print_lvl = 0;
    //    else if(!strcmp(argv[i], "--backtr"))
    //      params.backtr = 1;
    //    else
    //      extra_args[num_extra_args++] = argv[i];
    //  }

    if(params.print_lvl) tstart();

    psio_open(PSIF_CC_INFO, PSIO_OPEN_NEW);
}

void title(void)
{
    if(params.print_lvl) {
        outfile->Printf( "\n");
        outfile->Printf("\t**************************************************\n");
        outfile->Printf("\t* TRANSQT2: Program to transform integrals from  *\n");
        outfile->Printf("\t*           the SO basis to the MO basis.        *\n");
        outfile->Printf("\t*                                                *\n");
        outfile->Printf("\t*            Daniel, David, & Justin             *\n");
        outfile->Printf("\t**************************************************\n");
        outfile->Printf( "\n");
    }
}

void exit_io(void)
{
    psio_close(PSIF_CC_INFO,1);
    if(params.print_lvl) tstop();
}

} // namespace transqt2
} // namespace psi