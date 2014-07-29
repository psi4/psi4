/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
/*
** sort_oei(): Sort the one-electron integrals from QT Standard to CC
** ordering.
**
** There are separate, but similar versions of this routine for
** RHF/ROHF and UHF reference wave functions.
**
** TDC, 6/01
**
** Notes: 
**
** (1) The output of the one-electron integrals from transqt and
** the method by which we read them from disk is somewhat confusing
** due to the treatment of frozen orbitals.  From transqt both the
** bare one-electron integrals and frozen-core operator are given in
** the full MO basis (cf. transform_one.c of transqt), but, when we
** call iwl_rdone() here with non-zero values for nfzc and nfzv, the
** frozen indices are automagically filtered out (cf. iwl_rdone.c) and
** we are left with only active orbitals.  It is inconsistent to
** require iwl_rdone() to do this filtering, in my opinion.
**
** (2) When core orbitals are frozen, we must read the frozen-core
** operator instead of the bare one-electron integrals here in order
** to correctly build the Fock matrix later on.
*/

#include <cstdio>
#include <cstdlib>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void sort_oei_rhf(void);
void sort_oei_uhf(void);

void sort_oei(void)
{
  if(params.ref == 2) sort_oei_uhf();
  else sort_oei_rhf();
}

void sort_oei_uhf(void)
{
  int nirreps, nmo, nactive, ntri_all, ntri_act;
  int p, q, pq, pnew, qnew, psym, qsym;
  int *aocc, *bocc, *aocc_off, *bocc_off, *aocc_sym, *bocc_sym;
  int *avir, *bvir, *avir_off, *bvir_off, *avir_sym, *bvir_sym;
  int *cc_aocc, *cc_bocc, *cc_avir, *cc_bvir;
  double *a_oei, *b_oei, *tmp_oei;
  dpdfile2 hIJ, hij, hAB, hab, hIA, hia;

  nirreps = moinfo.nirreps; 
  nactive = moinfo.nactive;
  nmo = moinfo.nmo;
  aocc = moinfo.aocc;  bocc = moinfo.bocc;
  avir = moinfo.avir; bvir = moinfo.bvir;
  aocc_off = moinfo.aocc_off; bocc_off = moinfo.bocc_off;
  aocc_sym = moinfo.aocc_sym; bocc_sym = moinfo.bocc_sym;
  avir_off = moinfo.avir_off; bvir_off = moinfo.bvir_off;
  avir_sym = moinfo.avir_sym; bvir_sym = moinfo.bvir_sym;
  cc_aocc = moinfo.cc_aocc; cc_bocc = moinfo.cc_bocc;
  cc_avir = moinfo.cc_avir; cc_bvir = moinfo.cc_bvir;

  ntri_all = nmo * (nmo+1)/2;
  ntri_act = nactive * (nactive+1)/2;
  tmp_oei = init_array(ntri_all);
  a_oei = init_array(ntri_act);
  b_oei = init_array(ntri_act);
  iwl_rdone(PSIF_OEI, PSIF_MO_A_FZC, tmp_oei, ntri_all, 0, 0, outfile);
  filter(tmp_oei,a_oei,ioff,nmo,moinfo.nfzc,moinfo.nfzv);
  iwl_rdone(PSIF_OEI, PSIF_MO_B_FZC, tmp_oei, ntri_all, 0, 0, outfile);
  filter(tmp_oei,b_oei,ioff,nmo,moinfo.nfzc,moinfo.nfzv);
  free(tmp_oei);

  global_dpd_->file2_init(&hIJ, PSIF_CC_OEI, 0, 0, 0, "h(I,J)");
  global_dpd_->file2_init(&hij, PSIF_CC_OEI, 0, 2, 2, "h(i,j)");
  global_dpd_->file2_init(&hAB, PSIF_CC_OEI, 0, 1, 1, "h(A,B)");
  global_dpd_->file2_init(&hab, PSIF_CC_OEI, 0, 3, 3, "h(a,b)");
  global_dpd_->file2_init(&hIA, PSIF_CC_OEI, 0, 0, 1, "h(I,A)");
  global_dpd_->file2_init(&hia, PSIF_CC_OEI, 0, 2, 3, "h(i,a)");

  global_dpd_->file2_mat_init(&hIJ);
  global_dpd_->file2_mat_init(&hij);
  global_dpd_->file2_mat_init(&hAB);
  global_dpd_->file2_mat_init(&hab);
  global_dpd_->file2_mat_init(&hIA);
  global_dpd_->file2_mat_init(&hia);

  /* Loop over alpha QT indices and convert to CC ordering */
  for(p=0; p < nactive; p++) {
    for(q=0; q < nactive; q++) {
      pq = INDEX(p,q);

      if(aocc[p] && aocc[q]) {
	pnew = cc_aocc[p];  qnew = cc_aocc[q];
	psym = aocc_sym[pnew]; qsym = aocc_sym[qnew];
	pnew -= aocc_off[psym];  qnew -= aocc_off[qsym];
	if(psym == qsym) hIJ.matrix[psym][pnew][qnew] = a_oei[pq];
      }

      if(bocc[p] && bocc[q]) {
	pnew = cc_bocc[p];  qnew = cc_bocc[q];
	psym = bocc_sym[pnew]; qsym = bocc_sym[qnew];
	pnew -= bocc_off[psym];  qnew -= bocc_off[qsym];
	if(psym == qsym) hij.matrix[psym][pnew][qnew] = b_oei[pq];
      }

      if(avir[p] && avir[q]) {
	pnew = cc_avir[p];  qnew = cc_avir[q];
	psym = avir_sym[pnew]; qsym = avir_sym[qnew];
	pnew -= avir_off[psym];  qnew -= avir_off[qsym];
	if(psym == qsym) hAB.matrix[psym][pnew][qnew] = a_oei[pq];
      }

      if(bvir[p] && bvir[q]) {
	pnew = cc_bvir[p];  qnew = cc_bvir[q];
	psym = bvir_sym[pnew]; qsym = bvir_sym[qnew];
	pnew -= bvir_off[psym];  qnew -= bvir_off[qsym];
	if(psym == qsym) hab.matrix[psym][pnew][qnew] = b_oei[pq];
      }

      if(aocc[p] && avir[q]) {
	pnew = cc_aocc[p];  qnew = cc_avir[q];
	psym = aocc_sym[pnew]; qsym = avir_sym[qnew];
	pnew -= aocc_off[psym];  qnew -= avir_off[qsym];
	if(psym == qsym) hIA.matrix[psym][pnew][qnew] = a_oei[pq];
      }

      if(bocc[p] && bvir[q]) {
	pnew = cc_bocc[p];  qnew = cc_bvir[q];
	psym = bocc_sym[pnew]; qsym = bvir_sym[qnew];
	pnew -= bocc_off[psym];  qnew -= bvir_off[qsym];
	if(psym == qsym) hia.matrix[psym][pnew][qnew] = b_oei[pq];
      }

    }
  }

  global_dpd_->file2_mat_wrt(&hIJ);
  global_dpd_->file2_mat_wrt(&hij);
  global_dpd_->file2_mat_wrt(&hAB);
  global_dpd_->file2_mat_wrt(&hab);
  global_dpd_->file2_mat_wrt(&hIA);
  global_dpd_->file2_mat_wrt(&hia);

  global_dpd_->file2_mat_close(&hIJ);
  global_dpd_->file2_mat_close(&hij);
  global_dpd_->file2_mat_close(&hAB);
  global_dpd_->file2_mat_close(&hab);
  global_dpd_->file2_mat_close(&hIA);
  global_dpd_->file2_mat_close(&hia);

  global_dpd_->file2_close(&hIJ);
  global_dpd_->file2_close(&hij);
  global_dpd_->file2_close(&hAB);
  global_dpd_->file2_close(&hab);
  global_dpd_->file2_close(&hIA);
  global_dpd_->file2_close(&hia);

  free(a_oei);
  free(b_oei);
}

void sort_oei_rhf(void)
{
  int nirreps, nactive, nmo, ntri_all, ntri_act;
  int p, q, pq, pnew, qnew, psym, qsym;
  int *occ, *vir;
  int *cc_occ, *cc_vir;
  int *occ_sym, *vir_sym;
  int *occ_off, *vir_off;
  double *oei, *tmp_oei;
  dpdfile2 Hoo, Hov, Hvv;

  nirreps = moinfo.nirreps; 
  nactive = moinfo.nactive;
  nmo = moinfo.nmo;
  occ = moinfo.occ; 
  vir = moinfo.vir;
  cc_occ = moinfo.cc_occ; 
  cc_vir = moinfo.cc_vir;
  occ_sym = moinfo.occ_sym; 
  vir_sym = moinfo.vir_sym;
  occ_off = moinfo.occ_off; 
  vir_off = moinfo.vir_off;

  /* Grab the frozen-core opertor */
  ntri_all = nmo * (nmo + 1)/2;
  ntri_act = nactive * (nactive + 1)/2;

  tmp_oei = init_array(ntri_all);
  oei = init_array(ntri_act);
  iwl_rdone(PSIF_OEI, PSIF_MO_FZC, tmp_oei, ntri_all, 0, 0, outfile);
  filter(tmp_oei, oei, ioff, nmo, moinfo.nfzc, moinfo.nfzv);
  free(tmp_oei);

  if(params.print_lvl > 5) {
    psi::fprintf(outfile, "\n\tFrozen-Core Operator:\n");
    psi::fprintf(outfile,   "\t---------------------");
    print_array(oei, nactive, outfile);
  }

  global_dpd_->file2_init(&Hoo, PSIF_CC_OEI, 0, 0, 0, "h(i,j)");
  global_dpd_->file2_init(&Hvv, PSIF_CC_OEI, 0, 1, 1, "h(a,b)");
  global_dpd_->file2_init(&Hov, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");

  global_dpd_->file2_mat_init(&Hoo);
  global_dpd_->file2_mat_init(&Hvv);
  global_dpd_->file2_mat_init(&Hov);

  /* Loop over QT indices and convert to CC ordering */
  for(p=0; p < nactive; p++) {
    for(q=0; q < nactive; q++) {
      pq = INDEX(p, q);

      /* Check occ-occ class */
      if(occ[p] && occ[q]) {
	/* Get relative indices */
	pnew = cc_occ[p];  qnew = cc_occ[q];
	/* Get orbital symmetries */
	psym = occ_sym[pnew]; qsym = occ_sym[qnew];
	/* Shift symmetry-relative indices */
	pnew -= occ_off[psym]; qnew -= occ_off[qsym];
	/* Check orbital symmetry and put integral in place */
	if(psym == qsym) {
	  Hoo.matrix[psym][pnew][qnew] = oei[pq];
	}
      }

      /* Check vir-vir class */
      if(vir[p] && vir[q]) {
	/* Get relative indices */
	pnew = cc_vir[p]; qnew = cc_vir[q];
	/* Get orbital symmetries */
	psym = vir_sym[pnew]; qsym = vir_sym[qnew];
	/* Shift symmetry-relative indices */
	pnew -= vir_off[psym]; qnew -= vir_off[qsym];
	/* Check orbital symmetry and put integral in place */
	if(psym == qsym) {
	  Hvv.matrix[psym][pnew][qnew] = oei[pq];
	}
      }

      /* Check occ-vir class */
      if(occ[p] && vir[q]) {
	/* Get relative indices */
	pnew = cc_occ[p]; qnew = cc_vir[q];
	/* Get orbital symmetries */
	psym = occ_sym[pnew]; qsym = vir_sym[qnew];
	/* Shift symmetry-relative indices */
	pnew -= occ_off[psym]; qnew -= vir_off[qsym];
	/* Check orbital symmetry and put integral in place */
	if(psym == qsym) {
	  Hov.matrix[psym][pnew][qnew] = oei[pq];
	}
      }
    }
  }

  global_dpd_->file2_mat_wrt(&Hoo);
  global_dpd_->file2_mat_wrt(&Hvv);
  global_dpd_->file2_mat_wrt(&Hov);

  global_dpd_->file2_mat_close(&Hoo);
  global_dpd_->file2_mat_close(&Hvv);
  global_dpd_->file2_mat_close(&Hov);

  global_dpd_->file2_close(&Hoo);
  global_dpd_->file2_close(&Hvv);
  global_dpd_->file2_close(&Hov);

  free(oei);
}

}} // namespace psi::ccsort
