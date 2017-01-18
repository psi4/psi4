/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

#include <cmath>
#include <algorithm>

#include "psi4/libmoinfo/libmoinfo.h"
#include "transform.h"
#include "matrix.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "algebra_interface.h"
#include "blas.h"

#define CCTRANSFORM_USE_BLAS

#define MAX(i,j) ((i>j) ? i : j)
#define MIN(i,j) ((i>j) ? j : i)
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define four(i,j,k,l) INDEX(INDEX(i,j),INDEX(k,l))

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/psifiles.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

extern MOInfo *moinfo;

using namespace std;

void CCTransform::read_integrals_mrpt2(IntegralTransform *ints)
{
    read_oei_mo_integrals_mrpt2();
    read_tei_mo_integrals_mrpt2(ints);
}

/**
 * Read the one electron MO integrals from an iwl buffer assuming Pitzer ordering and store them in oei_mo
 */
void CCTransform::read_oei_mo_integrals_mrpt2()
{
    read_oei_so_integrals();
    transform_oei_so_integrals();
}

/**
 * Read the two electron MO integrals from an iwl buffer assuming Pitzer ordering and store them in the packed array tei_mo
 */
void CCTransform::read_tei_mo_integrals_mrpt2(IntegralTransform *ints)
{
#define ID(x) ints->DPD_ID(x)
    // Read all the (frozen + non-frozen) TEI in Pitzer order
    // and store them in a in-core block-matrix
    //   CCIndex* mo_indexing = blas->get_index("[n>=n]");

    dpdbuf4 I;

    dpd_set_default(ints->get_dpd_id());
    std::shared_ptr<PSIO> psio(_default_psio_lib_);
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    std::vector<int> mopi   = moinfo->get_mopi();
    std::vector<int> doccpi = moinfo->get_docc();
    std::vector<int> foccpi = moinfo->get_focc();
    int nirreps = moinfo->get_nirreps();
    std::vector<int> offsets(nirreps, 0);
    int offset = 0;
    for(int h = 1; h < nirreps; ++h){
        offset += mopi[h-1];
        offsets[h] = offset;
    }
    size_t elements = 0;

    // Loop over the DPD buffers for the two integral classes (exchange and then coulomb).

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[M,A]"), ID("[M,A]"),
                           ID("[M,A]"),  ID("[M,A]"), 0, "MO Ints (MA|MA)");
    for(int h = 0; h < I.params->nirreps; ++h){
        global_dpd_->buf4_mat_irrep_init(&I, h);
        global_dpd_->buf4_mat_irrep_rd(&I, h);
        for(int pq = 0; pq < I.params->rowtot[h]; ++pq){
            int p = I.params->roworb[h][pq][0];
            int q = I.params->roworb[h][pq][1];
            int psym = I.params->psym[p];
            int qsym = I.params->qsym[q];
            int prel = p - I.params->poff[psym];
            int qrel = q - I.params->qoff[qsym];
            for(int rs = 0; rs <= pq; ++rs){
                int r = I.params->colorb[h][rs][0];
                int s = I.params->colorb[h][rs][1];
                int rsym = I.params->rsym[r];
                int ssym = I.params->ssym[s];
                int rrel = r - I.params->roff[rsym];
                int srel = s - I.params->soff[ssym];
                int pidx = offsets[psym] + prel;
                int qidx = offsets[qsym] + qrel;
                int ridx = offsets[rsym] + rrel;
                int sidx = offsets[ssym] + srel;
                integral_map[four(pidx,qidx,ridx,sidx)] = I.matrix[h][pq][rs];
                ++elements;
            }
        }
        global_dpd_->buf4_mat_irrep_close(&I, h);
    }
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[M>=M]+"), ID("[E>=E]+"),
                           ID("[M>=M]+"),  ID("[E>=E]+"), 0, "MO Ints (MM|EE)");
    for(int h = 0; h < I.params->nirreps; ++h){
        global_dpd_->buf4_mat_irrep_init(&I, h);
        global_dpd_->buf4_mat_irrep_rd(&I, h);
        for(int pq = 0; pq < I.params->rowtot[h]; ++pq){
            int p = I.params->roworb[h][pq][0];
            int q = I.params->roworb[h][pq][1];
            int psym = I.params->psym[p];
            int qsym = I.params->qsym[q];
            int prel = p - I.params->poff[psym];
            int qrel = q - I.params->qoff[qsym];
            for(int rs = 0; rs < I.params->coltot[h]; ++rs){
                int r = I.params->colorb[h][rs][0];
                int s = I.params->colorb[h][rs][1];
                int rsym = I.params->rsym[r];
                int ssym = I.params->ssym[s];
                int rrel = r - I.params->roff[rsym];
                int srel = s - I.params->soff[ssym];
                int pidx = offsets[psym] + prel;
                int qidx = offsets[qsym] + qrel;
                int ridx = offsets[rsym] + foccpi[rsym] + doccpi[rsym] + rrel;
                int sidx = offsets[ssym] + foccpi[ssym] + doccpi[ssym] + srel;
                integral_map[four(pidx,qidx,ridx,sidx)] = I.matrix[h][pq][rs];
                ++elements;
            }
        }
        global_dpd_->buf4_mat_irrep_close(&I, h);
    }
    global_dpd_->buf4_close(&I);

    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    outfile->Printf("\n    CCTransform: read %lu non-zero integrals (MRPT2)", elements);
}

double CCTransform::tei_mrpt2(int p, int q, int r, int s)
{
    //   outfile->Printf("\n  (%2d %2d|%2d %2d) = %20.15f",p,q,r,s,integral_map[four(p,q,r,s)]);
    return(integral_map[four(p,q,r,s)]);
}

}} /* End Namespaces */
