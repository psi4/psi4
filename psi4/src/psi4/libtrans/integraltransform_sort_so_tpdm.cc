/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "integraltransform.h"
#include "psi4/psi4-dec.h"
#include "psi4/libdpd/dpd.h"

#include "psi4/psifiles.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/integral.h"
#ifdef INDEX2
#undef INDEX2
#endif
#define INDEX2(i, j) ((i) > (j) ? (i) * ((i) + 1) / 2 + (j) : (j) * ((j) + 1) / 2 + (i))

namespace psi {

void IntegralTransform::setup_tpdm_buffer(const dpdbuf4 *D) {
    auto PQIter = std::make_shared<SO_PQ_Iterator>(sobasis_);
    tpdm_buffer_sizes_.clear();
    size_t max_size = 0;
    for (PQIter->first(); PQIter->is_done() == false; PQIter->next()) {
        int p = PQIter->p();
        int q = PQIter->q();
        std::shared_ptr<SO_RS_Iterator> RSIter =
            std::make_shared<SO_RS_Iterator>(p, q, sobasis_, sobasis_, sobasis_, sobasis_);
        size_t count = 0;
        for (RSIter->first(); RSIter->is_done() == false; RSIter->next()) {
            int ish = RSIter->p();
            int jsh = RSIter->q();
            int ksh = RSIter->r();
            int lsh = RSIter->s();

            int n1 = sobasis_->nfunction(ish);
            int n2 = sobasis_->nfunction(jsh);
            int n3 = sobasis_->nfunction(ksh);
            int n4 = sobasis_->nfunction(lsh);

            // The starting orbital for each irrep can be grabbed from DPD
            int *sym_offsets = D->params->poff;
            for (int itr = 0; itr < n1; itr++) {
                int ifunc = sobasis_->function(ish) + itr;
                int isym = sobasis_->irrep(ifunc);
                int irel = sobasis_->function_within_irrep(ifunc);
                int iabs = sym_offsets[isym] + irel;
                for (int jtr = 0; jtr < n2; jtr++) {
                    int jfunc = sobasis_->function(jsh) + jtr;
                    int jsym = sobasis_->irrep(jfunc);
                    int jrel = sobasis_->function_within_irrep(jfunc);
                    int jabs = sym_offsets[jsym] + jrel;
                    for (int ktr = 0; ktr < n3; ktr++) {
                        int kfunc = sobasis_->function(ksh) + ktr;
                        int ksym = sobasis_->irrep(kfunc);
                        int krel = sobasis_->function_within_irrep(kfunc);
                        int kabs = sym_offsets[ksym] + krel;
                        for (int ltr = 0; ltr < n4; ltr++) {
                            int lfunc = sobasis_->function(lsh) + ltr;
                            int lsym = sobasis_->irrep(lfunc);
                            int lrel = sobasis_->function_within_irrep(lfunc);
                            int labs = sym_offsets[lsym] + lrel;

                            if (isym ^ jsym ^ ksym ^ lsym) continue;  // Not totally symmetric

                            if (ish == jsh) {
                                if (iabs < jabs) continue;

                                if (ksh == lsh) {
                                    if (kabs < labs) continue;
                                    if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                        if (ish == ksh)  // IIII case
                                            continue;
                                    }
                                }
                            } else {
                                if (ksh == lsh) {  // IJKK case
                                    if (kabs < labs) continue;
                                } else {  // IJIJ case
                                    if (ish == ksh && jsh == lsh && INDEX2(iabs, jabs) < INDEX2(kabs, labs)) continue;
                                }
                            }
                            ++count;
                        }
                    }
                }
            }
        }  // End rs iterator
        max_size = count > max_size ? count : max_size;
        tpdm_buffer_sizes_.push_back(count);
    }  // End pq iterator
    size_t num_pairs = tpdm_buffer_sizes_.size();
    psio_->write_entry(PSIF_AO_TPDM, "Num. Pairs", (char *)&num_pairs, sizeof(size_t));
    tpdm_buffer_ = new double[max_size];
    size_t *temp = new size_t[num_pairs];
    for (int i = 0; i < num_pairs; ++i) temp[i] = tpdm_buffer_sizes_[i];
    psio_->write_entry(PSIF_AO_TPDM, "TPDM Buffer Sizes", (char *)temp, num_pairs * sizeof(size_t));
    delete[] temp;
}

void IntegralTransform::sort_so_tpdm(const dpdbuf4 *D, int irrep, size_t first_row, size_t num_rows, bool first_run) {
    // The buffer needs to be set up if the pointer is still null
    if (tpdm_buffer_ == nullptr) setup_tpdm_buffer(D);

    size_t last_row = first_row + num_rows;
    size_t pq_pair_count = 0;
    auto PQIter = std::make_shared<SO_PQ_Iterator>(sobasis_);
    char *toc = new char[40];
    for (PQIter->first(); PQIter->is_done() == false; PQIter->next()) {
        int p = PQIter->p();
        int q = PQIter->q();
        
        sprintf(toc, "SO_TPDM_FOR_PAIR_%zd", pq_pair_count);
        size_t buffer_size = tpdm_buffer_sizes_[pq_pair_count];
        if (first_run)
            ::memset((void *)tpdm_buffer_, '\0', buffer_size * sizeof(double));
        else
            psio_->read_entry(PSIF_AO_TPDM, toc, (char *)tpdm_buffer_, buffer_size * sizeof(double));
        ++pq_pair_count;
        size_t index = 0;

        std::shared_ptr<SO_RS_Iterator> RSIter =
            std::make_shared<SO_RS_Iterator>(p, q, sobasis_, sobasis_, sobasis_, sobasis_);
        for (RSIter->first(); RSIter->is_done() == false; RSIter->next()) {
            int ish = RSIter->p();
            int jsh = RSIter->q();
            int ksh = RSIter->r();
            int lsh = RSIter->s();

            int n1 = sobasis_->nfunction(ish);
            int n2 = sobasis_->nfunction(jsh);
            int n3 = sobasis_->nfunction(ksh);
            int n4 = sobasis_->nfunction(lsh);

            // The starting orbital for each irrep can be grabbed from DPD
            int *sym_offsets = D->params->poff;

            for (int itr = 0; itr < n1; itr++) {
                int ifunc = sobasis_->function(ish) + itr;
                int isym = sobasis_->irrep(ifunc);
                int irel = sobasis_->function_within_irrep(ifunc);
                int iabs = sym_offsets[isym] + irel;
                for (int jtr = 0; jtr < n2; jtr++) {
                    int jfunc = sobasis_->function(jsh) + jtr;
                    int jsym = sobasis_->irrep(jfunc);
                    int jrel = sobasis_->function_within_irrep(jfunc);
                    int jabs = sym_offsets[jsym] + jrel;
                    for (int ktr = 0; ktr < n3; ktr++) {
                        int kfunc = sobasis_->function(ksh) + ktr;
                        int ksym = sobasis_->irrep(kfunc);
                        int krel = sobasis_->function_within_irrep(kfunc);
                        int kabs = sym_offsets[ksym] + krel;
                        for (int ltr = 0; ltr < n4; ltr++) {
                            int lfunc = sobasis_->function(lsh) + ltr;
                            int lsym = sobasis_->irrep(lfunc);
                            if (isym ^ jsym ^ ksym ^ lsym) continue;  // Not totally symmetric
                            int lrel = sobasis_->function_within_irrep(lfunc);
                            int labs = sym_offsets[lsym] + lrel;
                            int iiabs = iabs;
                            int jjabs = jabs;
                            int kkabs = kabs;
                            int llabs = labs;

                            int iiirrep = isym;
                            int jjirrep = jsym;
                            int kkirrep = ksym;
                            int llirrep = lsym;

                            int iirel = irel;
                            int jjrel = jrel;
                            int kkrel = krel;
                            int llrel = lrel;

                            if (ish == jsh) {
                                if (iabs < jabs) continue;

                                if (ksh == lsh) {
                                    if (kabs < labs) continue;
                                    if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                        if (ish == ksh)  // IIII case
                                            continue;
                                        else {  // IIJJ case
                                            SWAP_INDEX(ii, kk);
                                            SWAP_INDEX(jj, ll);
                                        }
                                    }
                                } else {  // IIJK case
                                    if (labs > kabs) {
                                        SWAP_INDEX(kk, ll);
                                    }
                                    if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                        SWAP_INDEX(ii, kk);
                                        SWAP_INDEX(jj, ll);
                                    }
                                }
                            } else {
                                if (ksh == lsh) {  // IJKK case
                                    if (kabs < labs) continue;
                                    if (iabs < jabs) {
                                        SWAP_INDEX(ii, jj);
                                    }
                                    if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                        SWAP_INDEX(ii, kk);
                                        SWAP_INDEX(jj, ll);
                                    }
                                } else {  // IJIJ case
                                    if (ish == ksh && jsh == lsh && INDEX2(iabs, jabs) < INDEX2(kabs, labs)) continue;
                                    // IJKL case
                                    if (iabs < jabs) {
                                        SWAP_INDEX(ii, jj);
                                    }
                                    if (kabs < labs) {
                                        SWAP_INDEX(kk, ll);
                                    }
                                    if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                        SWAP_INDEX(ii, kk);
                                        SWAP_INDEX(jj, ll);
                                    }
                                }
                            }

                            int ijsym = iiirrep ^ jjirrep;
                            size_t ijrow = D->params->rowidx[iiabs][jjabs];
                            size_t ijcol = D->params->colidx[iiabs][jjabs];
                            size_t klrow = D->params->rowidx[kkabs][llabs];
                            size_t klcol = D->params->colidx[kkabs][llabs];
                            // We know that ijkl is totally symmetric, so klsym
                            // must be the same as ijsym
                            if ((ijsym == irrep) && (ijrow >= first_row) && (ijrow < last_row)) {
                                tpdm_buffer_[index] += 0.5 * D->matrix[ijsym][ijrow - first_row][klcol];
                            }
                            if ((ijsym == irrep) && (klrow >= first_row) && (klrow < last_row)) {
                                tpdm_buffer_[index] += 0.5 * D->matrix[ijsym][klrow - first_row][ijcol];
                            }
                            ++index;
                        }
                    }
                }
            }
        }  // End rs iterator
        psio_->write_entry(PSIF_AO_TPDM, toc, (char *)tpdm_buffer_, buffer_size * sizeof(double));
    }  // End pq iterator
    delete[] toc;
}
}  // namespace psi
