/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libiwl/iwl.hpp"

using namespace psi;

/**
 * Transforms the kinetic energy plus nuclear attraction one-electron
 * integrals.  This function is currently limited to IWL input and output and
 * Pitzer ordering, regardless of how the parameters passed to the constructor.
 *
 * @param s1 - the MOSpace for the bra
 * @param s2 - the MOSpace for the ket
 *
 * N.B. This need not be called if a two-electron transformation is performed, as
 * the sort_so_tei routine will perform this transformation in addition to the
 * Fock matrix construction.
 */
void IntegralTransform::transform_T_plus_V(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2) {
    check_initialized();
    std::vector<double> soInts(nTriSo_), T(nTriSo_);
    if (print_ > 4) outfile->Printf("The SO basis kinetic energy integrals\n");
    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_T, T.data(), nTriSo_, 0, print_ > 4, "outfile");
    if (print_ > 4) outfile->Printf("The SO basis nuclear attraction integrals\n");
    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_V, soInts.data(), nTriSo_, 0, print_ > 4, "outfile");
    // Add the nuclear and kinetic energy integrals
    std::transform(soInts.begin(), soInts.end(), T.begin(), soInts.begin(), std::plus<double>());
    if (transformationType_ == TransformationType::Restricted) {
        transform_oei_restricted(s1, s2, soInts, PSIF_MO_OEI);
    } else {
        transform_oei_unrestricted(s1, s2, soInts, PSIF_MO_A_OEI, PSIF_MO_B_OEI);
    }
}

/**
 * Transforms the one-electron integrals.  This function is currently limited to
 * IWL input and output and Pitzer ordering, regardless of how the parameters passed
 * to the constructor.
 *
 * @param s1 - the MOSpace for the bra
 * @param s2 - the MOSpace for the ket
 * @param labels - array with the labels of the OEI to be transformed
 *
 * @note The labels array is assumed to have the labels for the OEI in this
 * order:
 * 0. labels[0] is the AO integrals label,
 * 1. labels[1] is the MO integrals label, for the restricted case
 * 2. labels[2] is the alpha MO integrals label, for the unrestricted case
 * 3. labels[3] is the beta MO integrals label, for the unrestricted case
 */
void IntegralTransform::transform_oei(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2,
                                      std::array<std::string, 4> labels) {
    check_initialized();
    std::vector<double> soInts(nTriSo_);
    if (print_ > 4) outfile->Printf("Grabbing " + labels[0] + "\n");
    IWL::read_one(psio_.get(), PSIF_OEI, labels[0].c_str(), soInts.data(), nTriSo_, 0, print_ > 4, "outfile");
    if (transformationType_ == TransformationType::Restricted) {
        transform_oei_restricted(s1, s2, soInts, labels[1].c_str());
    } else {
        transform_oei_unrestricted(s1, s2, soInts, labels[2].c_str(), labels[3].c_str());
    }
}

/**
 * Perform restricted transformation of the one-electron integrals. This function is currently limited to
 * IWL input and output and Pitzer ordering, regardless of how the parameters passed
 * to the constructor.
 *
 * @param s1 - the MOSpace for the bra
 * @param s2 - the MOSpace for the ket
 * @param soInts - buffer with the SO-basis integrals to be transformed
 * @param label - MO label
 */
void IntegralTransform::transform_oei_restricted(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2,
                                                 const std::vector<double> &soInts, std::string label) {
    std::vector<double> moInts(nTriMo_);
    std::vector<int> order(nmo_);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    std::iota(order.begin(), order.end(), 0);
    for (int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h) {
        double **pCa = Ca_->pointer(h);
        trans_one(sopi_[h], mopi_[h], const_cast<double *>(soInts.data()), moInts.data(), pCa, soOffset,
                  &(order[moOffset]));
        soOffset += sopi_[h];
        moOffset += mopi_[h];
    }
    if (print_ > 4) {
        outfile->Printf("The MO basis " + label + "\n");
        print_array(moInts.data(), nmo_, "outfile");
    }
    IWL::write_one(psio_.get(), PSIF_OEI, label.c_str(), nTriMo_, moInts.data());
}

void IntegralTransform::transform_oei_unrestricted(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2,
                                                   const std::vector<double> &soInts, std::string A_label,
                                                   std::string B_label) {
    std::vector<double> moInts(nTriMo_);
    std::vector<int> order(nmo_);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    std::iota(order.begin(), order.end(), 0);
    for (int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h) {
        double **pCa = Ca_->pointer(h);
        trans_one(sopi_[h], mopi_[h], const_cast<double *>(soInts.data()), moInts.data(), pCa, soOffset,
                  &(order[moOffset]));
        soOffset += sopi_[h];
        moOffset += mopi_[h];
    }
    if (print_ > 4) {
        outfile->Printf("The MO basis alpha " + A_label + "\n");
        print_array(moInts.data(), nmo_, "outfile");
    }
    IWL::write_one(psio_.get(), PSIF_OEI, A_label.c_str(), nTriMo_, moInts.data());

    for (int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h) {
        double **pCb = Cb_->pointer(h);
        trans_one(sopi_[h], mopi_[h], const_cast<double *>(soInts.data()), moInts.data(), pCb, soOffset,
                  &(order[moOffset]));
        soOffset += sopi_[h];
        moOffset += mopi_[h];
    }
    if (print_ > 4) {
        outfile->Printf("The MO basis beta " + B_label + "\n");
        print_array(moInts.data(), nmo_, "outfile");
    }
    IWL::write_one(psio_.get(), PSIF_OEI, B_label.c_str(), nTriMo_, moInts.data());
}

/**
 * Transforms a packed symmetric matrix.
 *
 * @param m - input matrix row dimension
 * @param n - output matrix row dimension
 * @param input - pointer to input integrals (the lower-triangle of a symmetric matrix)
 * @param pointer to output integrals (the lower-triangle of a symmetric matrix)
 * @param C transformation matrix (rectangular, m X n)
 * @param offset - the point in the full list of SOs where we want to start.  This is
 *                 useful for transforming integrals one irrep at a time and in this
 *                 case the offset would correspond to the position of the first
 *                 orbital in the current irrep.
 * @param order - a reordering array to change the order of the output
 * @param backtransform - whether this is a forward or backwards transformation
 * @param scale - the amount of the existing output buffer to mix into the result
 */

void IntegralTransform::trans_one(int m, int n, double *input, double *output, double **C, int offset, int *order,
                                  bool backtransform, double scale) {
    // TODO the order argument is actually not used right now.  I don't know that anybody will need it
    // so I haven't bothered so far...
    int dim = (m > n) ? m : n;
    double **TMP0 = block_matrix(dim, dim);
    double **TMP1 = block_matrix(dim, dim);

    for (int p = 0; p < m; ++p) {
        for (int q = 0; q <= p; ++q) {
            size_t pq = INDEX((p + offset), (q + offset));
            TMP0[p][q] = TMP0[q][p] = input[pq];
        }
    }
    int nc;
    if (backtransform) {
        nc = m;
        if (m && n) {
            C_DGEMM('n', 't', m, n, m, 1.0, TMP0[0], dim, C[0], nc, 0.0, TMP1[0], dim);
            C_DGEMM('n', 'n', n, n, m, 1.0, C[0], nc, TMP1[0], dim, 0.0, TMP0[0], dim);
        }
    } else {
        nc = n;
        if (m && n) {
            C_DGEMM('n', 'n', m, n, m, 1.0, TMP0[0], dim, C[0], nc, 0.0, TMP1[0], dim);
            C_DGEMM('t', 'n', n, n, m, 1.0, C[0], nc, TMP1[0], dim, 0.0, TMP0[0], dim);
        }
    }

    for (int p = 0; p < nc; ++p) {
        for (int q = 0; q <= p; ++q) {
            size_t P = order[p];
            size_t Q = order[q];
            size_t PQ = INDEX(P, Q);
            output[PQ] = scale * output[PQ] + TMP0[p][q];
        }
    }

    free_block(TMP0);
    free_block(TMP1);

    return;
}

/**
 * Generates the frozen core operator, Fock matrix and one electron integrals
 * in the MO basis by looping over the IWL SO integral file on disk.  The resulting
 * integrals are written to PSIF_SO_OEI and are labelled according to the macros in
 * psifiles.h.  Regardless of any parameters, all integrals are transformed and
 * only IWL format is used to store the results.
 */
void IntegralTransform::generate_oei() {
#if 0
    // Set aside some memory for the frozen core density and frozen core operator
    double *aFzcD  = init_array(nTriSo_);
    double *aFzcOp = init_array(nTriSo_);
    double *aD     = init_array(nTriSo_);
    double *aFock  = init_array(nTriSo_);
    double *aoH    = init_array(nTriSo_);
    double *bFzcD  = aFzcD;
    double *bFzcOp = aFzcOp;
    double *bD     = aD;
    double *bFock  = aFock;
    if(transformationType_ != Restricted){
        bFzcD  = init_array(nTriSo_);
        bFzcOp = init_array(nTriSo_);
        bD     = init_array(nTriSo_);
        bFock  = init_array(nTriSo_);
    }

    // Form the Density matrices
    for(int h = 0, soOffset = 0; h < nirreps_; ++h){
        double **pCa = Ca_->pointer(h);
        double **pCb = Cb_->pointer(h);
        for(int p = 0; p < sopi_[h]; ++p){
            for(int q = 0; q <= p; ++q){
                int pq = INDEX((p + soOffset), (q + soOffset));
                for(int i = 0; i < frzcpi_[h]; ++i)
                    aFzcD[pq] += pCa[p][i] * pCa[q][i];
                for(int i = 0; i < clsdpi_[h] + openpi_[h]; ++i)
                    aD[pq] += pCa[p][i] * pCa[q][i];
                if(transformationType_ != Restricted){
                    for(int i = 0; i < frzcpi_[h]; ++i)
                        bFzcD[pq] += pCb[p][i] * pCb[q][i];
                    for(int i = 0; i < clsdpi_[h]; ++i)
                        bD[pq] += pCb[p][i] * pCb[q][i];
                }
            }
        }
        soOffset += sopi_[h];
    }

    double *T = init_array(nTriSo_);
    if(print_>4) outfile->Printf( "The SO basis kinetic energy integrals\n");
    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_T,   T, nTriSo_, 0, print_ > 4, outfile);
    if(print_>4) outfile->Printf( "The SO basis nuclear attraction integrals\n");
    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_V, aoH, nTriSo_, 0, print_ > 4, outfile);

    for(int pq=0; pq < nTriSo_; ++pq){
        aoH[pq] += T[pq];
        aFzcOp[pq] = aoH[pq];
        aFock[pq]  = aoH[pq];
        if(transformationType_ != Restricted){
            bFock[pq]  = aoH[pq];
            bFzcOp[pq] = aoH[pq];
        }
    }
    free(T);

    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    int soIntFile = PSIF_SO_TEI;

    IWL *iwl = new IWL(psio_.get(), soIntFile, tolerance_, 1, 1);
    Label *lblptr = iwl->labels();
    Value *valptr = iwl->values();
    int lastbuf   = iwl->last_buffer();
    for(int index = iwl->index(); index < iwl->buffer_count(); ++index){
        int labelIndex = 4*index;
        int p = std::abs((int) lblptr[labelIndex++]);
        int q = (int) lblptr[labelIndex++];
        int r = (int) lblptr[labelIndex++];
        int s = (int) lblptr[labelIndex++];
        double value = (double) valptr[index];
        build_fzc_and_fock(p, q, r, s, value, aFzcD, bFzcD,
                           aFzcOp, bFzcOp, aD, bD, aFock, bFock);
    } /* end loop through current buffer */

    /* Now run through the rest of the buffers in the file */
    while(!lastbuf){
        iwl->fetch();
        lastbuf = iwl->last_buffer();
        for(int index = iwl->index(); index < iwl->buffer_count(); ++index){
            int labelIndex = 4*index;
            int p = std::abs((int) lblptr[labelIndex++]);
            int q = (int) lblptr[labelIndex++];
            int r = (int) lblptr[labelIndex++];
            int s = (int) lblptr[labelIndex++];
            double value = (double) valptr[index];
            build_fzc_and_fock(p, q, r, s, value, aFzcD, bFzcD,
                               aFzcOp, bFzcOp, aD, bD, aFock, bFock);
        } /* end loop through current buffer */
    } /* end loop over reading buffers */
    iwl->set_keep_flag(true);
    delete iwl;

    double *moInts = init_array(nTriMo_);
    int *order = init_int_array(nmo_);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    for(int n = 0; n < nmo_; ++n) order[n] = n;
    if(print_)
        outfile->Printf( "\tTransforming the one-electron integrals and constructing Fock matrices\n");
    if(transformationType_ == Restricted){
        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aoH, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis one-electron integrals\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_OEI, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aFzcOp, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis frozen core operator\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_FZC, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aFock, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis Fock operator\n");
            print_array(moInts, nmo_, outfile);
        }

        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_FOCK, nTriMo_, aFock);
    }else{
        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aoH, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis alpha one-electron integrals\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_A_OEI, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCb = Cb_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aoH, moInts, pCb, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis beta one-electron integrals\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_B_OEI, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aFzcOp, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis alpha frozen core operator\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_A_FZC, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCb = Cb_->pointer(h);
            trans_one(sopi_[h], mopi_[h], bFzcOp, moInts, pCb, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis beta frozen core operator\n");
            print_array(moInts, nmo_, outfile);
        }
        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aFock, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis alpha Fock operator\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_A_FOCK, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCb = Cb_->pointer(h);
            trans_one(sopi_[h], mopi_[h], bFock, moInts, pCb, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            outfile->Printf( "The MO basis beta Fock operator\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_B_FOCK, nTriMo_, moInts);
    }
    free(order);
    free(moInts);
    free(aFzcD);
    free(aFzcOp);
    free(aD);
    free(aoH);
    free(aFock);
    if(transformationType_ != Restricted){
        free(bFzcD);
        free(bFzcOp);
        free(bD);
        free(bFock);
    }

    dpd_set_default(currentActiveDPD);
#endif
}
