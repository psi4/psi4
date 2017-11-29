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

#include "fcidump_helper.h"

#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <string>

#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"

namespace psi {
namespace fcidump {
void fcidump_tei_helper(int nirrep, bool restricted, std::map<std::string, int> DPD_info, double ints_tolerance,
                        std::string fname) {
    outfile->Printf("Writing TEI integrals in FCIDUMP format to " + fname + "\n");
    // Append to the file created by the fcidump function Python-side
    FILE* intdump = fopen(fname.c_str(), "a");

    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(DPD_info["instance_id"]);

    _default_psio_lib_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    dpdbuf4 K;

    // RHF
    if (restricted) {
        /* Convert a molecular orbital index [0,1,...] to [1,2,...] (i.e. from zero-based to one-based). */
        auto mo_index = [](const int i) { return i + 1; };
        // We want only the permutationally unique integrals, see libtrans documentation for details
        // DPD_info["alpha_MO"] is DPD_ID("[A>=A]+")
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, DPD_info["alpha_MO"], DPD_info["alpha_MO"],
                               DPD_info["alpha_MO"], DPD_info["alpha_MO"], 0, "MO Ints (AA|AA)");
        detail::write_tei_to_disk(intdump, nirrep, K, ints_tolerance, mo_index, mo_index);
        global_dpd_->buf4_close(&K);

        // if (dump_dipoles) {
        //     // TODO: dipoles filename.
        //     // TODO: use IntegralTransform.  But, as I don't know how that
        //     // works, for now follow preppert.cc (part of ccresponse) and do
        //     // the transformation ourselves.
        //     std::string fname[] = {"DIPOLES_X", "DIPOLES_Y", "DIPOLES_Z"};
        //     MintsHelper mints(wfn->basisset(), Process::environment.options, 0);
        //     std::vector<std::shared_ptr<Matrix>> dipole = mints.so_dipole();
        //     FILE* dipoledump;
        //     double frz_contrib;
        //     Vector3 origin =
        //         Vector3(0, 0, 0);  // In serious trouble if being asked for properties after moving the molecule...
        //     SharedVector ndip = DipoleInt::nuclear_contribution(molecule, origin);
        //     for (int i = 0; i < 3; i++) {
        //         dipoledump = fopen(fname[i].c_str(), "w");
        //         write_oei_prop_to_disk(dipoledump, wfn, dipole[i], ints_tolerance, mo_index, &frz_contrib);
        //         fprintf(dipoledump, "%29.20E%4d%4d\n", ndip->get(i) + frz_contrib, 0, 0);
        //         fclose(dipoledump);
        //     }
        //     // BONUS: and quadrupole moments.  (Just zz for now.)
        //     std::vector<std::shared_ptr<Matrix>> trquad = mints.so_traceless_quadrupole();
        //     SharedVector nquad = QuadrupoleInt::nuclear_contribution(molecule, origin);
        //     {
        //         int ij = 5;
        //         dipoledump = fopen("TRQUAD_ZZ", "w");
        //         write_oei_prop_to_disk(dipoledump, wfn, trquad[ij], ints_tolerance, mo_index, &frz_contrib);
        //         fprintf(dipoledump, "%29.20E%4d%4d\n", nquad->get(0, ij) + frz_contrib, 0, 0);
        //         fclose(dipoledump);
        //     }
        // }
    } else {
        /* Convert an alpha spin-orbital index [0,1,...] to [1,3,...] (i.e. from
         * zero-based to one-based, with corresponding beta orbitals interwoven).
         */
        auto alpha_index = [](const int i) { return 2 * i + 1; };
        /* Convert a beta spin-orbital index [0,1,...] to [2,4,...] (i.e. from
         * zero-based to one-based, with corresponding alpha orbitals interwoven).
         */
        auto beta_index = [](const int i) { return 2 * (i + 1); };
        // We want only the permutationally unique integrals, see libtrans documentation for details
        // DPD_info["alpha_MO"] is DPD_ID("[A>=A]+"), while DPD_info["beta_MO"] is DPD_ID("[a>=a]+")
        // alpha-alpha
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, DPD_info["alpha_MO"], DPD_info["alpha_MO"],
                               DPD_info["alpha_MO"], DPD_info["alpha_MO"], 0, "MO Ints (AA|AA)");
        detail::write_tei_to_disk(intdump, nirrep, K, ints_tolerance, alpha_index, alpha_index);
        global_dpd_->buf4_close(&K);
        // beta-beta
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, DPD_info["beta_MO"], DPD_info["beta_MO"], DPD_info["beta_MO"],
                               DPD_info["beta_MO"], 0, "MO Ints (aa|aa)");
        detail::write_tei_to_disk(intdump, nirrep, K, ints_tolerance, beta_index, beta_index);
        global_dpd_->buf4_close(&K);
        // alpha-beta
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, DPD_info["alpha_MO"], DPD_info["beta_MO"],
                               DPD_info["alpha_MO"], DPD_info["beta_MO"], 0, "MO Ints (AA|aa)");
        detail::write_tei_to_disk(intdump, nirrep, K, ints_tolerance, alpha_index, beta_index);
        global_dpd_->buf4_close(&K);
    }
    _default_psio_lib_->close(PSIF_LIBTRANS_DPD, 1);

    fclose(intdump);
}

namespace detail {
void write_tei_to_disk(FILE* intdump, int nirrep, dpdbuf4& K, double ints_tolerance, OrbitalIndexing indx1,
                       OrbitalIndexing indx2) {
    for (int h = 0; h < nirrep; ++h) {
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for (int pq = 0; pq < K.params->rowtot[h]; ++pq) {
            int p = K.params->roworb[h][pq][0];
            int q = K.params->roworb[h][pq][1];
            for (int rs = 0; rs < K.params->coltot[h]; ++rs) {
                int r = K.params->colorb[h][rs][0];
                int s = K.params->colorb[h][rs][1];

                if (std::abs(K.matrix[h][pq][rs]) > ints_tolerance)
                    fprintf(intdump, "%28.20E%4d%4d%4d%4d\n", K.matrix[h][pq][rs], indx1(p), indx1(q), indx2(r),
                            indx2(s));
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
}

void write_oei_prop_to_disk(FILE* intdump, std::shared_ptr<Wavefunction> wfn, std::shared_ptr<Matrix> prop_ints,
                            double ints_tolerance, OrbitalIndexing indx, double* frz_contrib) {
    double** scf = wfn->Ca()->to_block_matrix();  // TODO: UHF
    int nso = wfn->nso();
    int nmo = wfn->nmo();
    Dimension frzcpi = wfn->frzcpi();
    Dimension active_mopi = wfn->nmopi() - frzcpi - wfn->frzvpi();
    int nirrep = wfn->nirrep();

    double** TMP1 = prop_ints->to_block_matrix();
    double** TMP2 = block_matrix(nso, nso);

    C_DGEMM('n', 'n', nso, nmo, nso, 1, TMP1[0], nso, scf[0], nmo, 0, TMP2[0], nso);
    C_DGEMM('t', 'n', nmo, nmo, nso, 1, scf[0], nmo, TMP2[0], nso, 0, TMP1[0], nmo);
    // TMP1 now holds the dipole integrals in the MO basis, ordered 1->nmo (in symmetry blocks).
    // We just want to print out the active orbitals...
    // Can't just loop over the two indices as we only know the
    // active orbitals per irrep.  Instead, loop over everything
    // and just print out non-zero integrals (bit slower as we
    // don't use symmetry, but this isn't a hotspot...)
    int ioff1 = 0;
    int nfrz1 = 0;
    for (int h1 = 0; h1 < nirrep; ++h1) {
        nfrz1 += frzcpi[h1];
        int ioff2 = ioff1;
        int nfrz2 = nfrz1;
        for (int h2 = h1; h2 < nirrep; ++h2) {
            for (int m1 = frzcpi[h1]; m1 < frzcpi[h1] + active_mopi[h1]; ++m1) {
                int m2_init = h1 == h2 ? m1 : frzcpi[h2];
                for (int m2 = m2_init; m2 < frzcpi[h2] + active_mopi[h2]; ++m2) {
                    int iorb1 = m1 + ioff1;
                    int iorb2 = m2 + ioff2;
                    double intgrl = TMP1[iorb1][iorb2];
                    if (std::abs(intgrl) > ints_tolerance)
                        fprintf(intdump, "%29.20E%4d%4d\n", intgrl, indx(iorb1 - nfrz1), indx(iorb2 - nfrz2));
                }
            }
            nfrz2 += frzcpi[h2];
            ioff2 += prop_ints->rowdim(h2);
        }
        ioff1 += prop_ints->rowdim(h1);
    }
    // The contribution of the frozen core orbitals to a one-body
    // expectation value is just \sum_i <i|O_1|i>.
    *frz_contrib = 0.0;
    ioff1 = 0;
    for (int h = 0; h < nirrep; ++h) {
        for (int m = 0; m < frzcpi[h]; ++m) {
            int iorb = m + ioff1;
            *frz_contrib += 2 * TMP1[iorb][iorb];  // 2* for RHF.
        }
        ioff1 += prop_ints->rowdim(h);
    }
}
}  // End namespace detail
}  // End namespace fcidump
}  // End namespace psi
