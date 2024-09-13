/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "psi4/libdpd/dpd.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace cctransort {

void memcheck(const int reference) {
    size_t irrep_size, size;
    dpdbuf4 Z;

    outfile->Printf("\n");

    if (reference == 0) {
        global_dpd_->buf4_init(&Z, 99, 0, 5, 5, 5, 5, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of <ab|cd> integrals: %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));

        global_dpd_->buf4_init(&Z, 99, 0, 10, 5, 10, 5, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of <ia|bc> integrals: %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));

        global_dpd_->buf4_init(&Z, 99, 0, 0, 5, 0, 5, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of tijab amplitudes:  %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));

    } else if (reference == 1) {
        global_dpd_->buf4_init(&Z, 99, 0, 5, 5, 5, 5, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of <ab|cd> integrals: %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));

        global_dpd_->buf4_init(&Z, 99, 0, 10, 5, 10, 5, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of <ia|bc> integrals: %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));

        global_dpd_->buf4_init(&Z, 99, 0, 0, 5, 0, 5, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of tIjAb amplitudes:  %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));

    } else if (reference == 2) {
        global_dpd_->buf4_init(&Z, 99, 0, 7, 7, 7, 7, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of <AB|CD> integrals: %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));
        global_dpd_->buf4_init(&Z, 99, 0, 17, 17, 17, 17, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of <ab|cd> integrals: %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));
        global_dpd_->buf4_init(&Z, 99, 0, 28, 28, 28, 28, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of <Ab|Cd> integrals: %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));

        global_dpd_->buf4_init(&Z, 99, 0, 20, 5, 20, 5, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of <IA|BC> integrals: %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));

        global_dpd_->buf4_init(&Z, 99, 0, 30, 15, 30, 15, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of <ia|bc> integrals: %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));

        global_dpd_->buf4_init(&Z, 99, 0, 24, 28, 24, 28, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of <Ia|Bc> integrals: %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));

        global_dpd_->buf4_init(&Z, 99, 0, 27, 29, 27, 29, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of <iA|bC> integrals: %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));

        global_dpd_->buf4_init(&Z, 99, 0, 22, 28, 22, 28, 0, "Just a template");
        size = 0;
        for (int h = 0; h < Z.params->nirreps; h++) {
            irrep_size = (size_t)Z.params->rowtot[h] * Z.params->coltot[h];
            size += irrep_size;
            outfile->Printf("\tSize of irrep %d of tIjAb amplitudes:  %10.3lf (MW) / %10.3lf (MB)\n", h,
                            irrep_size / 1e6, (irrep_size / 1e6) * sizeof(double));
        }
        global_dpd_->buf4_close(&Z);
        outfile->Printf("\tTotal:                                %10.3lf (MW) / %10.3lf (MB)\n\n", size / 1e6,
                        (size / 1e6) * sizeof(double));
    }
}

}  // namespace cctransort
}  // namespace psi
