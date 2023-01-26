/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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
#ifdef USING_gdma

#include "psi4/psi4-dec.h"
#include "psi4/libmints/matrix.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include <iostream>

extern "C" {
void run_gdma(const char* outfilename, const char* datfilename);
int get_nsites();
int get_order(int site);
double get_dma_value(int site, int addr);
double get_tot_value(int addr);
};

namespace psi {
namespace gdma_interface {

SharedWavefunction gdma_interface(SharedWavefunction ref_wfn, Options& options, const std::string& datfilename) {
    run_gdma(outfile_name.c_str(), datfilename.c_str());
    // Reopen the outfile
    if (outfile_name == "stdout") {
        outfile = std::make_shared<PsiOutStream>();
    } else {
        outfile = std::make_shared<PsiOutStream>(outfile_name, std::ostream::app);
    }
    int nsites = get_nsites();
    int maxorder = 0;
    for (int site = 1; site <= nsites; ++site) {
        int site_order = get_order(site);
        maxorder = maxorder > site_order ? maxorder : site_order;
    }
    int nvals = (maxorder + 1) * (maxorder + 1);
    auto dmavals = std::make_shared<Matrix>("Spherical Harmonic DMA for each site", nsites, nvals);
    for (int site = 1; site <= nsites; ++site) {
        int site_order = get_order(site);
        int site_nvals = (site_order + 1) * (site_order + 1);
        for (int n = 1; n <= site_nvals; ++n) {
            double val = get_dma_value(site, n);
            dmavals->set(site - 1, n - 1, val);
        }
    }
    auto totvals = std::make_shared<Matrix>("Total multipoles, translated to the origin", 1, nvals);
    for (int n = 1; n <= nvals; ++n) {
        double val = get_tot_value(n);
        totvals->set(0, n - 1, val);
    }

    Process::environment.arrays["DMA DISTRIBUTED MULTIPOLES"] = dmavals;
    Process::environment.arrays["DMA TOTAL MULTIPOLES"] = totvals;
    outfile->Printf(
        "\n  DMA results are available in the Python driver through the\n"
        "\t  variable('DMA DISTRIBUTED MULTIPOLES')\n"
        "  and\n"
        "\t  variable('DMA TOTAL MULTIPOLES')\n"
        "  commands.\n\n");

    return ref_wfn;
}

}  // namespace gdma_interface
}  // namespace psi
#endif
