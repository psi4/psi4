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

// Minimal hosting globals for the standalone libpsio round-trip executable.
// libpsio/libpsi4util/libciomr are not standalone archives -- they reference a
// few runtime globals normally defined in core.cc. Supply just enough to link
// without dragging in libdpd or the python bindings TU.

#include <cstring>
#include <memory>
#include <string>

#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {

std::string restart_id;
char *psi_file_prefix = strdup("psio_test");
std::shared_ptr<PsiOutStream> outfile;

// libciomr's DSYEV_ascending lands in the same unity TU as helpers libpsio
// pulls in; the standalone link references its LAPACK wrapper though the test
// never diagonalizes. A no-op leaf keeps the test self-contained.
void C_DSYEV(char, char, int, double *, int, double *, double *, int) {}

// PSIO::close calls global_dpd_->file4_cache_del_filenum on every file close.
// In a non-DPD context the cache is empty, so a minimal local DPD with a no-op
// method and a non-null instance keeps the call well-defined. Linked without
// libdpd, so no ODR clash with the real definition.
class DPD {
   public:
    void file4_cache_del_filenum(unsigned long);
};
void DPD::file4_cache_del_filenum(unsigned long) {}

namespace {
DPD stub_dpd_instance;
}
DPD *global_dpd_ = &stub_dpd_instance;

}  // namespace psi
