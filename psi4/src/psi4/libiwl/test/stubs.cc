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

// Minimal hosting globals for the standalone IWL round-trip test executable.
// libpsio and libpsi4util are not standalone static archives -- they reference
// runtime globals normally defined in core.cc (the python bindings TU) and a
// member function on libdpd's DPD class. This file supplies just enough to
// satisfy the link without dragging in libdpd or core.cc.

#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>

#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {

// libciomr / libpsio scratch-name globals; normally set by psi4 startup.
std::string restart_id;
char *psi_file_prefix = strdup("iwl_test");

// libpsi4util's `outfile` is reached only from print paths the test never
// triggers (printflag = 0 everywhere). The symbol still has to exist.
std::shared_ptr<PsiOutStream> outfile;

// libpsio's PSIO::close calls global_dpd_->file4_cache_del_filenum during
// every file close. In a non-DPD context the cache walk is a no-op anyway,
// but dereferencing a null global_dpd_ is UB. Provide a minimal local DPD
// class with the same name and the one referenced method as a no-op, plus a
// non-null instance to dispatch through. We link without libdpd so there is
// no ODR clash with the real DPD definition.
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
