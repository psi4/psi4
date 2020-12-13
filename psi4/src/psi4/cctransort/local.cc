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

/*! \file
    \ingroup CCTRANSORT
    \brief Enter brief description of file here
*/
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/local.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libdpd/dpd.h"

namespace psi {
namespace cctransort {

void localize_occupied(std::shared_ptr<Wavefunction> wfn) {
    auto docc = wfn->doccpi()[0];
    auto nfrz = wfn->frzcpi()[0];
    auto nso = wfn->nso();
    auto nmo = wfn->nmo();

    auto C_occ = std::make_shared<Matrix>(nso, docc-nfrz);
    
    std::shared_ptr<Matrix> C_full = wfn->Ca();;

    for(int i=0; i < nso; ++i) {
        for(int j=nfrz; j < docc; ++j)
            C_occ->set(i,j-nfrz, C_full->get(i,j));
    }

    outfile->Printf("Orbitals without localizing:\n");
    C_occ->print();
    // Build localizer obj using active occupied orbitals
    std::shared_ptr<Localizer> local_obj = Localizer::build("PIPEK_MEZEY", wfn->basisset(), C_occ);
    local_obj->localize();

    auto new_C_occ = std::make_shared<Matrix>(nso, docc-nfrz);
    new_C_occ->copy(local_obj->L()); 

    for(int i=0; i < nso; ++i) {
        for(int j=nfrz; j < docc; ++j)
            C_full->set(i,j, new_C_occ->get(i,j-nfrz));
    }

    outfile->Printf("Orbitals after localizing:\n");
    new_C_occ->print();

    wfn->Ca()->copy(C_full);
    wfn->Cb()->copy(C_full);

}

} //namespace cctransort
}  // namespace psi
