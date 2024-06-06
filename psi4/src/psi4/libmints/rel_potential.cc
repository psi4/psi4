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

#include "psi4/libmints/rel_potential.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <libint2/engine.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define RELVDEBUG 0

using namespace psi;

RelPotentialInt::RelPotentialInt(std::vector<SphericalTransform>& st, std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> bs2, int deriv)
    : OneBodyAOInt(st, bs1, bs2, deriv) {
    if (deriv > 0) {
        throw PSIEXCEPTION("RelPotentialInt: deriv > 0 is not supported.");
    }

    if (bs1 != bs2) {
        outfile->Printf("*********************************************************************************************************************\n");
        outfile->Printf("When computing potential integrals with different bra and ket basis, the atom definition is taken from the bra basis.\n");
        outfile->Printf("*********************************************************************************************************************\n");
    }

    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    // Setup the initial field of partial charges
    std::vector<std::pair<double, std::array<double, 3>>> params;
    for (int A = 0; A < bs1_->molecule()->natom(); A++) {
        params.push_back({
            (double)bs1_->molecule()->Z(A),
            {bs1_->molecule()->x(A), bs1_->molecule()->y(A), bs1_->molecule()->z(A)}});
    }

    set_chunks(4);
    engine0_ = std::make_unique<libint2::Engine>(libint2::Operator::opVop, max_nprim, max_am, 0);
    engine0_->set_params(params);
    // If you want derivatives of these integrals, just take the code from potential.cc's constructor
    // and change out the operator-type. We don't have these integrals now for want of a use case.

    buffer_ = nullptr;
    buffers_.resize(nchunk_);
}

void RelPotentialInt::set_charge_field(const std::vector<std::pair<double, std::array<double, 3>>>& Zxyz) {
    engine0_->set_params(Zxyz);
    if (engine1_) engine1_->set_params(Zxyz);
    if (engine2_) engine2_->set_params(Zxyz);
    Zxyz_ = Zxyz;
}

RelPotentialSOInt::RelPotentialSOInt(const std::shared_ptr<OneBodyAOInt>& aoint,
                                     const std::shared_ptr<IntegralFactory>& fact)
    : OneBodySOInt(aoint, fact) {
    natom_ = ob_->basis1()->molecule()->natom();
}

RelPotentialSOInt::RelPotentialSOInt(const std::shared_ptr<OneBodyAOInt>& aoint, const IntegralFactory* fact)
    : OneBodySOInt(aoint, fact) {
    natom_ = ob_->basis1()->molecule()->natom();
}
