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
#include "psi4/libmints/electricfield.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include <stdexcept>
#include <vector>
#include "psi4/libciomr/libciomr.h"

#include <libint2/engine.h>

using namespace psi;

ElectricFieldInt::ElectricFieldInt(std::vector<SphericalTransform>& spherical_transforms, std::shared_ptr<BasisSet> bs1,
                                   std::shared_ptr<BasisSet> bs2, int nderiv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv) {
    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    if (nderiv == 0) {
        engine0_ = std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_am, 1);
        set_chunks(9);
    } else {
        throw FeatureNotImplemented("LibMints", "ElectricFieldInts called with deriv > 0", __FILE__, __LINE__);
    }

    buffer_ = nullptr;
    buffers_.resize(nchunk_);
}

ElectricFieldInt::~ElectricFieldInt() { }

void ElectricFieldInt::set_origin(const Vector3 &origin) {
    set_chunks(9);
    engine0_->set_params(std::vector<std::pair<double, std::array<double,3>>>{{-1.0, {origin[0], origin[1], origin[2]}}});
}

void ElectricFieldInt::compute(std::vector<SharedMatrix> &result) {
    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();
    const auto bs1_equiv_bs2 = (bs1_ == bs2_);

    if (9 != (size_t)nchunk_) {
        outfile->Printf("there should be 9 chunks in ElectricFieldInt::compute().  You should call set_origin first.");
        throw SanityCheckError("OneBodyInt::compute(result): result incorrect length.", __FILE__, __LINE__);
    }

    // Check the individual matrices, we can only handle nirrep() == 1
    for (int chunk = 0; chunk < result.size(); ++chunk) {
        const auto a = result[chunk];
        if (a->nirrep() != 1) {
            throw SanityCheckError("OneBodyInt::compute(result): one or more of the matrices given has symmetry.",
                                   __FILE__, __LINE__);
        }
    }

    for (const auto &pair : shellpairs_) {
        int p1 = pair.first;
        int p2 = pair.second;

        const auto &s1 = bs1_->l2_shell(p1);
        const auto &s2 = bs2_->l2_shell(p2);
        int ni = bs1_->shell(p1).nfunction();
        int nj = bs2_->shell(p2).nfunction();
        int i_offset = bs1_->shell_to_basis_function(p1);
        int j_offset = bs2_->shell_to_basis_function(p2);

        // Compute the shell
        compute_pair(s1, s2);

        // For each integral that we got put in its contribution
        for (int r = 6; r < nchunk_; ++r) {
            const double *location = buffers_[r];
            for (int p = 0; p < ni; ++p) {
                for (int q = 0; q < nj; ++q) {
                    result[r-6]->add(0, i_offset + p, j_offset + q, *location);
                    if (bs1_equiv_bs2 && p1 != p2) {
                        result[r-6]->add(0, j_offset + q, i_offset + p, *location);
                    }
                    location++;
                }
            }
        }
    }
}

Vector3 ElectricFieldInt::nuclear_contribution(const Vector3& origin, std::shared_ptr<Molecule> mol) {
    int natom = mol->natom();

    double Ex = 0.0;
    double Ey = 0.0;
    double Ez = 0.0;
    for (int i = 0; i < natom; ++i) {
        double x = origin[0] - mol->x(i);
        double y = origin[1] - mol->y(i);
        double z = origin[2] - mol->z(i);
        double r2 = x * x + y * y + z * z;
        double r = sqrt(r2);
        // This allows us to compute the field at the nuclei correctly
        if (r < 1.0E-8) continue;

        Ex += mol->Z(i) * x / (r * r2);
        Ey += mol->Z(i) * y / (r * r2);
        Ez += mol->Z(i) * z / (r * r2);
    }

    Vector3 result(Ex, Ey, Ez);

    return result;
}


template <typename ContractionFunctor>
void ElectricFieldInt::compute_with_functor(ContractionFunctor functor, SharedMatrix coords) {
    set_chunks(3 * (coords->rowdim() + 2));
    std::vector<std::pair<double, std::array<double, 3>>> Zxyz;
    for (int site = 0; site < coords->rowdim(); ++site) {
        Zxyz.push_back({-1.0, {coords->get(0, site, 0), coords->get(0, site, 1), coords->get(0, site, 2)}});
    }
    bool bs1_is_bs2 = bs1_ == bs2_;
    engine0_->set_params(Zxyz);
    size_t npairs = shellpairs_.size();
    for (int shellpair = 0; shellpair < npairs; ++shellpair) {
        int idx1 = shellpairs_[shellpair].first;
        int idx2 = shellpairs_[shellpair].second;

        const auto &s1 = bs1_->l2_shell(idx1);
        const auto &s2 = bs2_->l2_shell(idx2);
        engine0_->compute(s1, s2);
        const auto &shell1 = bs1_->shell(idx1);
        const auto &shell2 = bs2_->shell(idx2);
        auto nbf1 = shell1.nfunction();
        auto nbf2 = shell2.nfunction();
        auto offset_1 = shell1.function_index();
        auto offset_2 = shell2.function_index();
        for (int site = 0; site < coords->rowdim(); ++site) {
            // The +2 here is because we're using an L2 derivative engine, which gives us the
            // bra atom derivatives, then the ket atom derivatives, then the derivatives w.r.t.
            // atomic charge positions, which are what we want.
            const double *bufferx = engine0_->results()[3 * (2+site) + 0];
            const double *buffery = engine0_->results()[3 * (2+site) + 1];
            const double *bufferz = engine0_->results()[3 * (2+site) + 2];
            for (int bf1 = 0; bf1 < nbf1; ++bf1) {
                for (int bf2 = 0; bf2 < nbf2; ++bf2) {
                    functor(bf1 + offset_1, bf2 + offset_2, site, *bufferx, *buffery, *bufferz);
                    ++bufferx;
                    ++buffery;
                    ++bufferz;
                }
            }
            if (s1 != s2 && bs1_is_bs2) {
                bufferx = engine0_->results()[3 * (2+site) + 0];
                buffery = engine0_->results()[3 * (2+site) + 1];
                bufferz = engine0_->results()[3 * (2+site) + 2];
                for (int bf1 = 0; bf1 < nbf1; ++bf1) {
                    for (int bf2 = 0; bf2 < nbf2; ++bf2) {
                        functor(bf2 + offset_2, bf1 + offset_1, site, *bufferx, *buffery, *bufferz);
                        ++bufferx;
                        ++buffery;
                        ++bufferz;
                    }
                }
            }
        }
    }
}

template void ElectricFieldInt::compute_with_functor(ContractOverDipolesFunctor, SharedMatrix);
template void ElectricFieldInt::compute_with_functor(ContractOverDensityFieldFunctor, SharedMatrix);

