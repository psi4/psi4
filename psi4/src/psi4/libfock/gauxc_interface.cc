/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "gauxc_interface.h"

#ifdef USING_gauxc

#include <unordered_map>
#include <vector>

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/psi4-dec.h"

namespace psi {
namespace gauxc_interface {

GauXC::PruningScheme to_gauxc_pruning_scheme(const std::string& psi4_name) {
    static const std::unordered_map<std::string, GauXC::PruningScheme> map = {
        {"ROBUST", GauXC::PruningScheme::Robust},
        {"TREUTLER", GauXC::PruningScheme::Treutler},
        {"NONE", GauXC::PruningScheme::Unpruned},
    };
    auto it = map.find(psi4_name);
    if (it == map.end()) {
        throw PSIEXCEPTION("Pruning scheme '" + psi4_name +
                           "' has no GauXC equivalent. Supported by the GauXC path: ROBUST, TREUTLER, NONE.");
    }
    return it->second;
}

GauXC::RadialQuad to_gauxc_radial_scheme(const std::string& psi4_name) {
    static const std::unordered_map<std::string, GauXC::RadialQuad> map = {
        {"TREUTLER", GauXC::RadialQuad::TreutlerAhlrichs},
        {"MURA", GauXC::RadialQuad::MuraKnowles},
        // The Murray, Handy, Laming literature reference is mentioned in cubature.cc
        // with association to this keyword
        {"EM", GauXC::RadialQuad::MurrayHandyLaming},
    };
    auto it = map.find(psi4_name);
    if (it == map.end()) {
        throw PSIEXCEPTION("Radial scheme '" + psi4_name +
                           "' has no GauXC equivalent. Supported by the GauXC path: TREUTLER, MURA, EM.");
    }
    return it->second;
}

Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> generate_permutation_matrix(
    const std::shared_ptr<BasisSet> psi4_basisset, int gauxc_max_am) {
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> permutation_matrix(psi4_basisset->nbf());

    // general array for how to reorder integrals
    std::vector<int> cca_integral_order(2 * gauxc_max_am + 1, 0);

    // p shells or larger
    for (size_t l = 1, idx = 1; l != gauxc_max_am; idx += 2, ++l) {
        cca_integral_order[idx] = l;
        cca_integral_order[idx + 1] = -l;
    }

    // actually construct permutation matrix
    for (int ish = 0, ibf = 0; ish != psi4_basisset->nshell(); ++ish) {
        auto& sh = psi4_basisset->shell(ish);
        auto am = sh.am();

        auto ibf_base = ibf;
        for (int ishbf = 0; ishbf != 2 * am + 1; ++ishbf, ++ibf) {
            permutation_matrix.indices()[ibf] = ibf_base + cca_integral_order[ishbf] + am;
        }
    }

    return permutation_matrix;
}

GauXC::Molecule psi4_to_gauxc_molecule(std::shared_ptr<Molecule> psi4_molecule) {
    GauXC::Molecule gauxc_molecule;

    for (size_t iatom = 0; iatom != psi4_molecule->natom(); ++iatom) {
        auto atomic_number = psi4_molecule->true_atomic_number(iatom);
        auto x_coord = psi4_molecule->x(iatom);
        auto y_coord = psi4_molecule->y(iatom);
        auto z_coord = psi4_molecule->z(iatom);

        gauxc_molecule.emplace_back(GauXC::AtomicNumber(atomic_number), x_coord, y_coord, z_coord);
    }

    return gauxc_molecule;
}

template <typename T>
GauXC::BasisSet<T> psi4_to_gauxc_basisset(std::shared_ptr<BasisSet> psi4_basisset, double basis_tol,
                                          bool force_cartesian) {
    using prim_array = typename GauXC::Shell<T>::prim_array;
    using cart_array = typename GauXC::Shell<T>::cart_array;

    GauXC::BasisSet<T> gauxc_basisset(psi4_basisset->nshell());

    for (size_t ishell = 0; ishell != psi4_basisset->nshell(); ++ishell) {
        auto psi4_shell = psi4_basisset->shell(ishell);

        const auto nprim = GauXC::PrimSize(psi4_shell.nprimitive());
        prim_array alpha;
        prim_array coeff;

        for (size_t iprim = 0; iprim != psi4_shell.nprimitive(); ++iprim) {
            alpha.at(iprim) = psi4_shell.exp(iprim);
            coeff.at(iprim) = psi4_shell.coef(iprim);
        }

        auto psi4_shell_center = psi4_shell.center();
        cart_array center = {psi4_shell_center[0], psi4_shell_center[1], psi4_shell_center[2]};

        gauxc_basisset[ishell] = GauXC::Shell(
            nprim, GauXC::AngularMomentum(psi4_shell.am()),
            (force_cartesian ? GauXC::SphericalType(false) : GauXC::SphericalType(!(psi4_shell.is_cartesian()))), alpha,
            coeff, center,
            false  // do not normalize shell via GauXC; it is normalized via Psi4
        );
    }

    for (auto& sh : gauxc_basisset) {
        sh.set_shell_tolerance(basis_tol);
    }

    return gauxc_basisset;
}

template GauXC::BasisSet<double> psi4_to_gauxc_basisset<double>(std::shared_ptr<BasisSet>, double, bool);

std::unique_ptr<GauXC::RuntimeEnvironment> make_runtime_environment(bool use_gpu, double gpu_mem_fraction) {
#ifdef GAUXC_HAS_DEVICE
    if (use_gpu) {
        return std::make_unique<GauXC::DeviceRuntimeEnvironment>(GAUXC_MPI_CODE(MPI_COMM_WORLD, ) gpu_mem_fraction);
    }
#else
    if (use_gpu) {
        throw PSIEXCEPTION(
            "GPU execution was requested, but this GauXC installation was built without device support. "
            "Rebuild Psi4 against a device-enabled GauXC.");
    }
#endif
    return std::make_unique<GauXC::RuntimeEnvironment>(GAUXC_MPI_CODE(MPI_COMM_WORLD));
}

}  // namespace gauxc_interface
}  // namespace psi

#endif  // USING_gauxc
