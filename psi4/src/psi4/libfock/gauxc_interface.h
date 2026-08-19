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

#ifndef GAUXC_INTERFACE_H
#define GAUXC_INTERFACE_H

#ifdef USING_gauxc

#include <memory>
#include <string>

#include <Eigen/Core>

#include <gauxc/basisset.hpp>
#include <gauxc/molecule.hpp>
#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/runtime_environment.hpp>

namespace psi {

class BasisSet;
class Molecule;

/// Conversions between Psi4 and GauXC conventions, shared by the snLinK
/// exchange path and the GauXC XC quadrature.
namespace gauxc_interface {

/// Maps the Psi4 DFT_PRUNING_SCHEME / SNLINK_PRUNING_SCHEME string onto the
/// GauXC enum. Throws PSIEXCEPTION if GauXC has no equivalent.
GauXC::PruningScheme to_gauxc_pruning_scheme(const std::string& psi4_name);

/// Maps the Psi4 DFT_RADIAL_SCHEME / SNLINK_RADIAL_SCHEME string onto the
/// GauXC enum. Throws PSIEXCEPTION if GauXC has no equivalent.
GauXC::RadialQuad to_gauxc_radial_scheme(const std::string& psi4_name);

/// Maps a Psi4 DFT_NUCLEAR_SCHEME value to the GauXC atomic partition weight scheme.
/// Throws for values GauXC does not implement.
GauXC::XCWeightAlg to_gauxc_weight_scheme(const std::string& psi4_name);

/// Permutation between Psi4's Gaussian-ordered and GauXC's CCA-ordered
/// spherical harmonics.
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> generate_permutation_matrix(
    const std::shared_ptr<BasisSet> psi4_basisset, int gauxc_max_am);

/// Psi4::Molecule -> GauXC::Molecule
GauXC::Molecule psi4_to_gauxc_molecule(std::shared_ptr<Molecule> psi4_molecule);

/// Psi4::BasisSet -> GauXC::BasisSet
template <typename T>
GauXC::BasisSet<T> psi4_to_gauxc_basisset(std::shared_ptr<BasisSet> psi4_basisset, double basis_tol,
                                          bool force_cartesian);

/// Builds a host or device runtime environment. Throws if a device is
/// requested but GauXC was built without device support.
std::unique_ptr<GauXC::RuntimeEnvironment> make_runtime_environment(bool use_gpu, double gpu_mem_fraction);

}  // namespace gauxc_interface
}  // namespace psi

#endif  // USING_gauxc
#endif  // GAUXC_INTERFACE_H
