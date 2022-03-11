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
#include "psi4/libmints/integral.h"
#include "psi4/libmints/shellrotation.h"
#include "psi4/libmints/cartesianiter.h"
#include "psi4/libmints/rel_potential.h"
#include "psi4/libmints/electricfield.h"
#include "psi4/libmints/tracelessquadrupole.h"
#include "psi4/libmints/multipolepotential.h"
#include "psi4/libmints/eri.h"
#include "psi4/libmints/multipoles.h"
#include "psi4/libmints/quadrupole.h"
#include "psi4/libmints/angularmomentum.h"
#include "psi4/libmints/nabla.h"
#include "psi4/libmints/dipole.h"
#include "psi4/libmints/electrostatic.h"
#include "psi4/libmints/kinetic.h"
#include "psi4/libmints/3coverlap.h"
#include "psi4/libmints/overlap.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/potentialint.h"
#ifdef USING_ecpint
#include "psi4/libmints/ecpint.h"
#endif
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/erd_eri.h"

#ifdef USING_simint
#include "psi4/libmints/siminteri.h"
#endif

using namespace psi;

IntegralFactory::IntegralFactory(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                 std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4) {
    set_basis(bs1, bs2, bs3, bs4);
}

IntegralFactory::IntegralFactory(std::shared_ptr<BasisSet> bs1) { set_basis(bs1, bs1, bs1, bs1); }

IntegralFactory::~IntegralFactory() {}

std::shared_ptr<BasisSet> IntegralFactory::basis1() const { return bs1_; }

std::shared_ptr<BasisSet> IntegralFactory::basis2() const { return bs2_; }

std::shared_ptr<BasisSet> IntegralFactory::basis3() const { return bs3_; }

std::shared_ptr<BasisSet> IntegralFactory::basis4() const { return bs4_; }

void IntegralFactory::set_basis(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                std::shared_ptr<BasisSet> bs3, std::shared_ptr<BasisSet> bs4) {
    bs1_ = bs1;
    bs2_ = bs2;
    bs3_ = bs3;
    bs4_ = bs4;

    // Use the max am from libint
    init_spherical_harmonics(8);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_overlap(int deriv) {
    return std::make_unique<OverlapInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_overlap(int deriv) {
    std::shared_ptr<OneBodyAOInt> ao_int(ao_overlap(deriv));
    return std::make_unique<OneBodySOInt>(ao_int, this);
}

std::unique_ptr<ThreeCenterOverlapInt> IntegralFactory::overlap_3c() { return std::make_unique<ThreeCenterOverlapInt>(bs1_, bs2_, bs3_); }

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_kinetic(int deriv) {
    return std::make_unique<KineticInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_kinetic(int deriv) {
    std::shared_ptr<OneBodyAOInt> ao_int(ao_kinetic(deriv));
    return std::make_unique<OneBodySOInt>(ao_int, this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_potential(int deriv) {
    return std::make_unique<PotentialInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_potential(int deriv) {
    std::shared_ptr<OneBodyAOInt> ao_int(ao_potential(deriv));
    return  std::make_unique<PotentialSOInt>(ao_int, this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_ecp(int deriv) {
#ifdef USING_ecpint
    return std::make_unique<ECPInt>(spherical_transforms_, bs1_, bs2_, deriv);
#else
    throw PSIEXCEPTION("ECP shells requested but libecpint addon not enabled. Re-compile with `-D ENABLE_ecpint=ON`.");
#endif
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_ecp(int deriv) {
#ifdef USING_ecpint
    std::shared_ptr<OneBodyAOInt> ao_int(ao_ecp(deriv));
    return  std::make_unique<ECPSOInt>(ao_int, this);
#else
    throw PSIEXCEPTION("ECP shells requested but libecpint addon not enabled. Re-compile with `-D ENABLE_ecpint=ON`.");
#endif
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_rel_potential(int deriv) {
    return  std::make_unique<RelPotentialInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_rel_potential(int deriv) {
    std::shared_ptr<OneBodyAOInt> ao_int(ao_rel_potential(deriv));
    return std::make_unique<RelPotentialSOInt>(ao_int, this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::electrostatic() { return std::make_unique<ElectrostaticInt>(spherical_transforms_, bs1_, bs2_, 0); }

std::unique_ptr<OneBodyAOInt> IntegralFactory::pcm_potentialint() { return  std::make_unique<PCMPotentialInt>(spherical_transforms_, bs1_, bs2_, 0); }

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_pseudospectral(double omega, int deriv) {
    return new PseudospectralInt(spherical_transforms_, bs1_, bs2_, omega, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_pseudospectral(double omega, int deriv) {
    std::shared_ptr<OneBodyAOInt> ao_int(ao_pseudospectral(omega, deriv));
    return new OneBodySOInt(ao_int, this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_dipole(int deriv) { return  std::make_unique<DipoleInt>(spherical_transforms_, bs1_, bs2_, deriv); }

std::unique_ptr<OneBodySOInt> IntegralFactory::so_dipole(int deriv) {
    std::shared_ptr<OneBodyAOInt> ao_int(ao_dipole(deriv));
    return  std::make_unique<OneBodySOInt>(ao_int, this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_nabla(int deriv) { return  std::make_unique<NablaInt>(spherical_transforms_, bs1_, bs2_, deriv); }

std::unique_ptr<OneBodySOInt> IntegralFactory::so_nabla(int deriv) {
    std::shared_ptr<OneBodyAOInt> ao_int(ao_nabla(deriv));
    return std::make_unique<OneBodySOInt>(ao_int, this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_angular_momentum(int deriv) {
    return std::make_unique<AngularMomentumInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_angular_momentum(int deriv) {
    std::shared_ptr<OneBodyAOInt> ao_int(ao_angular_momentum(deriv));
    return std::make_unique<OneBodySOInt>(ao_int, this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_quadrupole() { return  std::make_unique<QuadrupoleInt>(spherical_transforms_, bs1_, bs2_); }

std::unique_ptr<OneBodySOInt> IntegralFactory::so_quadrupole() {
    std::shared_ptr<OneBodyAOInt> ao_int(ao_quadrupole());
    return std::make_unique<OneBodySOInt>(ao_int, this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_multipoles(int order, int deriv) {
    return  std::make_unique<MultipoleInt>(spherical_transforms_, bs1_, bs2_, order, deriv);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_multipole_potential(int order, int deriv) {
    return  std::make_unique<MultipolePotentialInt>(spherical_transforms_, bs1_, bs2_, order, deriv);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_multipoles(int order, int deriv) {
    std::shared_ptr<OneBodyAOInt> ao_int(ao_multipoles(order, deriv));
    return  std::make_unique<OneBodySOInt>(ao_int, this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_traceless_quadrupole() {
    return  std::make_unique<TracelessQuadrupoleInt>(spherical_transforms_, bs1_, bs2_);
}

std::unique_ptr<OneBodySOInt> IntegralFactory::so_traceless_quadrupole() {
    std::shared_ptr<OneBodyAOInt> ao_int(ao_traceless_quadrupole());
    return std::make_unique<OneBodySOInt>(ao_int, this);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::electric_field(int deriv) {
    return  std::make_unique<ElectricFieldInt>(spherical_transforms_, bs1_, bs2_, deriv);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erd_eri(int deriv, bool use_shell_pairs, bool needs_exchange) {
    auto integral_package = Process::environment.options.get_str("INTEGRAL_PACKAGE");
    auto threshold = Process::environment.options.get_double("INTS_TOLERANCE");
#ifdef USING_simint
    if (deriv == 0 && integral_package == "SIMINT") return std::make_unique<SimintERI>(this, deriv, use_shell_pairs, needs_exchange);
#endif
    if (integral_package == "LIBINT2") return  std::make_unique<Libint2ERI>(this, threshold, deriv, use_shell_pairs, needs_exchange);
#ifdef USING_erd
    if (deriv == 0 && integral_package == "ERD") return  std::make_unique<ERDERI>(this, deriv, use_shell_pairs);
#endif
    if (deriv > 0 && integral_package != "LIBINT1")
        outfile->Printf("ERI derivative integrals only available using Libint");
    if (integral_package == "SIMINT" || integral_package == "ERD")
        outfile->Printf("Chosen integral package " + integral_package +
                        " unavailable.\nRecompile with the appropriate option set.\nFalling back to Libint");
#ifdef ENABLE_Libint1t
    return std::make_unique<ERI>(this, deriv, use_shell_pairs);
#endif
    throw PSIEXCEPTION("No ERI object to return.");
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::eri(int deriv, bool use_shell_pairs, bool needs_exchange) {
    auto integral_package = Process::environment.options.get_str("INTEGRAL_PACKAGE");
    auto threshold = Process::environment.options.get_double("INTS_TOLERANCE");
#ifdef USING_simint
    if (deriv == 0 && integral_package == "SIMINT") return  std::make_unique<SimintERI>(this, deriv, use_shell_pairs, needs_exchange);
#endif
    if (integral_package == "LIBINT2") return  std::make_unique<Libint2ERI>(this, threshold, deriv, use_shell_pairs, needs_exchange);
#ifdef USING_erd
    if (deriv == 0 && integral_package == "ERD") return std::make_unique<ERDERI>(this, deriv, use_shell_pairs);
#endif
    if (deriv > 0 && integral_package != "LIBINT1")
        outfile->Printf("ERI derivative integrals only available using Libint");
    if (integral_package == "SIMINT" || integral_package == "ERD")
        outfile->Printf("Chosen integral package " + integral_package +
                        " unavailable.\nRecompile with the appropriate option set.\nFalling back to Libint");
#ifdef ENABLE_Libint1t
    return std::make_unique<ERI>(this, deriv, use_shell_pairs);
#endif
    throw PSIEXCEPTION("No ERI object to return.");
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erf_eri(double omega, int deriv, bool use_shell_pairs, bool needs_exchange) {
    auto integral_package = Process::environment.options.get_str("INTEGRAL_PACKAGE");
    auto threshold = Process::environment.options.get_double("INTS_TOLERANCE");
    if (integral_package == "LIBINT2")
        return std::make_unique<Libint2ErfERI>(omega, this, threshold, deriv, use_shell_pairs, needs_exchange);
#ifdef ENABLE_Libint1t
    return std::make_unique<ErfERI>(omega, this, deriv, use_shell_pairs);
#endif
    throw PSIEXCEPTION("No ERI object to return.");
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erf_complement_eri(double omega, int deriv, bool use_shell_pairs, bool needs_exchange) {
    auto integral_package = Process::environment.options.get_str("INTEGRAL_PACKAGE");
    auto threshold = Process::environment.options.get_double("INTS_TOLERANCE");
    if (integral_package == "LIBINT2")
        return std::make_unique<Libint2ErfComplementERI>(omega, this, threshold, deriv, use_shell_pairs, needs_exchange);
#ifdef ENABLE_Libint1t
    return std::make_unique<ErfComplementERI>(omega, this, deriv, use_shell_pairs);
#endif
    throw PSIEXCEPTION("No ERI object to return.");
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::yukawa_eri(double zeta, int deriv, bool use_shell_pairs, bool needs_exchange) {
    auto threshold = Process::environment.options.get_double("INTS_TOLERANCE");
    return std::make_unique<Libint2YukawaERI>(zeta, this, threshold, deriv, use_shell_pairs, needs_exchange);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12(std::vector<std::pair<double, double>> exp_coeff, int deriv, bool use_shell_pairs) {
    return  std::make_unique<Libint2F12>(exp_coeff, this, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12_squared(std::vector<std::pair<double, double>> exp_coeff, int deriv,
                                           bool use_shell_pairs) {
    return  std::make_unique<Libint2F12Squared>(exp_coeff, this, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12g12(std::vector<std::pair<double, double>> exp_coeff, int deriv,
                                      bool use_shell_pairs) {
    return std::make_unique<Libint2F12G12>(exp_coeff, this, deriv, use_shell_pairs);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::f12_double_commutator(std::vector<std::pair<double, double>> exp_coeff, int deriv,
                                                     bool use_shell_pairs) {
    return std::make_unique<Libint2F12DoubleCommutator>(exp_coeff, this, deriv, use_shell_pairs);
}

void IntegralFactory::init_spherical_harmonics(int max_am) {
    spherical_transforms_.clear();
    ispherical_transforms_.clear();

    for (int i = 0; i <= max_am; ++i) {
        spherical_transforms_.push_back(SphericalTransform(i));
        ispherical_transforms_.push_back(ISphericalTransform(i));
    }
}

AOShellCombinationsIterator IntegralFactory::shells_iterator() {
    return AOShellCombinationsIterator(bs1_, bs2_, bs3_, bs4_);
}

AOShellCombinationsIterator* IntegralFactory::shells_iterator_ptr() {
    return new AOShellCombinationsIterator(bs1_, bs2_, bs3_, bs4_);
}

AOIntegralsIterator IntegralFactory::integrals_iterator(int p, int q, int r, int s) {
    return AOIntegralsIterator(bs1_->shell(p), bs2_->shell(q), bs3_->shell(r), bs4_->shell(s));
}

CartesianIter* IntegralFactory::cartesian_iter(int l) const { return new CartesianIter(l); }

RedundantCartesianIter* IntegralFactory::redundant_cartesian_iter(int l) const { return new RedundantCartesianIter(l); }

RedundantCartesianSubIter* IntegralFactory::redundant_cartesian_sub_iter(int l) const {
    return new RedundantCartesianSubIter(l);
}

ShellRotation IntegralFactory::shell_rotation(int am, SymmetryOperation& so, int pure) const {
    ShellRotation r(am, so, this, pure);
    return r;
}

SphericalTransformIter* IntegralFactory::spherical_transform_iter(int am, int inv, int subl) const {
    if (subl != -1) throw NOT_IMPLEMENTED_EXCEPTION();

    if (inv) {
        return new SphericalTransformIter(ispherical_transforms_[am]);
    }
    return new SphericalTransformIter(spherical_transforms_[am]);
}
