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
#include "psi4/libmints/integral.h"

#include <functional>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <vector>

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
#include "psi4/libmints/potential_erf.h"
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

#ifdef USING_simint
#include "psi4/libmints/siminteri.h"
#endif

#ifdef USING_libcint
#include "psi4/libmints/libcinteri.h"
#endif

using namespace psi;

namespace {

std::mutex engine_notes_mutex;
std::set<std::string> engine_notes;

/// Name a class of two-electron integrals, e.g., "ERI" or "erf ERI 1st deriv".
std::string eri_kind(const std::string &base, int deriv) {
    if (deriv == 0) return base;
    return base + (deriv == 1 ? " 1st deriv" : (deriv == 2 ? " 2nd deriv" : " deriv" + std::to_string(deriv)));
}

/// Report to the output file which library computes a class of two-electron integrals.
/// Printed once per class and engine (until IntegralFactory::reset_engine_notes, i.e.,
/// psi4.core.clean()), from the factory, once the object is built.
void note_engine(const std::string &kind, const std::string &engine) {
    std::lock_guard<std::mutex> lock(engine_notes_mutex);
    if (engine_notes.insert(kind + "|" + engine).second)
        outfile->Printf("  Two-electron integrals (%s) from %s.\n", kind.c_str(), engine.c_str());
}

/// The libraries that compute two-electron integrals.
enum class EriEngine { Libint2, Libcint, Simint };

/// How each engine constructs one class of two-electron integral; absent or empty where
/// the engine can't compute that class.
using EriMakers = std::map<EriEngine, std::function<std::unique_ptr<TwoBodyAOInt>()>>;

std::string engine_name(EriEngine engine) {
    switch (engine) {
        case EriEngine::Libint2:
            return "Libint2";
        case EriEngine::Libcint:
            return "libcint";
        case EriEngine::Simint:
            return "Simint";
    }
    return "";
}

/// Why Psi4 can't use the engine at all, or empty if it was built with it.
std::string engine_missing(EriEngine engine) {
#ifndef USING_libcint
    if (engine == EriEngine::Libcint) return "Psi4 wasn't built with libcint (`-D ENABLE_libcint=ON`)";
#endif
#ifndef USING_simint
    if (engine == EriEngine::Simint) return "Psi4 wasn't built with Simint (`-D ENABLE_simint=ON`)";
#endif
    return "";
}

/// The engines INTEGRAL_PACKAGE tries, in order. The first that can compute a class of
/// two-electron integrals supplies it; the _ONLY packages never fall back.
std::vector<EriEngine> eri_engine_chain(const std::string &integral_package) {
    if (integral_package == "LIBINT2") return {EriEngine::Libint2};
    if (integral_package == "LIBCINT") return {EriEngine::Libcint, EriEngine::Libint2};
    if (integral_package == "LIBCINT_ONLY") return {EriEngine::Libcint};
    if (integral_package == "SIMINT") return {EriEngine::Simint};
    throw PSIEXCEPTION("INTEGRAL_PACKAGE " + integral_package + " has no two-electron integral engines.");
}

/// Construct `kind` integrals from the first engine in INTEGRAL_PACKAGE's chain that can
/// compute them, and report to the output file which engine that is and why any earlier
/// ones were passed over. Throws if the package's own library is missing or no engine can.
std::unique_ptr<TwoBodyAOInt> build_from_chain(const std::string &kind, const EriMakers &makers) {
    const auto integral_package = Process::environment.options.get_str("INTEGRAL_PACKAGE");
    const auto chain = eri_engine_chain(integral_package);
    std::string passed;  // why earlier engines couldn't
    for (const auto engine : chain) {
        std::string why = engine_missing(engine);
        if (!why.empty() && engine == chain.front())
            throw PSIEXCEPTION("INTEGRAL_PACKAGE " + integral_package + " unavailable: " + why + ".");
        const auto maker = makers.find(engine);
        if (why.empty() && (maker == makers.end() || !maker->second))
            why = engine_name(engine) + " doesn't compute them";
        if (why.empty()) {
            auto ints = maker->second();
            note_engine(kind, engine_name(engine) + (passed.empty() ? "" : " (fallback: " + passed + ")"));
            return ints;
        }
        passed += (passed.empty() ? "" : "; ") + why;
    }
    std::string hint;
    if (const auto only = integral_package.rfind("_ONLY"); only != std::string::npos)
        hint = " INTEGRAL_PACKAGE " + integral_package.substr(0, only) + " falls back to other libraries.";
    throw PSIEXCEPTION("INTEGRAL_PACKAGE " + integral_package + " can't supply " + kind + " integrals: " + passed +
                       "." + hint);
}

}  // namespace

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

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_potential_erf(double omega, int deriv) {
    return std::make_unique<PotentialErfInt>(spherical_transforms_, bs1_, bs2_, omega, deriv);
}

std::unique_ptr<OneBodyAOInt> IntegralFactory::ao_potential_erf_complement(double omega, int deriv) {
    return std::make_unique<PotentialErfComplementInt>(spherical_transforms_, bs1_, bs2_, omega, deriv);
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

void IntegralFactory::reset_engine_notes() {
    std::lock_guard<std::mutex> lock(engine_notes_mutex);
    engine_notes.clear();
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::eri(int deriv, bool use_shell_pairs, bool needs_exchange) {
    auto threshold = Process::environment.options.get_double("INTS_TOLERANCE");
    EriMakers makers;
    makers[EriEngine::Libint2] = [&] {
        return std::make_unique<Libint2ERI>(this, threshold, deriv, use_shell_pairs, needs_exchange);
    };
#ifdef USING_libcint
    if (deriv == 0)
        makers[EriEngine::Libcint] = [&] {
            return std::make_unique<LibcintERI>(this, deriv, use_shell_pairs, needs_exchange);
        };
#endif
#ifdef USING_simint
    if (deriv == 0)
        makers[EriEngine::Simint] = [&] {
            return std::make_unique<SimintERI>(this, deriv, use_shell_pairs, needs_exchange);
        };
#endif
    return build_from_chain(eri_kind("ERI", deriv), makers);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erf_eri(double omega, int deriv, bool use_shell_pairs, bool needs_exchange) {
    auto threshold = Process::environment.options.get_double("INTS_TOLERANCE");
    EriMakers makers;
    makers[EriEngine::Libint2] = [&] {
        return std::make_unique<Libint2ErfERI>(omega, this, threshold, deriv, use_shell_pairs, needs_exchange);
    };
#ifdef USING_libcint
    if (deriv == 0)
        makers[EriEngine::Libcint] = [&] {
            return std::make_unique<LibcintErfERI>(omega, this, deriv, use_shell_pairs, needs_exchange);
        };
#endif
    return build_from_chain(eri_kind("erf ERI", deriv), makers);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erf_complement_eri(double omega, int deriv, bool use_shell_pairs, bool needs_exchange) {
    auto threshold = Process::environment.options.get_double("INTS_TOLERANCE");
    EriMakers makers;
    makers[EriEngine::Libint2] = [&] {
        return std::make_unique<Libint2ErfComplementERI>(omega, this, threshold, deriv, use_shell_pairs,
                                                         needs_exchange);
    };
#ifdef USING_libcint
    if (deriv == 0)
        makers[EriEngine::Libcint] = [&] {
            return std::make_unique<LibcintErfComplementERI>(omega, this, deriv, use_shell_pairs, needs_exchange);
        };
#endif
    return build_from_chain(eri_kind("erfc ERI", deriv), makers);
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
