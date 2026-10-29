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

#include <algorithm>
#include <array>
#include <functional>
#include <mutex>
#include <set>
#include <string>

#include <libint2/util/configuration.h>

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

#ifdef USING_cint
#include "psi4/libmints/libcinteri.h"
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

namespace {

/// Stop unless INTEGRAL_PACKAGE names an engine compiled into this build. AUTO and LIBINT2 always are.
void require_compiled_in(const std::string &package) {
    bool compiled = package == "AUTO" || package == "LIBINT2";
#ifdef USING_simint
    compiled |= package == "SIMINT";
#endif
#ifdef USING_cint
    compiled |= package == "LIBCINT";
#endif
    if (compiled) return;
    const std::string option = package == "LIBCINT" ? "cint" : "simint";  // the CMake ENABLE_<option>
    throw PSIEXCEPTION("INTEGRAL_PACKAGE " + package + " was requested, but Psi4 was built without it. " +
                       "Rebuild with `-D ENABLE_" + option + "=ON`.");
}

/// The engines IntegralFactory::eri() has built since the last reset.
std::mutex coulomb_eri_lock;
std::set<std::string> coulomb_eri_built;

void record_coulomb_eri_engine(const std::string &engine) {
    std::lock_guard<std::mutex> guard(coulomb_eri_lock);
    coulomb_eri_built.insert(engine);
}

/// Print a note once per process, not once per (per-thread) integral object.
void note_once(const std::string &note) {
    static std::mutex lock;
    static std::set<std::string> announced;
    std::lock_guard<std::mutex> guard(lock);
    if (announced.insert(note).second) outfile->Printf("  %s\n", note.c_str());
}

/// INTEGRAL_PACKAGE's engine provides no integrals of this kind, so Libint2 supplies them.
void note_libint2_fallback(const std::string &package, const std::string &missing) {
    if (package == "AUTO" || package == "LIBINT2") return;
    note_once("INTEGRAL_PACKAGE " + package + ": provides no " + missing + "; using Libint2 instead.");
}

using TwoBodyMaker = std::function<std::unique_ptr<TwoBodyAOInt>()>;

/// "Coulomb integrals", "Coulomb first-derivative integrals", ...
std::string integral_kind(const std::string &op, int deriv) {
    return op + (deriv == 1 ? " first-derivative" : deriv == 2 ? " second-derivative" : "") + " integrals";
}

/// The highest angular momentum of Libint2's Coulomb-type integrals over 2, 3 or 4 centers, by derivative order, or -1
/// where Libint2 has none. Libint2 (>= 2.13) names each component of its configuration by the angular momentum common
/// to all centers, e.g., eri_kkkk_d0 for 4-center integrals up to k.
int libint2_lmax(int centers, int deriv) {
    static const auto table = []() {
        const std::string am = "spdfghiklmnoqrtuvwxyz";
        std::array<std::array<int, 3>, 3> lmax;
        for (int c = 2; c <= 4; ++c)
            for (int d = 0; d <= 2; ++d) {
                lmax[c - 2][d] = -1;
                for (int l = 0; l < static_cast<int>(am.size()); ++l) {
                    const std::string component = "eri_" + std::string(c, am[l]) + "_d" + std::to_string(d);
                    if (libint2::supports(component)) lmax[c - 2][d] = l;
                }
            }
        return lmax;
    }();
    return table[centers - 2][deriv];
}

/// Libint2's integrals over these bases, or, under AUTO, libcint's when Libint2 lacks them: beyond the angular
/// momentum it was built for, or at a derivative order it was built without (Hessians on some platforms).
/// make_libcint is empty without libcint. Sets engine to the one used.
std::unique_ptr<TwoBodyAOInt> libint2_unless_beyond(const std::string &package, const std::string &op, int deriv,
                                                    const std::vector<const BasisSet *> &bases,
                                                    const TwoBodyMaker &make_libint2, const TwoBodyMaker &make_libcint,
                                                    std::string &engine) {
    // The density-fitting dummy basis is not a center.
    int centers = 0, l = 0;
    for (const auto *bs : bases) {
        if (bs->nshell() == 1 && bs->shell(0).am() == 0 && bs->shell(0).exp(0) == 0.0) continue;
        ++centers;
        l = std::max(l, bs->max_am());
    }
    const int lmax = centers >= 2 ? libint2_lmax(centers, deriv) : 0;
    // A Libint2 that reports no 4-center energy integrals reports nothing at all: let it try.
    if (l <= lmax || libint2_lmax(4, 0) < 0) {
        engine = "libint2";
        return make_libint2();
    }
    const std::string kind = std::to_string(centers) + "-center " + integral_kind(op, deriv);
    const std::string lacks = lmax < 0 ? "Libint2 was built without " + kind
                                       : "Libint2 was built for " + kind + " up to angular momentum " +
                                             std::to_string(lmax) + ", but these need " + std::to_string(l);
    if (package == "AUTO" && make_libcint) {
        note_once("INTEGRAL_PACKAGE AUTO: " + lacks + "; using libcint instead.");
        engine = "libcint";
        return make_libcint();
    }
    const std::string remedy = make_libcint ? "INTEGRAL_PACKAGE AUTO or LIBCINT computes them with libcint."
                                            : "Rebuild with `-D ENABLE_cint=ON` for libcint to compute them.";
    throw PSIEXCEPTION("INTEGRAL_PACKAGE " + package + ": " + lacks + ". " + remedy);
}

}  // namespace

std::unique_ptr<TwoBodyAOInt> IntegralFactory::eri(int deriv, bool use_shell_pairs, bool needs_exchange) {
    const auto package = Process::environment.options.get_str("INTEGRAL_PACKAGE");
    const auto threshold = Process::environment.options.get_double("INTS_TOLERANCE");
    require_compiled_in(package);
#ifdef USING_cint
    if (package == "LIBCINT") {
        record_coulomb_eri_engine("libcint");
        return std::make_unique<LibcintERI>(this, deriv, use_shell_pairs, needs_exchange);
    }
#endif
#ifdef USING_simint
    if (package == "SIMINT" && deriv == 0) {
        record_coulomb_eri_engine("simint");
        return std::make_unique<SimintERI>(this, deriv, use_shell_pairs, needs_exchange);
    }
#endif
    note_libint2_fallback(package, "derivative integrals");
    TwoBodyMaker make_libcint;
#ifdef USING_cint
    make_libcint = [&] { return std::make_unique<LibcintERI>(this, deriv, use_shell_pairs, needs_exchange); };
#endif
    std::string engine;
    auto ints = libint2_unless_beyond(
        package, "Coulomb", deriv, {bs1_.get(), bs2_.get(), bs3_.get(), bs4_.get()},
        [&] { return std::make_unique<Libint2ERI>(this, threshold, deriv, use_shell_pairs, needs_exchange); },
        make_libcint, engine);
    record_coulomb_eri_engine(engine);
    return ints;
}

void IntegralFactory::reset_coulomb_eri_engines() {
    std::lock_guard<std::mutex> guard(coulomb_eri_lock);
    coulomb_eri_built.clear();
}

std::string IntegralFactory::coulomb_eri_engines() {
    std::lock_guard<std::mutex> guard(coulomb_eri_lock);
    std::string engines;
    for (const char *engine : {"simint", "libint2", "libcint"})
        if (coulomb_eri_built.count(engine)) engines += (engines.empty() ? "" : ">") + std::string(engine);
    return engines;
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erf_eri(double omega, int deriv, bool use_shell_pairs, bool needs_exchange) {
    const auto package = Process::environment.options.get_str("INTEGRAL_PACKAGE");
    const auto threshold = Process::environment.options.get_double("INTS_TOLERANCE");
    require_compiled_in(package);
#ifdef USING_cint
    if (package == "LIBCINT")
        return std::make_unique<LibcintErfERI>(omega, this, deriv, use_shell_pairs, needs_exchange);
#endif
    note_libint2_fallback(package, "range-separated (erf) integrals");
    TwoBodyMaker make_libcint;
#ifdef USING_cint
    make_libcint = [&] { return std::make_unique<LibcintErfERI>(omega, this, deriv, use_shell_pairs, needs_exchange); };
#endif
    std::string engine;
    return libint2_unless_beyond(
        package, "range-separated (erf)", deriv, {bs1_.get(), bs2_.get(), bs3_.get(), bs4_.get()},
        [&] {
            return std::make_unique<Libint2ErfERI>(omega, this, threshold, deriv, use_shell_pairs, needs_exchange);
        },
        make_libcint, engine);
}

std::unique_ptr<TwoBodyAOInt> IntegralFactory::erf_complement_eri(double omega, int deriv, bool use_shell_pairs, bool needs_exchange) {
    const auto package = Process::environment.options.get_str("INTEGRAL_PACKAGE");
    const auto threshold = Process::environment.options.get_double("INTS_TOLERANCE");
    require_compiled_in(package);
#ifdef USING_cint
    if (package == "LIBCINT")
        return std::make_unique<LibcintErfComplementERI>(omega, this, deriv, use_shell_pairs, needs_exchange);
#endif
    note_libint2_fallback(package, "range-separated (erfc) integrals");
    TwoBodyMaker make_libcint;
#ifdef USING_cint
    make_libcint = [&] {
        return std::make_unique<LibcintErfComplementERI>(omega, this, deriv, use_shell_pairs, needs_exchange);
    };
#endif
    std::string engine;
    return libint2_unless_beyond(
        package, "range-separated (erfc)", deriv, {bs1_.get(), bs2_.get(), bs3_.get(), bs4_.get()},
        [&] {
            return std::make_unique<Libint2ErfComplementERI>(omega, this, threshold, deriv, use_shell_pairs,
                                                             needs_exchange);
        },
        make_libcint, engine);
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
