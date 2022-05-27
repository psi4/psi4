#include "jk.h"
#include "composite.h"

#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"
#include "psi4/libqt/qt.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <vector>
#include <algorithm>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

using namespace psi;

namespace psi {

SplitJKBase::SplitJKBase(std::shared_ptr<BasisSet> primary, Options& options) : primary_(primary), options_(options) {

    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
    bench_ = options_.get_int("BENCH");

    nthread_ = 1;
#ifdef _OPENMP
    nthread_ = Process::environment.get_n_threads();
#endif

    early_screening_ = false;
    lr_symmetric_ = true;
}

CompositeJK::CompositeJK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, 
                         std::string jtype, std::string ktype, Options& options)
                        : JK(primary), auxiliary_(auxiliary), jtype_(jtype), ktype_(ktype), options_(options) { common_init(); }

void CompositeJK::common_init() {

    nthread_ = 1;
#ifdef _OPENMP
    nthread_ = Process::environment.get_n_threads();
#endif

    if (jtype_ == "DIRECT_DF") {
        jalgo_ = std::make_shared<DirectDFJ>(primary_, auxiliary_, options_);
    } else {
        throw PSIEXCEPTION("J BUILD TYPE " + jtype_ + " IS NOT SUPPORTED IN COMPOSITE JK!");
    }

    if (ktype_ == "LINK") {
        kalgo_ = std::make_shared<LinK>(primary_, options_);
    } else if (ktype_ == "COSK") {
        kalgo_ = std::make_shared<COSK>(primary_, options_);
        early_screening_ = true;
    } else {
        throw PSIEXCEPTION("K BUILD TYPE " + ktype_ + " IS NOT SUPPORTED IN COMPOSITE JK!");
    }
}

size_t CompositeJK::memory_estimate() {
    return 0;   // only Direct-based integral algorithms are currently implemented in Composite JK
}

void CompositeJK::preiterations() { return; }
void CompositeJK::postiterations() { return; }

void CompositeJK::print_header() const {
    std::string screen_type = options_.get_str("SCREENING");
    if (print_) {
        outfile->Printf("  ==> CompositeJK: Mix/Match JK Builds <==\n\n");

        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        if (do_J_) outfile->Printf("    J Algorithm:       %11s\n", jtype_.c_str());
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        if (do_K_) outfile->Printf("    K Algorithm:       %11s\n", ktype_.c_str());
        outfile->Printf("    Integrals threads: %11d\n", nthread_);
        outfile->Printf("    Incremental Fock:  %11s\n", incfock_ ? "Yes" : "No");
        outfile->Printf("\n");
    }
    jalgo_->print_header();
    kalgo_->print_header();
}

void CompositeJK::compute_JK() {

    if (do_wK_) throw PSIEXCEPTION("Composite K algorithms do not currently support wK integrals!");

    if (incfock_) {
        timer_on("CompositeJK: INCFOCK Preprocessing");
        incfock_setup();
        int reset = options_.get_int("INCFOCK_FULL_FOCK_EVERY");
        double dconv = options_.get_double("INCFOCK_CONVERGENCE");
        double Dnorm = Process::environment.globals["SCF D NORM"];
        // Do IFB on this iteration?
        do_incfock_iter_ = (Dnorm >= dconv) && !initial_iteration_ && (incfock_count_ % reset != reset - 1);

        if (!initial_iteration_ && (Dnorm >= dconv)) incfock_count_ += 1;
        timer_off("CompositeJK: INCFOCK Preprocessing");
    }

    // Matrices to use/build depending on whether or not incremental Fock build is performed in the iteration
    std::vector<SharedMatrix>& D_ref = (do_incfock_iter_ ? delta_D_ao_ : D_ao_);
    std::vector<SharedMatrix>& J_ref = (do_incfock_iter_ ? delta_J_ao_ : J_ao_);
    std::vector<SharedMatrix>& K_ref = (do_incfock_iter_ ? delta_K_ao_ : K_ao_);

    if (do_J_) {
        jalgo_->set_lr_symmetric(lr_symmetric_);
        jalgo_->set_early_screening(early_screening_);
        jalgo_->build_G_component(D_ref, J_ref);
    }
    if (do_K_) {
        kalgo_->set_lr_symmetric(lr_symmetric_);
        kalgo_->set_early_screening(early_screening_);
        kalgo_->build_G_component(D_ref, K_ref);
    }

    if (incfock_) {
        timer_on("CompositeJK: INCFOCK Postprocessing");
        incfock_postiter();
        timer_off("CompositeJK: INCFOCK Postprocessing");
    }

    if (initial_iteration_) initial_iteration_ = false;
}

} // namespace psi