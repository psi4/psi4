#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include "integralfunctors.h"

#include <libmints/mints.h>
#include "rohf.h"
#include <psi4-dec.h>

#define _DEBUG

using namespace std;
using namespace psi;

namespace psi { namespace scf {

ROHF::ROHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : HF(options, psio, chkpt)
{
    common_init();
}

ROHF::ROHF(Options& options, shared_ptr<PSIO> psio)
    : HF(options, psio)
{
    common_init();
}

ROHF::~ROHF() {
}

void ROHF::common_init()
{
    Fc_      = SharedMatrix(factory_->create_matrix("F closed"));
    Fo_      = SharedMatrix(factory_->create_matrix("F open"));
    Fa_      = SharedMatrix(factory_->create_matrix("F effective (MO basis)"));
    Fb_      = Fa_;
    Feff_    = Fa_;
    Ca_      = SharedMatrix(factory_->create_matrix("Moleular orbitals"));
    Cb_      = Ca_;
    Dc_      = SharedMatrix(factory_->create_matrix("D closed"));
    Do_      = SharedMatrix(factory_->create_matrix("D open"));
    Kc_      = SharedMatrix(factory_->create_matrix("K closed"));
    Ko_      = SharedMatrix(factory_->create_matrix("K open"));
    Dc_old_  = SharedMatrix(factory_->create_matrix("D closed old"));
    Do_old_  = SharedMatrix(factory_->create_matrix("D open old"));
    Gc_      = SharedMatrix(factory_->create_matrix("G closed"));
    Go_      = SharedMatrix(factory_->create_matrix("G open"));
    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = epsilon_a_;

    fprintf(outfile, "  DIIS %s.\n\n", diis_enabled_ ? "enabled" : "disabled");
}

void ROHF::finalize()
{
    Fc_.reset();
    Fo_.reset();
    Dc_old_.reset();
    Do_old_.reset();
    Gc_.reset();
    Go_.reset();

    HF::finalize();
}

void ROHF::form_initial_C()
{
    Matrix temp;
    Vector values;
    factory_->create_matrix(temp);
    factory_->create_vector(values);

    // In ROHF the creation of the C matrix depends on the previous iteration's C
    // matrix. Here we use H to generate the first C.
    temp.copy(H_);
    temp.transform(Shalf_);
    temp.diagonalize(Ca_, values);
    find_occupation();
    temp.gemm(false, false, 1.0, Shalf_, Ca_, 0.0);
    Ca_->copy(temp);
    Cb_->copy(temp);

    if (print_ > 3)
        Ca_->print(outfile, "initial C");
}

void ROHF::save_density_and_energy()
{
    Dc_old_->copy(Dc_); // save previous density
    Do_old_->copy(Do_); // save previous density
    Eold_ = E_; // save previous energy
}

void ROHF::save_information()
{
    // Print the final docc vector
    char **temp2 = molecule_->irrep_labels();

    // Stupid hack...
    if (!psio_->open_check(32))
        psio_->open(32, PSIO_OPEN_OLD);

    // TODO: Delete this as soon as possible!!!
    // Can't believe I'm adding this...
    chkpt_->wt_nirreps(factory_->nirrep());
    chkpt_->wt_irr_labs(temp2);

    int nso = basisset_->nbf();

    fprintf(outfile, "\n  Final DOCC vector = (");
    for (int h=0; h<factory_->nirrep(); ++h) {
        fprintf(outfile, "%2d %3s ", doccpi_[h], temp2[h]);
    }
    fprintf(outfile, ")\n");

    fprintf(outfile, "  Final SOCC vector = (");
    for (int h=0; h<factory_->nirrep(); ++h) {
        fprintf(outfile, "%2d %3s ", soccpi_[h], temp2[h]);
    }
    fprintf(outfile, ")\n");

    int print_mos = options_.get_bool("PRINT_MOS");
    if (print_mos) {
        fprintf(outfile, "\n  Molecular orbitals:\n");
        Ca_->eivprint(epsilon_a_);
    }

    // Print out orbital energies.
    std::vector<std::pair<double, int> > pairs;
    for (int h=0; h<epsilon_a_->nirrep(); ++h) {
        for (int i=0; i<epsilon_a_->dimpi()[h]; ++i)
            pairs.push_back(make_pair(epsilon_a_->get(h, i), h));
    }
    sort(pairs.begin(), pairs.end());
    int ndocc = 0, nsocc = 0;
    for (int i=0; i<epsilon_a_->nirrep(); ++i) {
        ndocc += doccpi_[i];
        nsocc += soccpi_[i];
    }

    fprintf(outfile,
            "\n  Orbital energies (a.u.):\n    Doubly occupied orbitals\n      ");
    for (int i=1; i<=ndocc; ++i) {
        fprintf(outfile, "%12.6f %3s  ", pairs[i-1].first,
                temp2[pairs[i-1].second]);
        if (i % 4 == 0)
            fprintf(outfile, "\n      ");
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "\n    Singly occupied orbitals\n      ");
    for (int i=ndocc+1; i<=ndocc+nsocc; ++i) {
        fprintf(outfile, "%12.6f %3s  ", pairs[i-1].first,
                temp2[pairs[i-1].second]);
        if ((i-ndocc) % 4 == 0)
            fprintf(outfile, "\n      ");
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "\n    Unoccupied orbitals\n      ");
    for (int i=ndocc+nsocc+1; i<=nso; ++i) {
        fprintf(outfile, "%12.6f %3s  ", pairs[i-1].first,
                temp2[pairs[i-1].second]);
        if ((i-ndocc-nsocc) % 4 == 0)
            fprintf(outfile, "\n      ");
    }
    fprintf(outfile, "\n");

    for (int i=0; i<epsilon_a_->nirrep(); ++i)
        free(temp2[i]);
    free(temp2);

    chkpt_->wt_nso(nso);
    chkpt_->wt_nmo(nso);
    chkpt_->wt_nao(nso);
    chkpt_->wt_ref(2); // ROHF
    chkpt_->wt_etot(E_);
    chkpt_->wt_escf(E_);
    chkpt_->wt_eref(E_);
    chkpt_->wt_clsdpi(doccpi_);
    chkpt_->wt_orbspi(epsilon_a_->dimpi());
    chkpt_->wt_openpi(soccpi_);
    chkpt_->wt_phase_check(0);

    Feff_->save(psio_, 32);

    // Figure out frozen core orbitals
    int nfzc = molecule_->nfrozen_core();
    int nfzv = options_.get_int("FREEZE_VIRT");
    int *frzcpi = compute_fcpi(nfzc, epsilon_a_);
    int *frzvpi = compute_fvpi(nfzv, epsilon_a_);
    chkpt_->wt_frzcpi(frzcpi);
    chkpt_->wt_frzvpi(frzvpi);
    delete[](frzcpi);
    delete[](frzvpi);

    int nopenirreps = 0;
    for (int i=0; i<epsilon_a_->nirrep(); ++i)
        if (soccpi_[i])
            nopenirreps++;

    // This code currently only handles ROHF
    chkpt_->wt_iopen(nopenirreps * (nopenirreps + 1));

    // Write eigenvectors and eigenvalue to checkpoint
    double *values = epsilon_a_->to_block_vector();
    chkpt_->wt_evals(values);
    free(values);
    double **vectors = Ca_->to_block_matrix();
    chkpt_->wt_scf(vectors);
    free_block(vectors);

}

void ROHF::save_fock()
{
    if (initialized_diis_manager_ == false) {
        diis_manager_ = shared_ptr<DIISManager>(new DIISManager(max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError, DIISManager::OnDisk, psio_));
        diis_manager_->set_error_vector_size(1, DIISEntry::Matrix, Feff_.get());
        diis_manager_->set_vector_size(1, DIISEntry::Matrix, Feff_.get());
        initialized_diis_manager_ = true;
    }

    // Save the effective Fock, back transform to AO, and orthonormalize
//    diis_F_[current_diis_fock_]->copy(Feff_);
//    diis_F_[current_diis_fock_]->back_transform(C_);
//    diis_F_[current_diis_fock_]->transform(Sphalf_);

//    // Determine error matrix for this Fock
//    diis_E_[current_diis_fock_]->copy(Feff_);
//    diis_E_[current_diis_fock_]->zero_diagonal();
//    diis_E_[current_diis_fock_]->back_transform(C_);
//    diis_E_[current_diis_fock_]->transform(Sphalf_);

//#ifdef _DEBUG
//    if (debug_) {
//        fprintf(outfile, "  New error matrix:\n");
//        diis_E_[current_diis_fock_]->print(outfile);
//    }
//#endif
//    current_diis_fock_++;
//    if (current_diis_fock_ == min_diis_vectors_)
//        current_diis_fock_ = 0;
}

bool ROHF::diis()
{
    return diis_manager_->extrapolate(1, Feff_.get());
}

bool ROHF::test_convergency()
{
    double ediff = E_ - Eold_;

    if (fabs(ediff) < energy_threshold_)
        return true;
    else
        return false;
}

void ROHF::form_initialF()
{
    // Form the initial Fock matrix, closed and open variants
    Fc_->copy(H_);
    Fo_->copy(H_);
    Fo_->scale(0.5);

    // Transform the Focks
    Fc_->transform(Shalf_);
    Fo_->transform(Shalf_);

#ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "Initial closed Fock matrix:\n");
        Fc_->print(outfile);
        fprintf(outfile, "Initial open Fock matrix:\n");
        Fo_->print(outfile);
    }
#endif
}

void ROHF::form_F()
{
    SharedMatrix Fct(factory_->create_matrix("Fock closed transformed"));
    SharedMatrix Fot(factory_->create_matrix("Fock open transformed"));

    // Form Fc_ and Fo_. See derivation notebook for equations.
    Fc_->copy(H_);
    Fc_->add(Gc_);
    Fo_->copy(H_);
    Fo_->scale(0.5);
    Fo_->add(Go_);

    // Transform Fc_ and Fo_ to MO basis
    Fct->transform(Fc_, Ca_);
    Fot->transform(Fo_, Ca_);

    // Form the effective Fock matrix, too
    // The effective Fock matrix has the following structure
    //  closed   | open     | virtual
    //  Fc         2(Fc-Fo)   Fc
    //  2(Fc-Fo)   Fc         2Fo
    //  Fc         2Fo         Fc
    int *opi = Fc_->rowspi();
    Feff_->copy(Fct);
    for (int h=0; h<Feff_->nirrep(); ++h) {
        for (int i=doccpi_[h]; i<doccpi_[h]+soccpi_[h]; ++i) {
            // Set the open/closed portion
            for (int j=0; j<doccpi_[h]; ++j) {
                double val = 2.0 * (Fct->get(h, i, j) - Fot->get(h, i, j));
                Feff_->set(h, i, j, val);
                Feff_->set(h, j, i, val);
            }
            // Set the open/virtual portion
            for (int j=doccpi_[h]+soccpi_[h]; j<opi[h]; ++j) {
                double val = 2.0 * Fot->get(h, i, j);
                Feff_->set(h, i, j, val);
                Feff_->set(h, j, i, val);
            }
            // Set the open/open portion
            for (int j=doccpi_[h]; j<doccpi_[h]+soccpi_[h]; ++j) {
                double val = Fot->get(h, i, j);
                Feff_->set(h, i, j, val);
                Feff_->set(h, j, i, val);
            }
        }
    }

    if (debug_) {
        Fc_->print(outfile);
        Fo_->print(outfile);
        Fct->print(outfile);
        Fot->print(outfile);
        Feff_->print(outfile);
    }
}

void ROHF::form_C()
{
    SharedMatrix temp(factory_->create_matrix());
    SharedMatrix eigvec(factory_->create_matrix());

    // Obtain new eigenvectors
    Feff_->diagonalize(eigvec, epsilon_a_);
    find_occupation();

    if (debug_) {
        fprintf(outfile, "In ROHF::form_C:\n");
        Feff_->print();
        eigvec->eivprint(epsilon_a_);
    }
    temp->gemm(false, false, 1.0, Ca_, eigvec, 0.0);
    Ca_->copy(temp);

    if (debug_) {
        Ca_->print(outfile);
    }
}

void ROHF::form_D()
{
    int h, i, j, m;
    int *opi = Dc_->rowspi();
    int nirreps = Dc_->nirrep();
    double val;
    for (h=0; h<nirreps; ++h) {
        for (i=0; i<opi[h]; ++i) {
            for (j=0; j<opi[h]; ++j) {
                val = 0.0;
                for (m=0; m<doccpi_[h]; ++m)
                    val += Ca_->get(h, i, m) * Ca_->get(h, j, m);
                Dc_->set(h, i, j, val);

                val = 0.0;
                for (m=doccpi_[h]; m<doccpi_[h]+soccpi_[h]; ++m)
                    val += Ca_->get(h, i, m) * Ca_->get(h, j, m);
                Do_->set(h, i, j, val);
            }
        }
    }

    if (debug_) {
        fprintf(outfile, "in ROHF::form_D:\n");
        Dc_->print(outfile);
        Do_->print(outfile);
    }
}

double ROHF::compute_initial_E()
{
    SharedMatrix Ho(factory_->create_matrix());
    Ho->copy(H_);
    Ho->scale(0.5);

    H_->print();
    Dc_->print();
    Do_->print();

    return nuclearrep_ + Dc_->vector_dot(H_) + Do_->vector_dot(Ho);
}

double ROHF::compute_E() {
    SharedMatrix HFc(factory_->create_matrix());
    HFc->copy(H_);
    HFc->add(Fc_);
    SharedMatrix HFo(factory_->create_matrix());
    HFo->copy(H_);
    HFo->scale(0.5);
    HFo->add(Fo_);
    double Etotal = nuclearrep_ + Dc_->vector_dot(HFc) + Do_->vector_dot(HFo);
    return Etotal;
}

void ROHF::form_G()
{
    /*
     * If we define
     *
     * Jpq = [Dc_rs + 0.5 Do_rs](pq|rs), Kc_pq = Dc_rs(pr|qs), and Ko_pq = 0.5 Do_rs(pr|qs)
     *
     * then the contributions we want are
     *
     * Gc = 2J - Kc - Ko, and Go = J - 0.5Kc - Ko
     *
     * So, if we temporarily scale Do, we can make this look just like a UHF G build.  Nice.
     *
     * Addendum: To support DF/Direct algorithms wlog, the functor call builds (effectively)
     *  Go_ = 2*J
     *  Ka  = Kc + 2Ko
     *  Kb  = Kc
     * These are then unwound to make Andy's canonical products above
     */
    shared_ptr<Matrix> Da = shared_ptr<Matrix>(factory_->create_matrix("Da"));
    shared_ptr<Matrix> Db = shared_ptr<Matrix>(factory_->create_matrix("Db"));
    shared_ptr<Matrix> Ka = shared_ptr<Matrix>(factory_->create_matrix("Ka"));
    shared_ptr<Matrix> Kb = shared_ptr<Matrix>(factory_->create_matrix("Kb"));
    Da->copy(Do_);
    Da->add(Dc_);
    Db->copy(Dc_);

    J_Ka_Kb_Functor jk_builder(Go_, Ka, Kb, Da, Db, Ca_, Cb_, nalphapi_, nbetapi_);
    process_tei<J_Ka_Kb_Functor>(jk_builder);

    Go_->scale(0.5);
    Kc_->copy(Kb);
    Ko_->copy(Ka);
    Ko_->scale(-1.0);
    Ko_->add(Kb);
    Ko_->scale(-0.5);

    //Go_->print();
    //Ka->print();
    //Kb->print();

    Kc_->scale(0.5);
    Go_->subtract(Kc_);
    Gc_->copy(Go_);
    Gc_->scale(2.0);
    Gc_->subtract(Ko_);
    Go_->subtract(Ko_);

}

}}
