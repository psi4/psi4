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
    Fc_        = SharedMatrix(factory_.create_matrix("F closed"));
    Fo_        = SharedMatrix(factory_.create_matrix("F open"));
    Fa_        = SharedMatrix(factory_.create_matrix("F effective (MO basis)"));
    Feff_      = Fa_;
    Ca_        = SharedMatrix(factory_.create_matrix("Moleular orbitals"));
    Cb_        = Ca_;
    Dc_        = SharedMatrix(factory_.create_matrix("D closed"));
    Do_        = SharedMatrix(factory_.create_matrix("D open"));
    Dc_old_    = SharedMatrix(factory_.create_matrix("D closed old"));
    Do_old_    = SharedMatrix(factory_.create_matrix("D open old"));
    Gc_        = SharedMatrix(factory_.create_matrix("G closed"));
    Go_        = SharedMatrix(factory_.create_matrix("G open"));
    epsilon_a_ = SharedVector(factory_.create_vector());
    epsilon_b_ = epsilon_a_;

    pk_ = NULL;
    k_ = NULL;

    fprintf(outfile, "  DIIS %s.\n\n", diis_enabled_ ? "enabled" : "disabled");

    if (scf_type_ == "PK")
        allocate_PK();
}

void ROHF::finalize()
{
    if (pk_)
        delete[](pk_);
    if (k_)
        delete[](k_);
   
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

#ifdef _DEBUG
    Ca_->print(outfile, "initial C");
#endif
}

double ROHF::compute_energy()
{
    bool converged = false, diis_iter=false;
    int iter = 0;

    // Do the initial work to give the iterations a starting point.
    form_H();

    if (scf_type_ == "PK")
        form_PK();
    //else if (scf_type_ == "DF")
    //    form_B();

    //if (scf_type_ == "PK")
    //    form_PK();

    form_Shalf();
    // Check to see if there are MOs already in the checkpoint file.
    // If so, read them in instead of forming them.
    if (load_or_compute_initial_C())
        fprintf(outfile, "  Read in previous MOs from file32.\n\n");

    fprintf(outfile, "                                  Total Energy            Delta E              Density RMS\n\n");
    do {
        iter++;

        Dc_old_->copy(Dc_); // save previous density
        Do_old_->copy(Do_); // save previous density
        Eold_ = E_; // save previous energy

        if (scf_type_ == "PK")
            form_G_from_PK();
        else{

        }
        //else if (scf_type_ == "DIRECT")
        //    form_G_from_direct_integrals();
        //else if (scf_type_ == "DF" || scf_type_ == "CD" || scf_type_ == "1C_CD")
        //    form_G_from_RI();
        //else if (scf_type_ == "OUT_OF_CORE")
        //    form_G();

        form_F(); // Forms: Fc_, Fo_, Feff_

        if (diis_enabled_)
            save_fock(); // Save the effective Fock for diis

        // Compute total energy
        E_ = compute_E();

        if (diis_enabled_ == true && iter >= min_diis_vectors_ && iter % 6 == 0) {
            diis();
            diis_iter = true;
        } else {
            diis_iter = false;
        }
        fprintf(outfile,
                "  @ROHF iteration %3d energy: %20.14f    %20.14f %s\n",
                iter, E_, E_ - Eold_, diis_iter == false ? " " : "DIIS");
        fflush(outfile);

        form_C(); 	// Uses Feff_ to form C_.
        //	find_occupation(_F);
        form_D();

        converged = test_convergency();
    } while (!converged && iter < maxiter_);
    //if (scf_type_ == "DF" || scf_type_ == "CD" || scf_type_ == "1C_CD")
    //{
    //    free_B();
    //}

    // Return the final ROHF energy
    if (converged) {
        fprintf(outfile, "\n  Energy converged.\n");
        save_information();
    } else {
        fprintf(outfile, "\n  Failed to converge.\n");
        E_ = 0.0;
    }
    
    finalize();

    return E_;
}

void ROHF::save_information()
{
    // Print the final docc vector
    char **temp2 = molecule_->irrep_labels();

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

void ROHF::allocate_PK()
{
    // Figure out how many pair combinations yield A1 symmetry (done in above loop)
    //   num_pair_combinations_of_A1 = ioff[_opi[0]] + ioff[_opi[1]] + ioff[_opi[2]] + ...
    // Allocate memory for the PK matrix (using a vector)
    if (pk_size_ < (memory_ / sizeof(double) / 2)) {
        pk_ = new double[pk_size_];
        k_ = new double[pk_size_];

        if (pk_ == NULL || k_ == NULL) {
            fprintf(outfile, "  Insufficient free system memory for in-core PK implementation.\n");
            fprintf(outfile, "  Switching to out-of-core algorithm.\n");
            scf_type_ = "OUT_OF_CORE";
        } else {
            // Zero out PK and K
            memset(pk_, 0, pk_size_*sizeof(double));
            memset(k_, 0, pk_size_*sizeof(double));

            fprintf(outfile,
                "  Allocated %lu elements (%lu pairs) for PK. (%5f MiB)\n",
                (unsigned long)pk_size_, (unsigned long)pk_pairs_, pk_size_ * 8.0 / 1048576.0);
            fprintf(outfile,
                "  Allocated %lu elements (%lu pairs) for K.  (%5f MiB)\n\n",
                (unsigned long)pk_size_, (unsigned long)pk_pairs_, pk_size_ * 8.0 / 1048576.0);
        }
    } else {
        fprintf(outfile,
                "  Insufficient memory for in-core PK implementation.\n");
        fprintf(outfile,
                "  Would need %lu elements of double memory. (%5f MiB)\n",
                (unsigned long)pk_size_*2, pk_size_ * 8.0 / 1048576.0 * 2.0);
        fprintf(outfile, "  Switching to out-of-core algorithm.\n");
        scf_type_ = "OUT_OF_CORE";
    }
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
#ifdef _DEBUG
    if (debug_) {
        Fc_->print(outfile);
        Fo_->print(outfile);
        Fct->print(outfile);
        Fot->print(outfile);
        Feff_->print(outfile);
    }
#endif
}

void ROHF::form_C()
{
    SharedMatrix temp(factory_->create_matrix());
    SharedMatrix eigvec(factory_->create_matrix());

    // Obtain new eigenvectors
    Feff_->diagonalize(eigvec, epsilon_a_);
    find_occupation();

#ifdef _DEBUG
    if (debug_) {
        eigvec->eivprint(epsilon_a_);
    }
#endif
    temp->gemm(false, false, 1.0, Ca_, eigvec, 0.0);
    Ca_->copy(temp);

#ifdef _DEBUG
    if (debug_) {
        Ca_->print(outfile);
    }
#endif
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

#ifdef _DEBUG
    if (debug_) {
        Dc_->print(outfile);
        Do_->print(outfile);
    }
#endif
}

double ROHF::compute_initial_E() {
    SharedMatrix Ho(factory_->create_matrix());
    Ho->copy(H_);
    Ho->scale(0.5);

    double Etotal = nuclearrep_ + Dc_->vector_dot(H_) + Do_->vector_dot(Ho);
    fprintf(outfile, "\n  Initial ROHF energy: %20.14f\n\n", Etotal);
    fflush(outfile);
    return Etotal;
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

void ROHF::form_PK() {
    // struct iwlbuf ERIIN;
    int ilsti, nbuf;
    int i, j, k, l;
    int ii, jj, kk, ll;
    int is, js, ks, ls;
    int fi;
    size_t bra, ket, braket;
    int idx;
    int counter = 0;
    double value;

    // PK zeroed out during allocation
    fprintf(outfile, "  Forming PK and K matrices.\n");
    fflush(outfile);

    IWL ERIIN(psio_.get(), PSIF_SO_TEI, 0.0, 1, 1);

    do {
        ilsti = ERIIN.last_buffer();
        nbuf = ERIIN.buffer_count();

        fi = 0;
        for (idx=0; idx<nbuf; ++idx) {
            i = ERIIN.labels()[fi] > 0 ? ERIIN.labels()[fi] : -ERIIN.labels()[fi];
            j = ERIIN.labels()[fi+1];
            k = ERIIN.labels()[fi+2];
            l = ERIIN.labels()[fi+3];
            value = ERIIN.values()[idx];
            fi += 4;

            // Get the symmetries
            is = so2symblk_[i];
            js = so2symblk_[j];
            ks = so2symblk_[k];
            ls = so2symblk_[l];

            // Get the offset of the SO index in its symblock
            ii = so2index_[i];
            jj = so2index_[j];
            kk = so2index_[k];
            ll = so2index_[l];

            // J
            if ((is == js) && (ks == ls)) {
                bra = INDEX2(ii, jj);
                ket = INDEX2(kk, ll);
                // pk_symoffset_ corrects for the symmetry offset in the pk_ vector
                braket = INDEX2(bra + pk_symoffset_[is], ket + pk_symoffset_[ks]);
                pk_[braket] += value;

                // K/2 (2nd sort)
                if ((ii != jj) && (kk != ll)) {
                    if ((is == ls) && (js == ks)) {
                        bra = INDEX2(ii, ll);
                        ket = INDEX2(jj, kk);
                        braket = INDEX2(bra + pk_symoffset_[is], ket + pk_symoffset_[js]);
                        if ((ii == ll) || (jj == kk)) {
                            pk_[braket] -= 0.5 * value;
                            k_[braket] -= 0.5 * value;
                        } else {
                            pk_[braket] -= 0.25 * value;
                            k_[braket] -= 0.25 * value;
                        }
                    }
                }
            }

            // K/2 (1st sort)
            if ((is == ks) && (js == ls)) {
                bra = INDEX2(ii, kk);
                ket = INDEX2(jj, ll);
                braket = INDEX2(bra + pk_symoffset_[is], ket + pk_symoffset_[js]);
                if ((ii == kk) || (jj == ll)) {
                    pk_[braket] -= 0.5 * value;
                    k_[braket] -= 0.5 * value;
                } else {
                    pk_[braket] -= 0.25 * value;
                    k_[braket] -= 0.25 * value;
                }
            }
            counter++;
        }

        if (!ilsti)
            ERIIN.fetch();
    } while (!ilsti);

    // Going out of scope will close the buffer
    // iwl_buf_close(&ERIIN, 1);

    // After stage two is complete, the elements of P must be halved for the case IJ=KL.
    for (size_t ij=0; ij < pk_pairs_; ++ij) {
        pk_[INDEX2(ij,ij)] *= 0.5;
        k_[INDEX2(ij,ij)] *= 0.5;
    }

    fprintf(outfile, "  Processed %d two-electron integrals.\n", counter);
#ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "pk_:\n");
        print_array(pk_, pk_pairs_, outfile);
        fprintf(outfile, "k_:\n");
        print_array(k_, pk_pairs_, outfile);
    }
#endif
}

void ROHF::form_G_from_PK()
{
    int nirreps = factory_->nirrep();
    int *opi = factory_->rowspi();
    size_t ij;
    double *Do_vector = new double[pk_pairs_];
    double *Dc_vector = new double[pk_pairs_];
    double *Gc_vector = new double[pk_pairs_];
    double *Go_vector = new double[pk_pairs_];

    Gc_->zero();
    Go_->zero();

    memset(Do_vector, 0, sizeof(double) * pk_pairs_);
    memset(Dc_vector, 0, sizeof(double) * pk_pairs_);
    memset(Gc_vector, 0, sizeof(double) * pk_pairs_);
    memset(Go_vector, 0, sizeof(double) * pk_pairs_);

    ij=0;
    for (int h=0; h<nirreps; ++h) {
        for (int p=0; p<opi[h]; ++p) {
            for (int q=0; q<=p; ++q) {
                if (p != q) {
                    Dc_vector[ij] = 2.0 * Dc_->get(h, p, q);
                    Do_vector[ij] = 2.0 * Do_->get(h, p, q);
                } else {
                    Dc_vector[ij] = Dc_->get(h, p, q);
                    Do_vector[ij] = Do_->get(h, p, q);
                }
                ij++;
            }
        }
    }

#ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "PK: ij = %lu\n", (unsigned long)ij);
        fflush(outfile);
        Dc_->print(outfile);
        fprintf(outfile, "PK: Dc vector:\n");
        for (ij=0; ij<pk_pairs_; ++ij)
            fprintf(outfile, "PK: Dc vector [%lu] = %20.14f\n", (unsigned long)ij, Dc_vector[ij]);
        Do_->print(outfile);
        fprintf(outfile, "PK: Do vector:\n");
        for (ij=0; ij<pk_pairs_; ++ij)
            fprintf(outfile, "PK: Do vector [%lu] = %20.14f\n", (unsigned long)ij, Do_vector[ij]);
    }
#endif

    /*
     * This code goes through the densities (Dc_ and Do_), PK, and K to form
     * two G matrices. One G matrix is for Fc_ and the other for Fo_.
     * See derivation notebook for equations.
     */
    double Gc_pq, Dc_pq;
    double Go_pq, Do_pq;
    double* Dc_rs;
    double* Gc_rs;
    double* Do_rs;
    double* Go_rs;
    int pq, rs;
    double* PK_block = pk_;
    double* K_block = k_;
    int ts_pairs = pk_pairs_;
    for (pq = 0; pq < ts_pairs; ++pq) {
        Gc_pq = 0.0;
        Dc_pq = Dc_vector[pq];
        Dc_rs = &Dc_vector[0];
        Gc_rs = &Gc_vector[0];
        Go_pq = 0.0;
        Do_pq = Do_vector[pq];
        Do_rs = &Do_vector[0];
        Go_rs = &Go_vector[0];
        for (rs = 0; rs <= pq; ++rs) {
            // D_{rs}^{c} * PK_{pqrs}         Also found in RHF
            Gc_pq  += *PK_block * (*Dc_rs);
            *Gc_rs += *PK_block * Dc_pq;
            // D_{rs}^{o} * PK_{pqrs} / 2     Yes, open D adds to closed G
            Gc_pq  += *PK_block * (*Do_rs) * 0.5;
            *Gc_rs += *PK_block * Do_pq    * 0.5;
            // D_{rs}^{c} * PK_{pqrs} / 2     Yes, closed D adds to open G
            Go_pq  += *PK_block * (*Dc_rs) * 0.5;
            *Go_rs += *PK_block * Dc_pq    * 0.5;
            // D_{rs}^{o} * (PK_{pqrs} + K_{pqrs}) / 4
            Go_pq  += (*PK_block + *K_block) * (*Do_rs) * 0.25;
            *Go_rs += (*PK_block + *K_block) * Do_pq    * 0.25;
            ++Dc_rs;
            ++Gc_rs;
            ++Do_rs;
            ++Go_rs;
            ++PK_block;
            ++K_block;
        }
        Gc_vector[pq] += Gc_pq;
        Go_vector[pq] += Go_pq;
    }

    // Convert G to a matrix
    ij = 0;
    for (int h = 0; h < nirreps; ++h) {
        for (int p = 0; p < opi[h]; ++p) {
            for (int q = 0; q <= p; ++q) {
                Gc_->set(h, p, q, 2.0 * Gc_vector[ij]);
                Gc_->set(h, q, p, 2.0 * Gc_vector[ij]);
                Go_->set(h, p, q, 2.0 * Go_vector[ij]);
                Go_->set(h, q, p, 2.0 * Go_vector[ij]);
                ij++;
            }
        }
    }

#ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "Gc from PK:\n");
        Gc_->print(outfile);
        fprintf(outfile, "Go from PK:\n");
        Go_->print(outfile);
    }
#endif

    delete[](Dc_vector);
    delete[](Do_vector);
    delete[](Gc_vector);
    delete[](Go_vector);
}


}}
