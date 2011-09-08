/*! \defgroup LMP2 lmp2: LMP2 Evaluation of Energy */

/*!
 ** \file
 ** \ingroup LMP2
 ** \LMP2 evaluation of energy
 */

#include <liblmp2_solver/lmp2.h>
#include <libqt/qt.h>


namespace boost {
template<class T> class shared_ptr;
}

namespace psi{

class BasisSet;
class Options;
class PSIO;
class Chkpt;

namespace lmp2 {

#ifdef HAVE_MADNESS

void LMP2::projection()
{

    timer_on("Build Projectors");
    for (int ij=0, v=0; ij < ij_pairs_; ij++, v++) {
        if (me_ == v%nproc_) {
            if (comm_ == "MADNESS") {
#ifdef HAVE_MADNESS
                task(me_, &LMP2::build_projection_matrix, ij);
//                madworld_->taskq.add(*this, &LMP2::build_projection_matrix, ij);
#else
                throw PSIEXCEPTION("PSI4 was not built with MADNESS. "
                                   "Please change your COMMUNICATOR env to "
                                   "MPI or LOCAL, or rebuild PSI4 with MADNESS.\n");
#endif
            }
            else {
                build_projection_matrix(ij);
            }
        }
    }

    Communicator::world->sync();

    timer_off("Build Projectors");

}

void LMP2::init_global_mat()
{
    Rt_full_ = boost::shared_ptr<Matrix>(nso_nso_->create_matrix("Full Virtual Space Projector"));
    for (int i=0; i < ij_pairs_per_proc_; i++) {
        S_virt_.push_back( boost::shared_ptr<Matrix> (new Matrix()) );
        F_virt_.push_back( boost::shared_ptr<Matrix> (new Matrix()) );
        W_.push_back( boost::shared_ptr<Matrix> (new Matrix()) );
        evals_.push_back( boost::shared_ptr<Vector> (new Vector()) );
    }

    pair_domain_nr_len_.resize(ij_pairs_);

    /* Compute the complete virtual space projector */
    // This is the only place the D is used, so we can get D from the wavefunction in the future
    Rt_full_->identity();
    Rt_full_->gemm(false,false,-1.0, D_AO_, S_, 1.0);

    /* Transform the overlap and fock matrices from the AO to PAO basis */
    S_->transform(Rt_full_);
    F_AO_->transform(Rt_full_);

    S_->set_name("Overlap Matrix (PAO)");
    F_AO_->set_name("Fock Matrix (PAO)");

    if (print_ >= 2) {
        Rt_full_->print();
        S_->print();
        F_AO_->print();
    }

}



int LMP2::build_projection_matrix(const int &ij)
{

    int ij_loc = ij_local_[ij];
    std::stringstream ij_val;
    ij_val << ij;


    // Build the virtual overlap and fock matrices
    S_virt_[ij_loc]->init(nirreps_, &pair_domain_len_[ij], &pair_domain_len_[ij]); //, "S_virt[" + ij_val.str() + "]");
    F_virt_[ij_loc]->init(nirreps_, &pair_domain_len_[ij], &pair_domain_len_[ij]); //, "F_virt[" + ij_val.str() + "]");

    int K=0;
    for (int r=0; r < natom_; r++) {
        if (pair_domain_[ij][r]) {
            for (int k=ao_start_[r]; k < ao_stop_[r]; k++,K++) {
                int L=0;
                for (int s=0; s < natom_; s++) {
                    if (pair_domain_[ij][s]) {
                        for (int l=ao_start_[s]; l < ao_stop_[s]; l++,L++) {
                            S_virt_[ij_loc]->set(0, K, L, S_->get(0, k, l));
                            F_virt_[ij_loc]->set(0, K, L, F_AO_->get(0, k, l));
                        }
                    }
                }
            }
        }
    }

    SharedMatrix evec_st = SharedMatrix(new Matrix("evec_st[" + ij_val.str() + "]",
                                                   nirreps_, &pair_domain_len_[ij],
                                                   &pair_domain_len_[ij]));
    SharedVector eval = SharedVector(new Vector("eval[" + ij_val.str() + "]",
                                                nirreps_, &pair_domain_len_[ij]));

    S_virt_[ij_loc]->diagonalize(evec_st, eval, static_cast<Matrix::DiagonalizeOrder>(1) );

    if (print_ > 5) {
#ifdef HAVE_MADNESS
        if (comm_ == "MADNESS") {
            print_mutex->lock();
        }
#endif

        S_virt_[ij_loc]->print();
        F_virt_[ij_loc]->print();
        evec_st->eivprint(eval);

#ifdef HAVE_MADNESS
        if (comm_ == "MADNESS") {
            print_mutex->unlock();
        }
#endif
    }

    // Build Xt
    pair_domain_nr_len_[ij] = pair_domain_len_[ij];
    for (int k=0; k < eval->dim(0); k++) {
        if (eval->get(0,k) <= 1e-6) {
            pair_domain_nr_len_[ij]--;
        }
    }

    // These are needed to build the projection matrix
    SharedMatrix Xt = SharedMatrix(new Matrix("Xt[" + ij_val.str() + "]",
                                              nirreps_, &pair_domain_len_[ij],
                                              &pair_domain_nr_len_[ij]));
    SharedMatrix evec_fbar = SharedMatrix(new Matrix("evec_fbar[" + ij_val.str() + "]",
                                                     nirreps_, &pair_domain_nr_len_[ij],
                                                     &pair_domain_nr_len_[ij]));

    // Build Xt
    int I=0;
    for (int k=0; k < pair_domain_len_[ij]; k++) {
        if (eval->get(0,k) > 1e-6) {
            for (int l=0; l < pair_domain_len_[ij]; l++)
                Xt->set(0, l, I, (evec_st->get(0,l,k)/sqrt(eval->get(0,k))));
            I++;
        }
    }
    if (print_ > 5) {
#ifdef HAVE_MADNESS
        if (comm_ == "MADNESS") {
            print_mutex->lock();
        }
#endif

        Xt->print();

#ifdef HAVE_MADNESS
        if (comm_ == "MADNESS") {
            print_mutex->unlock();
        }
#endif
    }

    evals_[ij_loc]->init(nirreps_, &pair_domain_nr_len_[ij], "Evals[" + ij_val.str() + "]");
    SharedMatrix F_bar = boost::shared_ptr<Matrix> (F_virt_[ij_loc]->clone());

    F_bar->transform(Xt);
    F_bar->diagonalize(evec_fbar, evals_[ij_loc]); //, static_cast<Matrix::DiagonalizeOrder>(1));

    if (print_ > 5) {
#ifdef HAVE_MADNESS
        if (comm_ == "MADNESS") {
            print_mutex->lock();
        }
#endif

        F_bar->set_name("F_bar[" + ij_val.str() + "]");
        F_bar->eivprint(evals_[ij_loc]);
        evals_[ij_loc]->print();

#ifdef HAVE_MADNESS
        if (comm_ == "MADNESS") {
            print_mutex->unlock();
        }
#endif
    }

    W_[ij_loc]->init(nirreps_, &pair_domain_len_[ij], &pair_domain_nr_len_[ij],
                     "W[" + ij_val.str() + "]");
    W_[ij_loc]->gemm(false, false, 1.0, Xt, evec_fbar, 0.0);

    if (print_ > 4) {
#ifdef HAVE_MADNESS
        if (comm_ == "MADNESS") {
            print_mutex->lock();
        }
#endif
        W_[ij_loc]->print();
#ifdef HAVE_MADNESS
        if (comm_ == "MADNESS") {
            print_mutex->unlock();
        }
#endif
    }

    fflush(outfile);

    return 0;

}

#endif // have_madness

}}
