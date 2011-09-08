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

#ifdef HAVE_MADNESS

madness::Future<int> LMP2::build_error(const int &ij)
{
//    error_[dmat2_][ij]->zero();
    error_[dmat2_][ij]->copy(T2_amp_[div_][ij]);
    error_[dmat2_][ij]->subtract(T2_amp_[dmat1_][ij]);
    return madness::Future<int>(0);
}

madness::Void LMP2::build_Bmat(const int &check, const int &ij)
{
    int thread = Communicator::world->thread_id(pthread_self());
    for(int p=0; p < matsize_; p++) {
        for(int q=0; q < matsize_; q++) {
            Bmat_[thread]->add(0, p, q, error_[p][ij]->vector_dot(error_[q][ij]));
        }
    }
    return madness::None;
}

#endif

void LMP2::perform_diis(const int &iter) {

    if (iter < diis_start_) {
        for (int ij=0; ij < ij_pairs_; ij++) {
            if (me_ == ij_owner_[ij]) {
                task(me_, &LMP2::build_error, ij);
            }
        }
        Communicator::world->sync();
    }
    else {
        // *** Compute the B matrix ***
        int bmat_size = matsize_+1;

        for (int i=0; i < nthread_; i++) {
            Bmat_.push_back(SharedMatrix(new Matrix("B matrix", nirreps_, &bmat_size, &bmat_size)));
            Bmat_[i]->zero();
        }

        for (int ij=0; ij < ij_pairs_; ij++) {
            if (me_ == ij_owner_[ij]) {
                madness::Future<int> check = task(ij_owner_[ij], &LMP2::build_error, ij);
                task(me_, &LMP2::build_Bmat, check, ij);
            }
        }

        Communicator::world->sync();

        for (int i=1; i < nthread_; i++) {
            Bmat_[0]->add(Bmat_[i]);
        }

        Communicator::world->sum(&(Bmat_[0]->pointer(0)[0][0]), bmat_size*bmat_size);

        double **bmat = Bmat_[0]->pointer();
        for (int i=0; i < matsize_; i++) {
            bmat[matsize_][i] = bmat[i][matsize_] = -1.0;
        }

//        Communicator::world->sync();
//        if (me_ == 0) Bmat[0]->print();
//        Communicator::world->sync();

        Communicator::world->bcast(&(Bmat_[0]->pointer(0)[0][0]), bmat_size*bmat_size);


        std::vector<double> co(matsize_+1);
        co[matsize_] = -1.0;

        std::vector<int> work(matsize_+1);

        // use DGESV to compute the solution to Ax=B
        C_DGESV(matsize_+1, 1, &(bmat[0][0]), matsize_+1, &(work[0]), &(co[0]), matsize_+1);

        // compute the new MP2 coefficients
        for (int ij = 0; ij < ij_pairs_; ij++) {
            if (me_ == ij_owner_[ij]) {
                task(ij_owner_[ij], &LMP2::build_T2_ext, ij, co);
            }
        }

        Communicator::world->sync();

        Bmat_.clear();
        co.clear();
        work.clear();

    }

}


madness::Void LMP2::build_T2_ext(const int &ij, const std::vector<double> &co) {
    T2_ext_[nmat_][ij]->zero();
    for (int a=0; a < pair_domain_len_[ij]; a++) {
        for (int b=0; b < pair_domain_len_[ij]; b++) {
            for (int k=0; k < matsize_; k++) {
                int c_ext;
                if (k == max_diis_vectors_ - 1) c_ext = k-k + 1;
                else if (k == max_diis_vectors_ - 2) c_ext = k - k;
                else c_ext = k + 2;
                T2_ext_[nmat_][ij]->add(0, a, b, (co[k] * T2_amp_[c_ext][ij]->get(0,a,b)));
            }
        }
    }
    return madness::None;
}

#endif // have_madness

}} // End of psi and LMP2 naglob.mespaces
