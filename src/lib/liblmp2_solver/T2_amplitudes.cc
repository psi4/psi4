/*! \defgroup LMP2 lmp2: LMP2 Evaluation of Energy */

/*!
 ** \file
 ** \ingroup LMP2
 ** \LMP2 evaluation of energy
 */

#include <liblmp2_solver/lmp2.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>

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

int LMP2::compute_T2_energy(const int &iter) {
#ifdef HAVE_MADNESS

    int conv = 0;
    int terms = 0;

    if (iter > 1 && diis_ == 1) {
        for (int ij=0; ij < ij_pairs_; ij++) {
            if (me_ == ij_owner_[ij]) {
////                madness::Future<SharedMatrix> F_sum = task(me_, &LMP2::build_F_sum, ij, iter);
//                SharedMatrix F_sum = build_F_sum(ij, iter);
//                madworld_->taskq.fence();

                madness::Future<SharedMatrix> F_sum = task(me_, &LMP2::build_F_sum, ij, iter);

                madness::Future<SharedMatrix> Rtilde = task(me_, &LMP2::build_rtilde, F_sum, ij, iter);

                madness::Future<SharedMatrix> T2 = task(me_, &LMP2::amplitudes_T2, Rtilde, ij, iter);

                task(me_, &LMP2::store_T2, T2, ij);
            }
        }

        madworld_->taskq.fence();
        perform_diis(iter);

        for (int ij=0; ij < ij_pairs_; ij++) {
            terms += pair_domain_len_[ij] * pair_domain_len_[ij];
            if (me_ == ij_owner_[ij]) {

//                if (iter < diis_start_) {
//                    madness::Future<double> elmp2 = task(me_, &LMP2::energy, T2_amp_[div_][ij], ij, iter);
//                    task(me_, &LMP2::Elmp2_local_sum, elmp2);
//                    if (iter > 0) {
//                        madness::Future<double> drms = task(me_, &LMP2::T2_rms, T2_amp_[div_][ij],
//                                                            T2_amp_[dmat1_][ij], ij);
//                        task(me_, &LMP2::Drms_local_sum, drms);
//                    }
//                }
//                else {
//                    madness::Future<double> elmp2 = task(me_, &LMP2::energy, T2_ext_[nmat_][ij], ij, iter);
//                    task(me_, &LMP2::Elmp2_local_sum, elmp2);
//                    if (iter > 0) {
//                        madness::Future<double> drms = task(me_, &LMP2::T2_rms, T2_ext_[nmat_][ij],
//                                                            T2_ext_[omat_][ij], ij);
//                        task(me_, &LMP2::Drms_local_sum, drms);
//                    }

//                }

                madness::Future<double> elmp2;
                madness::Future<double> drms;
                if (iter < diis_start_) {
                    elmp2 = task(me_, &LMP2::energy, T2_amp_[div_][ij], ij, iter);
                    drms = task(me_, &LMP2::T2_rms, T2_amp_[div_][ij],
                                T2_amp_[dmat1_][ij], ij);
                }
                else {
                    elmp2 = task(me_, &LMP2::energy, T2_ext_[nmat_][ij], ij, iter);
                    drms = task(me_, &LMP2::T2_rms, T2_ext_[nmat_][ij],
                                T2_ext_[omat_][ij], ij);
                }

                task(me_, &LMP2::Elmp2_local_sum, elmp2);
                task(me_, &LMP2::Drms_local_sum, drms);

            }
        }
    }
    else {
        std::map<int, madness::Future<SharedMatrix> > F_sum;
        for (int ij=0; ij < ij_pairs_; ij++) {
            terms += pair_domain_len_[ij] * pair_domain_len_[ij];
            if (me_ == ij_owner_[ij]) {
                F_sum.insert(std::pair<int, madness::Future<SharedMatrix> >(ij, task(me_, &LMP2::build_F_sum, ij, iter)));
            }
        }
        for (int ij=0; ij < ij_pairs_; ij++) {
            if (me_ == ij_owner_[ij]) {

//                madness::Future<SharedMatrix> F_sum = task(me_, &LMP2::build_F_sum, ij, iter);

                madness::Future<SharedMatrix> Rtilde = task(me_, &LMP2::build_rtilde, F_sum[ij], ij, iter);

                madness::Future<SharedMatrix> T2 = task(me_, &LMP2::amplitudes_T2, Rtilde, ij, iter);

                madness::Future<double> elmp2 = task(me_, &LMP2::energy, T2, ij, iter);
                task(me_, &LMP2::Elmp2_local_sum, elmp2);
                if (iter > 0) {
                    madness::Future<double> drms = task(me_, &LMP2::T2_rms, T2_amp_[div_][ij],
                                                        T2_amp_[dmat1_][ij], ij);
                    task(me_, &LMP2::Drms_local_sum, drms);
                }
                task(me_, &LMP2::store_T2, T2, ij);

            }
        }
    }

    Communicator::world->sync();
    Communicator::world->sum(&Elmp2_, 1);

    if (iter > 0) {
        Communicator::world->sum(&Drms_T2_, 1);

        Drms_T2_ = sqrt(Drms_T2_/terms);
        Delta_Elmp2_ = Elmp2_ - Elmp2_old_;

        if (fabs(Delta_Elmp2_) < econv_ && fabs(Drms_T2_) < rmsconv_ || iter >= maxiter_) {
          conv = 1;
          if (iter >= maxiter_)
            if (me_ == 0)
              fprintf(outfile, "LMP2 has not converged in the maximum number of iterations.\n maxiter = %d\n", maxiter_);
        }
        else conv = 0;
    }
    return conv;
#else
    throw PSIEXCEPTION("PSI4 was not build with MADNESS.  Either "
                       "change your COMMUNICATOR environment variable, "
                       "or rebuild PSI4 with MADNESS\n");
#endif
}

#ifdef HAVE_MADNESS
//madness::Future<SharedMatrix> LMP2::get_old_T2(const int &ij) {
//    return madness::Future<SharedMatrix> (T2_amp_[dmat1_][ij]);
//}

SharedMatrix LMP2::get_old_T2(const int &ij) {
    return T2_amp_[dmat1_][ij];
}

//madness::Void LMP2::print_matrix(SharedMatrix mat)
//{
//    mat->print();
//    return madness::None;
//}

SharedMatrix LMP2::build_F_sum(const int &ij, const int &iter) {

    if (iter > 0) {
        madness::TaskAttributes attr;
        attr.set_highpriority(true);

        std::vector< madness::Future<SharedMatrix> > F_sum_kj = madness::future_vector_factory<SharedMatrix>(ndocc_);
        std::vector< madness::Future<SharedMatrix> > F_sum_ik = madness::future_vector_factory<SharedMatrix>(ndocc_);

        for (int k=0; k < ndocc_; k++) {
            int kj = ij_map_neglect_[ij_kj_map_[ij][k]];
            int ik = ij_map_neglect_[ij_ik_map_[ij][k]];

            if (pairdom_exist_[ij_kj_map_[ij][k]]) {
                F_sum_kj[k] = task(ij_owner_[kj], &LMP2::get_old_T2, kj, attr);
            }
            if (pairdom_exist_[ij_ik_map_[ij][k]]) {
                F_sum_ik[k] = task(ij_owner_[ik], &LMP2::get_old_T2, ik, attr);
            }
        }

        return task(me_, &LMP2::F_sum_add, F_sum_kj, F_sum_ik, ij);
    }
    else {
        return SharedMatrix(new Matrix());
    }
}

SharedMatrix LMP2::F_sum_add(const std::vector<madness::Future<SharedMatrix> > &T2_kj,
                             const std::vector<madness::Future<SharedMatrix> > &T2_ik,
                             const int &ij)
{
    SharedMatrix F_sum(new Matrix(nirreps_, &nso_, &nso_));
    double **f = F_sum->pointer();
    double **f_lo = F_LO_->pointer();

    int i = ij_i_map_[ij];
    int j = ij_j_map_[ij];

    for (int k=0; k < ndocc_; k++) {
        int kj = ij_map_neglect_[ij_kj_map_[ij][k]];
        int ik = ij_map_neglect_[ij_ik_map_[ij][k]];

        double **t_kj = T2_kj[k].get()->pointer();
        double **t_ik = T2_ik[k].get()->pointer();

        if (pairdom_exist_[ij_kj_map_[ij][k]]) {
            if (k <= j) {
                for (int r=0, L=0; r < natom_; r++) {
                    if (pair_domain_[kj][r]) {
                        for (int l = ao_start_[r]; l < ao_stop_[r]; l++, L++) {
                            for (int s=0, M=0; s < natom_; s++) {
                                if (pair_domain_[kj][s]) {
                                    for (int m = ao_start_[s]; m < ao_stop_[s]; m++, M++) {
                                        f[l][m] -= f_lo[i][k] * t_kj[M][L];
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else {
                for (int r=0, L=0; r < natom_; r++) {
                    if (pair_domain_[kj][r]) {
                        for (int l = ao_start_[r]; l < ao_stop_[r]; l++, L++) {
                            for (int s=0, M=0; s < natom_; s++) {
                                if (pair_domain_[kj][s]) {
                                    for (int m = ao_start_[s]; m < ao_stop_[s]; m++, M++) {
                                        f[l][m] -= f_lo[i][k] * t_kj[L][M];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (pairdom_exist_[ij_ik_map_[ij][k]]) {
            if (k >= i) {
                for (int r=0, L=0; r < natom_; r++) {
                    if (pair_domain_[ik][r]) {
                        for (int l = ao_start_[r]; l < ao_stop_[r]; l++, L++) {
                            for (int s=0, M=0; s < natom_; s++) {
                                if (pair_domain_[ik][s]) {
                                    for (int m = ao_start_[s]; m < ao_stop_[s]; m++, M++) {
                                        f[l][m] -= f_lo[k][j] * t_ik[M][L];
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else {
                for (int r=0, L=0; r < natom_; r++) {
                    if (pair_domain_[ik][r]) {
                        for (int l = ao_start_[r]; l < ao_stop_[r]; l++, L++) {
                            for (int s=0, M=0; s < natom_; s++) {
                                if (pair_domain_[ik][s]) {
                                    for (int m = ao_start_[s]; m < ao_stop_[s]; m++, M++) {
                                        f[l][m] -= f_lo[k][j] * t_ik[L][M];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return F_sum;
}

//madness::Void LMP2::F_sum_add_kj(SharedMatrix F_sum, const SharedMatrix T2_kj,
//                                 const int &k, const int &ij)
//{
//    F_mutex_->lock();
//        double **f = F_sum->pointer();
//        double **t = T2_kj->pointer();
//        double **f_lo = F_LO_->pointer();
//    F_mutex_->unlock();


//    int i = ij_i_map_[ij];
//    int j = ij_j_map_[ij];
//    int kj = ij_map_neglect_[ij_kj_map_[ij][k]];

//    if (k <= j) {
//        for (int r=0, L=0; r < natom_; r++) {
//            if (pair_domain_[kj][r]) {
//                for (int l = ao_start_[r]; l < ao_stop_[r]; l++, L++) {
//                    for (int s=0, M=0; s < natom_; s++) {
//                        if (pair_domain_[kj][s]) {
//                            for (int m = ao_start_[s]; m < ao_stop_[s]; m++, M++) {
//                                f[l][m] -= f_lo[i][k] * t[M][L];
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    else {
//        for (int r=0, L=0; r < natom_; r++) {
//            if (pair_domain_[kj][r]) {
//                for (int l = ao_start_[r]; l < ao_stop_[r]; l++, L++) {
//                    for (int s=0, M=0; s < natom_; s++) {
//                        if (pair_domain_[kj][s]) {
//                            for (int m = ao_start_[s]; m < ao_stop_[s]; m++, M++) {
////                                F_sum->add(0, l, m, -(F_LO_->get(0,i,k) * T2_kj->get(0,L,M)));
//                                f[l][m] -= f_lo[i][k] * t[L][M];
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    return madness::None;
//}
//madness::Void LMP2::F_sum_add_ik(SharedMatrix F_sum, const SharedMatrix T2_ik,
//                                 const int &k, const int &ij)
//{
//    F_mutex_->lock();
//        double **f = F_sum->pointer();
//        double **t = T2_ik->pointer();
//        double **f_lo = F_LO_->pointer();
//    F_mutex_->unlock();

//    int i = ij_i_map_[ij];
//    int j = ij_j_map_[ij];
//    int ik = ij_map_neglect_[ij_ik_map_[ij][k]];

//    if (k >= i) {
//        for (int r=0, L=0; r < natom_; r++) {
//            if (pair_domain_[ik][r]) {
//                for (int l = ao_start_[r]; l < ao_stop_[r]; l++, L++) {
//                    for (int s=0, M=0; s < natom_; s++) {
//                        if (pair_domain_[ik][s]) {
//                            for (int m = ao_start_[s]; m < ao_stop_[s]; m++, M++) {
//                                f[l][m] -= f_lo[k][j] * t[M][L];
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    else {
//        for (int r=0, L=0; r < natom_; r++) {
//            if (pair_domain_[ik][r]) {
//                for (int l = ao_start_[r]; l < ao_stop_[r]; l++, L++) {
//                    for (int s=0, M=0; s < natom_; s++) {
//                        if (pair_domain_[ik][s]) {
//                            for (int m = ao_start_[s]; m < ao_stop_[s]; m++, M++) {
////                                F_sum->add(0, l, m, -(F_LO_->get(0,k,j) * T2_ik->get(0,L,M)));
//                                f[l][m] -= f_lo[k][j] * t[L][M];
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    return madness::None;
//}


SharedMatrix LMP2::build_rtilde(const SharedMatrix F_sum, const int &ij, const int &iter) {

    int ij_loc = ij_local_[ij];
    SharedMatrix T2(new Matrix(pair_domain_len_[ij], pair_domain_len_[ij]));
    T2->copy(eri_ij_[ij]);

    if(iter > 0) {

        SharedMatrix Temp1 = SharedMatrix(T2->clone());

        Temp1->gemm(false, false, 1.0, F_virt_[ij_loc], T2_amp_[dmat1_][ij], 0.0);
        T2->gemm(false, false, 1.0, Temp1, S_virt_[ij_loc], 1.0);

        Temp1->gemm(false, false, 1.0, S_virt_[ij_loc], T2_amp_[dmat1_][ij], 0.0);
        T2->gemm(false, false, 1.0, Temp1, F_virt_[ij_loc], 1.0);

        Temp1->zero();

        double **temp1 = block_matrix(pair_domain_len_[ij], nso_);
        double **Sij = block_matrix(pair_domain_len_[ij], nso_);

        int L=0;
        for (int r = 0; r < natom_; r++) {
            if (pair_domain_[ij][r]) {
                for (int l = ao_start_[r]; l < ao_stop_[r]; l++, L++) {
                    for (int m = 0; m < nso_; m++) {
                        Sij[L][m] = S_->get(0,l,m);
                    }
                }
            }
        }

        C_DGEMM('n', 'n', pair_domain_len_[ij], nso_, nso_, 1, &(Sij[0][0]),
                nso_, &(F_sum->pointer(0)[0][0]), nso_, 0, &(temp1[0][0]), nso_);
        C_DGEMM('n', 't', pair_domain_len_[ij], pair_domain_len_[ij], nso_, 1, &(temp1[0][0]),
                nso_, &(Sij[0][0]), nso_, 0, &(Temp1->pointer(0)[0][0]), pair_domain_len_[ij]);

        free_block(temp1);
        free_block(Sij);

//        T2->subtract(temp1);  ?????
        T2->add(Temp1);

    }

    return T2;

}

SharedMatrix LMP2::amplitudes_T2(const SharedMatrix T2, const int &ij, const int &iter) {

    std::stringstream ij_val;
    int ij_loc = ij_local_[ij];
    int i = ij_i_map_[ij];
    int j = ij_j_map_[ij];

    if (print_ > 3) {
        ij_val << ij;
    }

    if (print_ > 4) {
        T2->set_name("Rtilde[" + ij_val.str() + "]");
        T2->print();
    }

    T2->transform(W_[ij_loc]);
    T2->set_name("Rbar[" + ij_val.str() + "]");

    SharedMatrix Denominator = SharedMatrix(new Matrix("Denominator[" + ij_val.str() + "]",
                                                       nirreps_, &pair_domain_nr_len_[ij],
                                                       &pair_domain_nr_len_[ij]));

    for (int a = 0; a < pair_domain_nr_len_[ij]; a++) {
        for (int b = 0; b < pair_domain_nr_len_[ij]; b++) {
            Denominator->set( 0, a, b,
                             (evals_[ij_loc]->get(0, a) + evals_[ij_loc]->get(0, b) -
                              F_LO_->get(0,i,i) - F_LO_->get(0,j,j)) );
        }
    }
    T2->scale(-1.0);
    T2->apply_denominator(Denominator);

    if (print_ > 4) {
        T2->set_name("Tbar[" + ij_val.str() + "]");
        T2->print();
    }

    T2->back_transform(W_[ij_loc]);

    if (iter > 0)
        T2->add(T2_amp_[dmat1_][ij]);

    if (print_ > 3) {
        T2->set_name("T2 Amplitudes[" + ij_val.str() + "]");
        T2->print();
    }

    return T2;

}
#endif

//void amplitudes_T2_old(const int &iter) {

//    std::vector< std::vector<int> > ij_kj_map(lmp2_info.ij_pairs_, std::vector<int>(lmp2_info.ndocc_, 0));
//    std::vector< std::vector<int> > ij_ik_map(lmp2_info.ij_pairs_, std::vector<int>(lmp2_info.ndocc_, 0));

//    double ****Replicate_kj = (double ****) malloc(lmp2_info.ij_pairs_per_proc_ * sizeof (double ***));
//    double ****Replicate_ik = (double ****) malloc(lmp2_info.ij_pairs_per_proc_ * sizeof (double ***));
//    SharedMatrix Denominator = SharedMatrix(new Matrix("Denominator"));

//    // Set up the ij to kj/ik maps
//    for (int ij = 0; ij < lmp2_info.ij_pairs_; ij++) {
//        int i = lmp2_info.ij_i_map_[ij];
//        int j = lmp2_info.ij_j_map_[ij];

//        for (int k = 0; k < lmp2_info.ndocc_; k++) {
//            if (k > j)
//                ij_kj_map[ij][k] = (k * (k + 1)) / 2 + j;
//            else
//                ij_kj_map[ij][k] = (j * (j + 1)) / 2 + k;

//            if (i > k)
//                ij_ik_map[ij][k] = (i * (i + 1)) / 2 + k;
//            else
//                ij_ik_map[ij][k] = (k * (k + 1)) / 2 + i;
//        }
//    }

//    timer_on("T2 Replication");

//    // This will replicate the required amplitudes
//    if(iter > 0) {
//        for (int ij = 0; ij < lmp2_info.ij_pairs_; ij++) {
//            int i = lmp2_info.ij_i_map_[ij];
//            int j = lmp2_info.ij_j_map_[ij];

//            if (lmp2_info.me_ == lmp2_info.ij_owner_[ij]) {
//                Replicate_kj[lmp2_info.ij_local_[ij]] = (double ***) malloc(lmp2_info.ndocc_ * sizeof (double ***));
//                Replicate_ik[lmp2_info.ij_local_[ij]] = (double ***) malloc(lmp2_info.ndocc_ * sizeof (double ***));
//            }
//            // Get the required amplitudes from their respective proc and put them in a temporary array
//            // If ij == ik/kj then do DCOPY
//            for (int k = 0; k < lmp2_info.ndocc_; k++) {

//                int kj = lmp2_info.ij_map_neglect_[ij_kj_map[ij][k]];
//                int ik = lmp2_info.ij_map_neglect_[ij_ik_map[ij][k]];

//                if (lmp2_info.pairdom_exist_[ij_kj_map[ij][k]]) {
//                    if ((lmp2_info.me_ == lmp2_info.ij_owner_[ij]) & (lmp2_info.me_ == lmp2_info.ij_owner_[kj])) {
//                        Replicate_kj[lmp2_info.ij_local_[ij]][k] = block_matrix(lmp2_info.pair_domain_len_[kj], lmp2_info.pair_domain_len_[kj]);
//                        C_DCOPY(lmp2_info.pair_domain_len_[kj] * lmp2_info.pair_domain_len_[kj], lmp2_info.T2_[lmp2_info.dmat1_]->get_pointer(lmp2_info.ij_local_[kj]), 1,
//                                &(Replicate_kj[lmp2_info.ij_local_[ij]][k][0][0]), 1);
//                    }
//                    else if ((lmp2_info.me_ == lmp2_info.ij_owner_[ij]) & (lmp2_info.me_ != lmp2_info.ij_owner_[kj])) {
//                        Replicate_kj[lmp2_info.ij_local_[ij]][k] = block_matrix(lmp2_info.pair_domain_len_[kj], lmp2_info.pair_domain_len_[kj]);
//                        Communicator::World->recv(&(Replicate_kj[lmp2_info.ij_local_[ij]][k][0][0]), lmp2_info.pair_domain_len_[kj] * lmp2_info.pair_domain_len_[kj], lmp2_info.ij_owner_[kj]);
//                    }
//                    else if ((lmp2_info.me_ != lmp2_info.ij_owner_[ij]) & (lmp2_info.me_ == lmp2_info.ij_owner_[kj])) {
//                        Communicator::World->task(&(lmp2_info.T2_[lmp2_info.dmat1_]->pointer(lmp2_info.ij_local_[kj])[0][0]), lmp2_info.pair_domain_len_[kj] * lmp2_info.pair_domain_len_[kj], lmp2_info.ij_owner_[ij]);
//                    }
//                }

//                if (ik != kj) {
//                    if (lmp2_info.pairdom_exist_[ij_ik_map[ij][k]]) {
//                        if ((lmp2_info.me_ == lmp2_info.ij_owner_[ij]) & (lmp2_info.me_ == lmp2_info.ij_owner_[ik])) {
//                            Replicate_ik[lmp2_info.ij_local_[ij]][k] = block_matrix(lmp2_info.pair_domain_len_[ik], lmp2_info.pair_domain_len_[ik]);
//                            C_DCOPY(lmp2_info.pair_domain_len_[ik] * lmp2_info.pair_domain_len_[ik], lmp2_info.T2_[lmp2_info.dmat1_]->get_pointer(lmp2_info.ij_local_[ik]), 1,
//                                    &(Replicate_ik[lmp2_info.ij_local_[ij]][k][0][0]), 1);
//                        }
//                        else if ((lmp2_info.me_ == lmp2_info.ij_owner_[ij]) & (lmp2_info.me_ != lmp2_info.ij_owner_[ik])) {
//                            Replicate_ik[lmp2_info.ij_local_[ij]][k] = block_matrix(lmp2_info.pair_domain_len_[ik], lmp2_info.pair_domain_len_[ik]);
//                            Communicator::World->recv(&(Replicate_ik[lmp2_info.ij_local_[ij]][k][0][0]), lmp2_info.pair_domain_len_[ik] * lmp2_info.pair_domain_len_[ik], lmp2_info.ij_owner_[ik]);

//                        }
//                        else if ((lmp2_info.me_ != lmp2_info.ij_owner_[ij]) & (lmp2_info.me_ == lmp2_info.ij_owner_[ik])) {
//                            Communicator::World->task(&(lmp2_info.T2_[lmp2_info.dmat1_]->pointer(lmp2_info.ij_local_[ik])[0][0]), lmp2_info.pair_domain_len_[ik] * lmp2_info.pair_domain_len_[ik], lmp2_info.ij_owner_[ij]);
//                        }
//                    }
//                }
//            }
//        }
//        Communicator::World->sync();

//        // This copies rep_kj to rep_ik when kj == ik, this reduces communication
//        for (int ij = 0; ij < lmp2_info.ij_pairs_; ij++) {
//            int i = lmp2_info.ij_i_map_[ij];
//            int j = lmp2_info.ij_j_map_[ij];

//            for (int k = 0; k < lmp2_info.ndocc_; k++) {
//                int kj = lmp2_info.ij_map_neglect_[ij_kj_map[ij][k]];
//                int ik = lmp2_info.ij_map_neglect_[ij_ik_map[ij][k]];

//                if (lmp2_info.me_ == lmp2_info.ij_owner_[ij]) {
//                    if (kj == ik) {
//                        if (lmp2_info.pairdom_exist_[ij_ik_map[ij][k]]) {
//                            Replicate_ik[lmp2_info.ij_local_[ij]][k] = block_matrix(lmp2_info.pair_domain_len_[ik], lmp2_info.pair_domain_len_[ik]);

//                            C_DCOPY(lmp2_info.pair_domain_len_[ik] * lmp2_info.pair_domain_len_[ik], &(Replicate_kj[lmp2_info.ij_local_[ij]][k][0][0]), 1,
//                                    &(Replicate_ik[lmp2_info.ij_local_[ij]][k][0][0]), 1);
//                        }
//                    }
//                }
//            }
//        }
//    }

//    timer_off("T2 Replication");


//    timer_on("T2 Computation");

//    lmp2_info.T2_[lmp2_info.div_]->copy(lmp2_info.Ktilde_);
//    if (iter > 0) {

//        SharedMatrix Temp1 = SharedMatrix (lmp2_info.T2_[lmp2_info.div_]->clone());

//        Temp1->gemm(false, false, 1.0, lmp2_info.F_virt_, lmp2_info.T2_[lmp2_info.dmat1_], 0.0);
//        lmp2_info.T2_[lmp2_info.div_]->gemm(false, false, 1.0, Temp1, lmp2_info.S_virt_, 1.0);

//        Temp1->gemm(false, false, 1.0, lmp2_info.S_virt_, lmp2_info.T2_[lmp2_info.dmat1_], 0.0);
//        lmp2_info.T2_[lmp2_info.div_]->gemm(false, false, 1.0, Temp1, lmp2_info.F_virt_, 1.0);

//        if (lmp2_info.print_ > 5) {
//            lmp2_info.T2_[lmp2_info.div_]->set_name("Rtilde after transforms");
//            lmp2_info.T2_[lmp2_info.div_]->print();
//        }

//        Temp1->zero();

//        // Build S*sum[F(ik)T(kj) + F(kj)T(ik)]S
//        for (int ij = 0; ij < lmp2_info.ij_pairs_; ij++) {
//            int i = lmp2_info.ij_i_map_[ij];
//            int j = lmp2_info.ij_j_map_[ij];

//            if (lmp2_info.me_ != lmp2_info.ij_owner_[ij])
//                continue;

//            double **F_sum = block_matrix(lmp2_info.nso_, lmp2_info.nso_);


//            for (int k = 0; k < lmp2_info.ndocc_; k++) {
//                int kj = lmp2_info.ij_map_neglect_[ij_kj_map[ij][k]];
//                int ik = lmp2_info.ij_map_neglect_[ij_ik_map[ij][k]];

//                if (lmp2_info.pairdom_exist_[ij_kj_map[ij][k]]) {
//                    for (int r=0, L=0; r < lmp2_info.natom_; r++) {
//                        if (lmp2_info.pair_domain_[kj][r]) {
//                            for (int l = lmp2_info.ao_start_[r]; l < lmp2_info.ao_stop_[r]; l++, L++) {
//                                for (int s=0, M=0; s < lmp2_info.natom_; s++) {
//                                    if (lmp2_info.pair_domain_[kj][s]) {
//                                        for (int m = lmp2_info.ao_start_[s]; m < lmp2_info.ao_stop_[s]; m++, M++) {
//                                            if (k <= j)
//                                                F_sum[l][m] -= lmp2_info.F_LO_->get(0,i,k) * Replicate_kj[lmp2_info.ij_local_[ij]][k][M][L];
//                                            else
//                                                F_sum[l][m] -= lmp2_info.F_LO_->get(0,i,k) * Replicate_kj[lmp2_info.ij_local_[ij]][k][L][M];
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//                if (lmp2_info.pairdom_exist_[ij_ik_map[ij][k]]) {
//                    for (int r = 0, L=0; r < lmp2_info.natom_; r++) {
//                        if (lmp2_info.pair_domain_[ik][r]) {
//                            for (int l = lmp2_info.ao_start_[r]; l < lmp2_info.ao_stop_[r]; l++, L++) {
//                                for (int s = 0, M=0; s < lmp2_info.natom_; s++) {
//                                    if (lmp2_info.pair_domain_[ik][s]) {
//                                        for (int m = lmp2_info.ao_start_[s]; m < lmp2_info.ao_stop_[s]; m++, M++) {
//                                            if (k >= i)
//                                                F_sum[l][m] -= lmp2_info.F_LO_->get(0,k,j) * Replicate_ik[lmp2_info.ij_local_[ij]][k][M][L];
//                                            else
//                                                F_sum[l][m] -= lmp2_info.F_LO_->get(0,k,j) * Replicate_ik[lmp2_info.ij_local_[ij]][k][L][M];
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }

//            if (lmp2_info.print_ > 6 & lmp2_info.me_ == 0) {
//                fprintf(outfile, "\tF_sum ij: %d\n", ij);
//                print_mat(F_sum, lmp2_info.nso_, lmp2_info.nso_, outfile);
//            }
//            double **temp1 = block_matrix(lmp2_info.pair_domain_len_[ij], lmp2_info.nso_);
//            double **Sij = block_matrix(lmp2_info.pair_domain_len_[ij], lmp2_info.nso_);

//            int L=0;
//            for (int r = 0; r < lmp2_info.natom_; r++) {
//                if (lmp2_info.pair_domain_[ij][r]) {
//                    for (int l = lmp2_info.ao_start_[r]; l < lmp2_info.ao_stop_[r]; l++, L++) {
//                        for (int m = 0; m < lmp2_info.nso_; m++) {
//                            Sij[L][m] = lmp2_info.S_->get(0,l,m);
//                        }
//                    }
//                }
//            }

//            if (lmp2_info.print_ > 5 & lmp2_info.me_ == 0) {
//                fprintf(outfile, "\tSij: %d\n", ij);
//                print_mat(Sij, lmp2_info.pair_domain_len_[ij], lmp2_info.nso_, outfile);

//                fprintf(outfile, "\tF_sum: %d\n", ij);
//                print_mat(F_sum, lmp2_info.nso_, lmp2_info.nso_, outfile);
//            }

//            C_DGEMM('n', 'n', lmp2_info.pair_domain_len_[ij], lmp2_info.nso_, lmp2_info.nso_, 1, &(Sij[0][0]),
//                    lmp2_info.nso_, &(F_sum[0][0]), lmp2_info.nso_, 0, &(temp1[0][0]), lmp2_info.nso_);
//            C_DGEMM('n', 't', lmp2_info.pair_domain_len_[ij], lmp2_info.pair_domain_len_[ij], lmp2_info.nso_, 1, &(temp1[0][0]),
//                    lmp2_info.nso_, &(Sij[0][0]), lmp2_info.nso_, 0, Temp1->get_pointer(lmp2_info.ij_local_[ij]), lmp2_info.pair_domain_len_[ij]);

//            free_block(temp1);
//            free_block(Sij);
//            free_block(F_sum);

//        }
//        if (lmp2_info.print_ > 5) {
//            Temp1->set_name("Xt");
//            Temp1->print();
//        }

//        lmp2_info.T2_[lmp2_info.div_]->add(Temp1);
//    }
//    if (lmp2_info.print_ > 4) {
//        lmp2_info.T2_[lmp2_info.div_]->set_name("Rtilde");
//        lmp2_info.T2_[lmp2_info.div_]->print();
//    }

//    lmp2_info.T2_[lmp2_info.div_]->transform(lmp2_info.W_);
//    lmp2_info.T2_[lmp2_info.div_]->set_name("Rbar");

//    Denominator->copy(lmp2_info.T2_[lmp2_info.div_]);
//    Denominator->set_name("Denominator");
//    for (int ij = 0; ij < lmp2_info.ij_pairs_; ij++) {
//        if (lmp2_info.me_ != lmp2_info.ij_owner_[ij])
//            continue;

//        int i = lmp2_info.ij_i_map_[ij];
//        int j = lmp2_info.ij_j_map_[ij];

//        for (int a = 0; a < lmp2_info.pair_domain_nr_len_[lmp2_info.ij_local_[ij]]; a++) {
//            for (int b = 0; b < lmp2_info.pair_domain_nr_len_[lmp2_info.ij_local_[ij]]; b++) {
//                Denominator->set( lmp2_info.ij_local_[ij], a, b,
//                        (lmp2_info.evals_->get(lmp2_info.ij_local_[ij], a) + lmp2_info.evals_->get(lmp2_info.ij_local_[ij], b) -
//                         lmp2_info.F_LO_->get(0,i,i) - lmp2_info.F_LO_->get(0,j,j)) );
//            }
//        }
//    }
//    lmp2_info.T2_[lmp2_info.div_]->scale(-1.0);
//    lmp2_info.T2_[lmp2_info.div_]->apply_denominator(Denominator);
//    if (lmp2_info.print_ > 4) {
//        lmp2_info.T2_[lmp2_info.div_]->set_name("Tbar");
//        lmp2_info.T2_[lmp2_info.div_]->print();
//    }

//    lmp2_info.T2_[lmp2_info.div_]->back_transform(lmp2_info.W_);

//    if (iter > 0)
//        lmp2_info.T2_[lmp2_info.div_]->add(lmp2_info.T2_[lmp2_info.dmat1_]);

//    if (lmp2_info.print_ > 3) {
//        lmp2_info.T2_[lmp2_info.div_]->set_name("T2 Amplitudes");
//        lmp2_info.T2_[lmp2_info.div_]->print();
//    }

//    // Free the memory used by the Replicated T2 Amplitudes
//    if (iter > 0) {
//        for (int ij=0; ij < lmp2_info.ij_pairs_per_proc_; ij++) {
//            for (int k = 0; k < lmp2_info.ndocc_; k++) {
//                free_block(Replicate_kj[ij][k]);
//                free_block(Replicate_ik[ij][k]);
//            }
//        free(Replicate_kj[ij]);
//        free(Replicate_ik[ij]);
//        }
//    }
//    free(Replicate_kj);
//    free(Replicate_ik);

//    timer_off("T2 Computation");



//}

#endif // have_madnesss

}} // End of psi and lmp2 namespace
