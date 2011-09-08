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


void LMP2::integral_direct_transformation() {
    set_mn_integral_maps();
    setup_eri_matrices();

//    Integrals_ = distributed_container(*Communicator::world->get_madworld());
//    Integrals_.clear();
//    distributed_container::iterator it;

//    for (int M=0, MN=0; M < nshell_; M++) {
//        // We need to do some screening here
//        for (int N=0; N < nshell_; N++, MN++) {
//            // We need to do some screening here
//            if (me_ == Integrals_.owner(dist_key(MN))) {
//                Integrals_.replace(dist_key(MN), dist_container());

//                it = Integrals_.find(dist_key(MN));
//                it->second.init_wfn(wfn_, basisset_, molecule_, eri_);
//                it->second.init_MN(MN, ndocc_, C_, 0.0);
//            }
//        }
//    }
//    // Make sure everything is initialized
//    madworld_->gop.fence();

//    for (int M=0, MN=0; M < nshell_; M++) {
//        // We need to do some screening here
//        for (int N=0; N < nshell_; N++, MN++) {
//            // We need to do some screening here
//            if (me_ == Integrals_.owner(dist_key(MN))) {
//                Integrals_.task(dist_key(MN), &dist_container::allocate_RS_block);
//            }
//        }
//    }
//    // Make sure the RS blocks for each MN pair are set up
//    madworld_->gop.fence();


//    for (int M=0, MN=0; M < nshell_; M++) {
//        for (int N=0; N < nshell_; N++, MN++) {
//            if (me_ == Integrals_.owner(dist_key(MN))) {
//                Integrals_.task(dist_key(MN), &dist_container::first_half_transformation, M, N);
//            }
//        }
//    }
//    // Need to make sure that all of the integrals have been half transformed
//    madworld_->gop.fence();

//    for (int M=0, MN=0; M < nshell_; M++) {
//        for (int N=0; N < nshell_; N++, MN++) {
//            if (me_ == 0) {
//                it = Integrals_.find(dist_key(MN));
//                it->second.print_mn();
//            }
//        }
//    }




    timer_on("First_Half_Transformation");

    int MN=0;
    for (MN_iter_ = MN_Pairs_.begin(); MN_iter_ != MN_Pairs_.end(); MN_iter_++, MN++) {

        if (comm_ == "MADNESS") {
#ifdef HAVE_MADNESS
            if (me_ == MN_Owner_[MN]) {
                task(me_, &LMP2::first_half_integral_transformation,
                     MN_iter_->second[0], MN_iter_->second[2],
                     MN_iter_->second[1], MN_iter_->second[3], MN);
            }

#else
            throw PSIEXCEPTION("PSI4 was not built with MADNESS. "
                               "Please change your COMMUNICATOR env to "
                               "MPI or LOCAL, or rebuild PSI4 with MADNESS.\n");
#endif
        }
        else {
            if (me_ == MN_Owner_[MN]) {
                first_half_integral_transformation(MN_iter_->second[0], MN_iter_->second[2],
                                                   MN_iter_->second[1], MN_iter_->second[3], MN);
            }
        }
    }

    Communicator::world->sync();


    timer_off("First_Half_Transformation");


    redistribute_integrals();

    Communicator::world->sync();


    eri_2_MN_.clear();

//    for (int ij=0; ij < ij_pairs_; ij++) {
//        task(ij_owner_[ij], &LMP2::print_eri, ij, "eri_2");
//    }

    Communicator::world->sync();


    timer_on("Second_Half_Transformation");

    for (int ij=0; ij < ij_pairs_; ij++) {
        if (me_ == ij_owner_[ij]) {
            task(me_, &LMP2::second_half_transformation, ij);
        }
    }
    Communicator::world->sync();


//    for (int ij=0; ij < ij_pairs_; ij++) {
//        task(ij_owner_[ij], &LMP2::print_eri, ij, "k_tilde");
//    }

    Communicator::world->sync();

    timer_off("Second_Half_Transformation");


}

int LMP2::print_eri(const int &ij, std::string name) {
    if (me_ == ij_owner_[ij]) {
        std::stringstream ij_val;
        ij_val << ij;

        name += "[";
        name += ij_val.str();
        name += "]";

        eri_ij_[ij]->set_name(name.c_str());
        eri_ij_[ij]->print();
    }
}

void LMP2::set_mn_integral_maps() {
    timer_on("MN_shell_sort");
    // This creates a map of MN pairs in decreasing AM size
    for (int M = 0; M < nshell_; M++) {
        std::vector<int> Shell(4,0);
        Shell[0] = M;
        Shell[1] = basisset_->shell(M)->nfunction();
        for (int N = 0; N < nshell_; N++) {
            Shell[2] = N;
            Shell[3] = basisset_->shell(N)->nfunction();
            MN_Pairs_.insert(Shell_Pair(Shell[1]+Shell[3], Shell));
        }
        Shell.clear();
    }
    timer_off("MN_shell_sort");

    maxshell_ = MN_Pairs_.begin()->second[1];

    mn_pairs_per_proc_ = 0;
    for (int MN=0, local=0; MN < MN_Pairs_.size(); MN++) {
        MN_Owner_.push_back(MN%nproc_);

        if (me_ == MN_Owner_[MN])
            MN_local_.insert(std::pair<int,int>(MN, local));
        if (MN%nproc_ == nproc_-1) local++;
    }
    mn_pairs_per_proc_ = MN_local_.size();

    madworld_->gop.fence();


}

void LMP2::setup_eri_matrices() {
    for (int MN=0; MN < MN_Pairs_.size(); MN++)  {
        // We need screeing here.  This may need to be pushed into first_half_transformation to do the screening.
        if (me_ == MN_Owner_[MN])
            eri_2_MN_.insert(std::pair<int,SharedMatrix>(MN, SharedMatrix(new Matrix())) );
    }
    for (int ij=0; ij < ij_pairs_; ij++) {
        // We need screeing here.  This may need to be pushed into first_half_transformation to do the screening.
        if (me_ == ij_owner_[ij])
            eri_ij_.insert(std::pair<int,SharedMatrix>(ij,SharedMatrix(new Matrix(nso_,nso_))));
    }
}

int LMP2::first_half_integral_transformation(const int &M, const int &N,
                                             const int &numm, const int &numn,
                                             const int &MN)
{

    int thread = Communicator::world->thread_id( pthread_self() );
    int mn_size = numm*numn;

    int mn_loc = MN_local_[MN];
    std::stringstream mn_val;
    mn_val << MN;

    const double *buffer = eri_[thread]->buffer();

    int occ_squared = ndocc_ * ndocc_;
    eri_2_MN_[MN]->init(nirreps_, &occ_squared, &mn_size);

    std::vector< SharedMatrix > eri_1;
    std::stringstream occ_val;
    for (int i=0; i < ndocc_; i++) {
        occ_val << i;
        eri_1.push_back( boost::shared_ptr<Matrix>(
                            new Matrix(nirreps_, &mn_size, &nso_)) );
    }

//    int max_per_shell = basisset_->max_function_per_shell();
//    SharedMatrix temp_buf(new Matrix("temp buff", max_per_shell*max_per_shell*max_per_shell, nso_));
//    double **Tp = temp_buf->pointer();

//    SharedMatrix temp2_buf(new Matrix("temp buff", max_per_shell*max_per_shell*max_per_shell, ndocc_));
//    double **T2p = temp2_buf->pointer();

//    std::vector< SharedMatrix > temp_eri_ij;
//    for (int i=0; i < ndocc_; i++) {
//        occ_val << i;
//        temp_eri_ij.push_back( SharedMatrix(
//                            new Matrix(nirreps_, &ndocc_, &mn_size)) );
//    }



    // **** First Quarter Integral Transformation
    for (int R = 0; R < nshell_; R++) {
        int numr = basisset_->shell(R)->nfunction();
        int abs_r = basisset_->shell(R)->function_index();
        int r_size = abs_r+numr;

        for (int S = 0; S < nshell_; S++) {
            int nums = basisset_->shell(S)->nfunction();
            int abs_s = basisset_->shell(S)->function_index();

            eri_[thread]->compute_shell(M, R, N, S);

//            for (int m=0, index=0; m<numm; m++) {
//                for (int r=0; r<numr; r++) {
//                    for (int n=0; n<numn; n++) {
//                        for (int s=0; s<nums; s++, index++) {
//                            Tp[m*numn*numr+n*numr+r][s+abs_s] = buffer[index];
//                        }
//                    }
//                }
//            }
//        } // end of S

//        if (numm*numn*numr > 30) {
//            C_DGEMM('N','N', numm*numn*numr, ndocc_, nso_, 1.0,
//                    Tp[0], nso_, C_->pointer()[0], ndocc_, 0.0, T2p[0], ndocc_);
//        }
//        else {
//            for (int mnr=0; mnr < numm*numn*numr; mnr++) {
//                C_DGEMV('T',nso_,ndocc_,1.0,C_->pointer()[0],ndocc_,
//                        Tp[mnr], 1, 0.0, T2p[mnr], 1);
//            }
//        }

//        for (int m=0, index=0; m<numm; m++) {
//            for (int n=0; n<numn; n++) {
//                for (int r=0; r<numr; r++, index++) {
//                    for (int j=0; j < ndocc_; j++) {
//                        eri_1[j]->set(0,m*numn +n, r+abs_r, T2p[index][j]);
//                    }
//                }
//            }
//        }

            // **** First Quarter Integral Transformation **** //
            for (int j=0; j < ndocc_; j++) {
                for (int m=0, index=0; m < numm; m++) {
                    for (int r=abs_r; r < r_size; r++) {
                        for (int n=0; n < numn; n++, index+=nums) {
                            int mn = m*numn + n;
                            eri_1[j]->add(0, mn, r, C_DDOT(nums, &(C_->pointer(0)[abs_s][j]), nso_, const_cast<double*>(&(buffer[index])), 1));
                        }
                    }
                }
            }
        } // End of S loop
    } // End of R loop

    if (print_ > 6) {
#ifdef HAVE_MADNESS
        if (comm_ == "MADNESS") {
            print_mutex->lock();
        }
#endif
        for (int i=0; i < ndocc_; i++) {
            eri_1[i]->set_name("eri_1[" + occ_val.str() + "]");
            fprintf(outfile, "eri_1[%d]", i);
            eri_1[i]->print();
        }
#ifdef HAVE_MADNESS
        if (comm_ == "MADNESS") {
            print_mutex->unlock();
        }
#endif
    }

    // **** Second Quarter Integral Transformation **** //
//    for (int i=0; i < ndocc_; i++) {
//        C_DGEMM('T','T',ndocc_, numm*numn, nso_, 1.0, C_->pointer()[0], ndocc_,
//                eri_1[i]->pointer()[0], nso_, 0.0, temp_eri_ij[i]->pointer()[0], numm*numn);
//    }
//    // unpack
//    for (int ij = 0; ij < ij_pairs_; ij++) {
//        int i = ij_i_map_[ij];
//        int j = ij_j_map_[ij];
//        for (int mn=0; mn < mn_size; mn++) {
//            eri_2_MN_[mn_loc]->set(0, ij, mn, temp_eri_ij[i]->get(0,j,mn));
//        }
//    }
    for (int ij = 0; ij < ij_pairs_; ij++) {
        int i = ij_i_map_[ij];
        int j = ij_j_map_[ij];
        for (int mn=0; mn < mn_size; mn++) {
            eri_2_MN_[MN]->add( 0, ij, mn, C_DDOT(nso_, &(C_->pointer(0)[0][i]), nso_, &(eri_1[j]->pointer(0)[mn][0]), 1) );
        }
    }

    eri_1.clear();

    if (print_ > 6) {
#ifdef HAVE_MADNESS
        if (comm_ == "MADNESS") {
            print_mutex->lock();
        }
#endif
        eri_2_MN_[MN]->set_name("ERI_2_MN[" + mn_val.str() + "]");
        eri_2_MN_[MN]->print();
#ifdef HAVE_MADNESS
        if (comm_ == "MADNESS") {
            print_mutex->unlock();
        }
#endif
    }

    return 0;
}

int LMP2::second_half_transformation(const int &ij) {

    SharedMatrix Rt = SharedMatrix(new Matrix(nirreps_, &nso_, &pair_domain_len_[ij]));

    for (int k=0, b=0; k < natom_; k++) {
        if (pair_domain_[ij][k]) {
            for (int l=ao_start_[k]; l < ao_stop_[k]; l++, b++) {
                C_DCOPY(nso_, &(Rt_full_->pointer(0)[0][l]),
                        nso_, &(Rt->pointer(0)[0][b]),
                        pair_domain_len_[ij]);
            }
        }
    }

    eri_ij_[ij]->transform(Rt);

    return 0;
}

std::vector<double> LMP2::copy_eri_mn(const int &MN, const int &ij,
                                      const int &numm, const int &numn)
{
    std::vector<double> temp(numm*numn);
    C_DCOPY(numm*numn, &(eri_2_MN_[MN]->pointer(0)[ij][0]), 1, &(temp[0]), 1);
    return temp;
}

#ifdef HAVE_MADNESS
madness::Future< std::vector<double> > LMP2::copy_eri_mn_mad(const int &MN, const int &ij,
                                      const int &numm, const int &numn)
{
    std::vector<double> temp(numm*numn);
    C_DCOPY(numm*numn, &(eri_2_MN_[MN]->pointer(0)[ij][0]), 1, &(temp[0]), 1);
    return madness::Future< std::vector<double> > (temp);
}

madness::Void LMP2::redist_madness_integrals(const int &MN, const int &ij,
                                             const int &numm, const int &numn,
                                             const int &abs_m, const int &abs_n)
{
    task(ij_owner_[ij], &LMP2::redist_integrals, copy_eri_mn(MN, ij, numm, numn), ij, numm, numn, abs_m, abs_n);
    return madness::None;
}
#endif

int LMP2::redist_integrals(const std::vector<double> &eri_mn,
                           const int &ij, const int &numm, const int &numn,
                           const int &abs_m, const int &abs_n)
{
    for (int m=abs_m, mn=0; m < abs_m+numm; m++) {
        for (int n=abs_n; n < abs_n+numn; n++, mn++) {
            eri_ij_[ij]->set(0, m, n, eri_mn[mn]);
        }
    }
    return 0;
}




void LMP2::redistribute_integrals()
{

    int MN=0;
    for (MN_iter_ = MN_Pairs_.begin(); MN_iter_ != MN_Pairs_.end(); MN_iter_++, MN++) {
        int M = MN_iter_->second[0];
        int numm = MN_iter_->second[1];
        int N = MN_iter_->second[2];
        int numn = MN_iter_->second[3];
        int abs_m = basisset_->shell(M)->function_index();
        int abs_n = basisset_->shell(N)->function_index();

        if (comm_ == "MADNESS") {
#ifdef HAVE_MADNESS

            if (me_ == MN_Owner_[MN]) {
                for (int ij=0; ij < ij_pairs_; ij++) {
                    task(MN_Owner_[MN], &LMP2::redist_madness_integrals, MN, ij, numm, numn, abs_m, abs_n);
                }
            }

#else
            throw PSIEXCEPTION("PSI4 was not built with MADNESS. "
                               "Please change your COMMUNICATOR env to "
                               "MPI or LOCAL, or rebuild PSI4 with MADNESS.\n");
#endif
        }
        else {
            // This compiles but I haven't tested it yet
            for (int ij=0; ij < ij_pairs_; ij++) {
                std::vector<double> eri_mn;
                eri_mn.resize(numm*numn);
                memset(&eri_mn[0], 0, numm*numn*sizeof(double));

                if (me_ == MN_Owner_[MN] && me_ == ij_owner_[ij]) {
                    eri_mn = copy_eri_mn(MN, ij, numm, numn);
                    redist_integrals(eri_mn, ij, numm, numn, abs_m, abs_n);
                }
                else if (me_ != MN_Owner_[MN] && me_ == ij_owner_[ij]) {
                    Communicator::world->recv(&(eri_mn[0]), numm*numn, MN_Owner_[MN]);
                    redist_integrals(eri_mn, ij, numm, numn, abs_m, abs_n);
                }
                else if (me_ == MN_Owner_[MN] && me_ != ij_owner_[ij]) {
                    Communicator::world->send(&((copy_eri_mn(MN, ij, MN_iter_->second[1],
                                                            MN_iter_->second[3]))[0]),
                                              numm*numn, ij_owner_[ij]);
                }
            }
        }
    }

}

#endif // have madness

}}
