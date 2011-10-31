/*! \defgroup LMP2 lmp2: LMP2 Evaluation of Energy */

/*!
 ** \file
 ** \ingroup LMP2
 ** \LMP2 evaluation of energy
 */

#include "dist_container.h"
#include <sstream>
#include <string>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{

#ifdef HAVE_MADNESS


typedef std::map<int, std::vector<double> > block;
typedef std::pair<int, std::vector<double> > shell_block;
typedef block::iterator shell_block_it;

dist_container::dist_container() :
    me_(Communicator::world->me()), nproc_(Communicator::world->nproc()),
    nthread_(Communicator::world->nthread()),
    comm_(Communicator::world->communicator()),
    exists_(false)
{ }

void dist_container::init_wfn(boost::shared_ptr<Wavefunction> wfn,
                              boost::shared_ptr<BasisSet> basisset,
                              boost::shared_ptr<Molecule> molecule,
                              std::vector<boost::shared_ptr<TwoElectronInt> > &eri)
{
    wfn_ = wfn;
    basisset_ = basisset;
    molecule_ = molecule;
    eri_ = eri;
}

void dist_container::init_MN(const int &MN, const int &ndocc,
                            SharedMatrix C, const double bound)
{
    ndocc_ = ndocc;
    MN_ = MN;
    thresh_ = bound;
    C_ = C;
    exists_ = true;
}



madness::Void dist_container::allocate_RS_block()
{
    if (exists_) {
        // Run through all R and S for a given M and N
        for (int R=0, RS=0; R < basisset_->nshell(); R++) {
            // We need to do some screening here
            for (int S=0; S < basisset_->nshell(); S++, RS++) {
                // We need to do some screening here
                RS_map_.insert(std::pair<int, std::pair<int,int> >(RS, std::pair<int,int>(R,S)));
            }
        }
    }
    return madness::None;
}

void dist_container::print_mn()
{
    if (exists_) {
        shell_block_it it;
        fprintf(outfile, "\n *** MN (%d) Owner (%d) ***\n", MN_, me_);
        for (it = eri_2_mn_.begin(); it != eri_2_mn_.end(); it++) {
            fprintf (outfile,"\nRS %5d",it->first);
            for (int i=0; i < it->second.size(); i++) {
                fprintf (outfile,"\n%5d",i);
                fprintf (outfile,"%12.7f",it->second[i]);
            }
            fprintf (outfile,"\n");
        }
    }
    else {
        fprintf(outfile, "\n *** MN (%d) Owner (%d) ***\n", MN_, me_);
        fprintf(outfile, "Does Not Exist\n");
    }
}

madness::Void dist_container::first_half_transformation(const int &M, const int &N)
{
    if (exists_) {

        int thread = Communicator::world->thread_id( pthread_self() );
        const double *buffer = eri_[thread]->buffer();


        int numm = basisset_->shell(M)->nfunction();
        int numn = basisset_->shell(N)->nfunction();
        int mn_size = numm*numn;
        int nirreps = wfn_->nirrep();
        int nso = wfn_->nso();

        std::map<int,SharedMatrix> eri_1;

        std::map<int,std::pair<int,int> >::iterator it;
        for (it = RS_map_.begin(); it != RS_map_.end(); it++) {
            int R = it->second.first;
            int S = it->second.second;

            int numr = basisset_->shell(R)->nfunction();
            int nums = basisset_->shell(S)->nfunction();

            int abs_r = basisset_->shell(R)->function_index();
            int abs_s = basisset_->shell(S)->function_index();

            int r_size = abs_r+numr;

            eri_[thread]->compute_shell(M, R, N, S);

            // **** First Quarter Integral Transformation **** //
            for (int j=0; j < ndocc_; j++) {
                // We need some screening here
                eri_1.insert(std::pair<int,SharedMatrix>(j,SharedMatrix(new Matrix(nirreps, &mn_size, &nso))));
            }

            std::map<int,SharedMatrix>::iterator it;
            for (it = eri_1.begin(); it != eri_1.end(); it++) {
                int j = it->first;
                for (int m=0, index=0; m < numm; m++) {
                    for (int r=abs_r; r < r_size; r++) {
                        for (int n=0; n < numn; n++, index+=nums) {
                            int mn = m*numn + n;
                            it->second->add(0, mn, r, C_DDOT(nums, &(C_->pointer(0)[abs_s][j]), nso, const_cast<double*>(&(buffer[index])), 1));
                        }
                    }
                }
            }
        } // End of RS

        for (int i=0, ij=0; i < ndocc_; i++) {
            // We need some screening here
            for (int j=0; j <= i; j++, ij++) {
                // We need some screening here
                IJ_I_J_map_.insert(std::pair<int, std::pair<int,int> >(ij, std::pair<int,int>(i,j)));
                eri_2_mn_.insert(std::pair<int, std::vector<double> >(ij,  std::vector<double>(numm*numn) ));
            }
        }

        std::map<int, std::pair<int,int> >::iterator ij_iter;
        for (ij_iter = IJ_I_J_map_.begin(); ij_iter != IJ_I_J_map_.end(); ij_iter++) {
            int ij = ij_iter->first;
            int i = IJ_I_J_map_[ij].first;
            int j = IJ_I_J_map_[ij].second;
            for (int mn=0; mn < mn_size; mn++) {
                eri_2_mn_[ij][mn] += C_DDOT(nso, &(C_->pointer(0)[0][i]), nso, &(eri_1[j]->pointer(0)[mn][0]), 1);
            }
        }
        eri_1.clear();

    } // End of exists

    return madness::None;
}


//std::string distributed_container::to_string(const int &val)
//{
//    std::stringstream strm;
//    strm <<  val;
//    return strm.str();
//}

//int distributed_container::INDEX(const int &i, const int &j)
//{
//    return (i>j) ? (pair_ioff_[i] + j) : (pair_ioff_[j] + i);
//}


//void distributed_container::init(const std::vector<int> &index,
//                                 const std::vector<int> &dist,
//                                 boost::shared_ptr<Wavefunction> wfn,
//                                 boost::shared_ptr<BasisSet> basisset,
//                                 boost::shared_ptr<Molecule> molecule,
//                                 std::vector< boost::shared_ptr<TwoElectronInt> > &eri,
//                                 const double bound, const bool nosymm)
//{

//    wfn_ = wfn;
//    basisset_ = basisset;
//    molecule_ = molecule;
//    eri_ = eri;

//    bound_ = bound;

//    distribution_.clear();
//    dist_index_size_.clear();
//    index_size_.clear();
//    mat_.clear();

//    indices_ = index.size();

//    if (indices_ != 4) throw PSIEXCEPTION("distributed_container currently only supports 4 index objects.\n");

//    for (int i=0; i < dist.size(); i++) {
//        distribution_.push_back(dist[i]);
//    }

////    for (int i=0; i < dist.size(); i++)
////        std::cout << "distribution_[" << i << "] = " << distribution_[i] << std::endl;


//    for (int i=0; i < index.size(); i++) {
//        bool distribute = 0;
//        for (int j=0; j < distribution_.size(); j++) {
//           if (i == distribution_[j]-1) distribute = 1;
//        }

//        if (distribute)
//            dist_index_size_.push_back(index[i]);
//        else
//            index_size_.push_back(index[i]);
//    }

////    for (int i=0; i < dist_index_size_.size(); i++)
////        std::cout << "dist_index_size_[" << i << "] = " << dist_index_size_[i] << std::endl;
////    for (int i=0; i < index_size_.size(); i++)
////        std::cout << "index_size_[" << i << "] = " << index_size_[i] << std::endl;


//    if (distribution_[0] > indices_) throw PSIEXCEPTION("Your first distribution index does not match your"
//                                                        " distributed container template parameter.\n");
//    else if (distribution_[1] > indices_) throw PSIEXCEPTION("Your second distribution index does not match your"
//                                                             " distributed container template parameter.\n");
//    if (index_size_[0] != index_size_[1]) throw PSIEXCEPTION("dist_containers currently only supports square matrices.\n");

//    setup_distribution_maps(nosymm);

//    row_size_.resize(pairs_per_proc_);
//    col_size_.resize(pairs_per_proc_);

//}

//void distributed_container::setup_distribution_maps(const bool nosymm)
//{
//    int dist1 = distribution_[0]-1;
//    int dist2 = distribution_[1]-1;

//    if (nosymm)
//        total_pairs_ = dist_index_size_[dist1] * dist_index_size_[dist1];
//    else
//        total_pairs_ = (dist_index_size_[dist1] * dist_index_size_[dist1] - dist_index_size_[dist1] ) / 2 + dist_index_size_[dist1];

////    std::cout << "total_pairs = " << total_pairs_ << std::endl;


//    pairs_per_proc_ = 0;
//    owner_.clear();
//    for (int i=0, count=0; i < total_pairs_; i++) {
//        owner_.push_back(i%nproc_);
//        if (me_ == owner_[i]) {
//            global_local_map_.insert(std::pair<int,int>(i,count));
//            pairs_per_proc_++;
//        }
//        if (i%nproc_ == nproc_-1) count++;
//    }

////    for (int i=0; i < total_pairs_; i++) {
////        if (me_ == 0)
////            std::cout << "owner[" << i << "] = " << owner_[i] << std::endl;
////    }
////    Communicator::world->sync();
////    if (me_ == 0) {
////        for (int i=0; i < total_pairs_; i++) {
////            if (me_ == owner_[i])
////            std::cout << "proc " << me_ << ": global_local_map_[" << i << "] = " << global_local_map_[i] << std::endl;
////        }
////    }
////    Communicator::world->sync();
////    if (me_ == 1) {
////        for (int i=0; i < total_pairs_; i++) {
////            if (me_ == owner_[i])
////            std::cout << "proc " << me_ << ": global_local_map_[" << i << "] = " << global_local_map_[i] << std::endl;
////        }
////    }

//    pair_ioff_.clear();
//    pair_ioff_.push_back(0);
//    for(int i=1; i < index_size_[0]*index_size_[1]; i++)
//        pair_ioff_.push_back(pair_ioff_[i-1]+i);

//    original_offset_ = 0;
//}

//void distributed_container::allocate_shell_block(const int &i, const int &j,
//                                                 const int &k, const int &l)
//{


//    int ik = INDEX(i,k);
//    int jl = INDEX(j,l);

//    if (me_ == owner_[ik]) {

//        int size_i = basisset_->shell(i)->nfunction();
//        int size_j = basisset_->shell(j)->nfunction();
//        int size_k = basisset_->shell(k)->nfunction();
//        int size_l = basisset_->shell(l)->nfunction();

//        mat_offset_.insert(std::pair< std::pair<int,int>, int >
//                           (std::pair<int,int>(ik,jl),original_offset_));
//        original_offset_ += size_i*size_k*size_j*size_l;

//        mat_.resize(mat_.size() + size_i*size_k*size_j*size_l);
//        mat_exist_.insert(std::pair<int,int>(ik,true));

////        mat_.insert(std::pair< int, Matrix > (ik, Matrix(size_i*size_k, size_j*size_l) ));

////        mat_[ik].print();
//    }
//}

//void distributed_container::zero_mat() {
//    ::memset(static_cast<void*> (&mat_[0]), 0, mat_.size() * sizeof(double));
//}

//void distributed_container::compute_integrals(const int &i, const int &k)
//{
//    zero_mat();
//    int ik = INDEX(i,k);


//    Matrix test;
//    if (me_ == owner_[ik]) {
//        /// may run into a proble with mat_exist if we haven't initialized it in allocat due to screening
//        if (mat_exist_[ik]) {
//            task(me_, &distributed_container::compute_shell, i, k);
//        }
//    }

//}

//madness::Void distributed_container::compute_shell(const int &I, const int &K)
//{

//    int IK = INDEX(I,K);

//    int thread = Communicator::world->thread_id( pthread_self() );
//    int numi = basisset_->shell(I)->nfunction();
//    int numk = basisset_->shell(K)->nfunction();

//    const double *buffer = eri_[thread]->buffer();


//    for (int J=0; J < basisset_->nshell(); J++) {
//        int numj = basisset_->shell(J)->nfunction();
//        for (int L=0; L < basisset_->nshell(); L++) {
//            int numl = basisset_->shell(L)->nfunction();
//            int JL = INDEX(J,L);
//            int offset = mat_offset_[std::pair<int,int>(IK,JL)];

//            eri_[thread]->compute_shell(I, J, K, L);

//            for (int i=0, index=0, mat_index=offset; i < numi; i++) {
//                for (int j=0; j < numj; j++) {
//                    for (int k=0; k < numk; k++) {
//                        for (int l=0; l < numl; l++, index++, mat_index++) {

//                            mat_[mat_index] = buffer[index];
////                            mat_[IK].set(0,ik,jl,buffer[index]);

//                        } // end of l
//                    } // end of k
//                } // end of j
//            } // end of i

//        } // end of L
//    } // end of J


//    return madness::None;
//}


#endif // have_madness

} // End of psi namespace
