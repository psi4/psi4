/*! \file
    \ingroup Distributed Container
    \brief Enter brief description of file here
*/
#ifndef _psi_src_lib_DIST_CONTAINER_H
#define _psi_src_lib_DIST_CONTAINER_H

#include <string>
#include <sstream>

#include <libmints/mints.h>
#include <libparallel/parallel.h>
#include <libqt/qt.h>
#include "psi4-dec.h"

#ifdef HAVE_MADNESS
    #include <world/worldobj.h>
    #include <world/worlddc.h>
#endif


namespace boost {
    template<class T> class shared_ptr;
}

namespace psi{

#ifdef HAVE_MADNESS

class dist_key {

private:
    int k;

public:
    dist_key() : k(-1) {}

    dist_key(int k) : k(k) {}

    madness::hashT hash() const {
        return k;
    }

    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & k;
    }

    bool operator==(const dist_key& b) const {
        return k==b.k;
    }
};

class dist_container {

    typedef std::map<int, std::vector<double> > block;
private:
    block eri_2_mn_;
    block eri_2_ij_;
    int MN_;
    int me_;
    int nproc_;
    int nthread_;
    int ndocc_;
    std::string comm_;
    double thresh_;
    bool exists_;

    std::pair<int,bool> RS_exist_;
    std::map<int,std::pair<int,int> > RS_map_;
    std::vector<int> IJ_I_map_;
    std::map<int, std::pair<int,int> > IJ_I_J_map_;

    SharedMatrix C_;

    // Need to serialize Wavefunction, BasisSet, Molecule, and TwoElectronInt
    boost::shared_ptr<Wavefunction> wfn_;
    boost::shared_ptr<BasisSet> basisset_;
    boost::shared_ptr<Molecule> molecule_;
    std::vector< boost::shared_ptr<TwoElectronInt> > eri_;


public:
    dist_container();

    void init_wfn(boost::shared_ptr<Wavefunction> wfn,
                  boost::shared_ptr<BasisSet> basisset,
                  boost::shared_ptr<Molecule> molecule,
                  std::vector<boost::shared_ptr<TwoElectronInt> > &eri);

    void init_MN(const int &MN, const int &ndocc,
                 SharedMatrix C, const double bound);

    block get_mn() const { return eri_2_mn_; }

    int pair() const { return MN_; }
    int owner() const { return me_; }

    void print_mn();

    madness::Void first_half_transformation(const int &M, const int &N);

    madness::Void allocate_RS_block();

    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & eri_2_mn_ & MN_ & me_ & exists_;
    }
};


typedef madness::WorldContainer<psi::dist_key,psi::dist_container> distributed_container;

class key
{
  private:
    int k;

  public:
    key() : k(-1) {}

    key(int k) : k(k) {}

    madness::hashT hash() const
    { return k; }

    template <typename Archive>
    inline void serialize(const Archive& ar)
    { ar & k; }

    bool operator==(const key& b) const
    { return k == b.k; }

};

class value
{
  private:
    double val_;

  public:
    value() : val_(-1) { }

    value(double k) : val_(k) { }

    double get_value() const
    {
        return val_;
    }

    madness::Void print(const int &me)
    {
        std::cout << "owner = " << me << ": value = " << val_ << std::endl;
    }

    template <typename Archive>
    inline void serialize(const Archive& ar)
    {
        ar & val_;
    }
};

typedef madness::WorldContainer<key,value> dc;

//class distributed_container :
//        public madness::WorldObject<distributed_container>
//{

//private:

//    boost::shared_ptr<Wavefunction> wfn_;
//    boost::shared_ptr<BasisSet> basisset_;
//    boost::shared_ptr<Molecule> molecule_;
//    std::vector< boost::shared_ptr<TwoElectronInt> > eri_;


//    /// Each process owns a unique matrix
//    std::vector<double> mat_;
//    std::map<std::pair<int,int>, int> mat_offset_;
//    std::map<int,bool> mat_exist_;
//    int original_offset_;
////    std::map<int, Matrix > mat_;
////    std::vector<Matrix> mat_;

//    /// The cutoff value for screening
//    double bound_;

//    /// A vector containing the sizes of the distributed indices
//    std::vector<int> dist_index_size_;
//    /// A vector containing the sizes of the non-distributed indices
//    std::vector<int> index_size_;

//    /// The total number of indices
//    size_t indices_;

//    /// The 1st index to distribute over
//    std::vector<int> distribution_;

//    // ****  These are for parallel computing ***
//    /// The process ID
//    int me_;
//    /// The total number of processes
//    int nproc_;
//    /// The number of threads
//    int nthread_;
//    /// The communicator type
//    std::string comm_;
//    /// The madness world communicator
//#ifdef HAVE_MADNESS
//     SharedMadWorld madworld_;
//#endif


//    /// index lookup array
//    std::vector<int> pair_ioff_;

//    /// The total number of pairs that will be distributed
//    int total_pairs_;
//    /// The number of pairs per proc
//    int pairs_per_proc_;
//    /// This is a vector of who owns each matrix
//    std::vector<int> owner_;
//    /// global to local map
//    std::map<int,int> global_local_map_;
//    /// local to global map
//    std::vector<int> local_global_map_;
//    /// This is a vector of the size of the row for each matrix
//    std::vector<int> row_size_;
//    /// This is a vector of the size of the row for each matrix
//    std::vector<int> col_size_;

//    /// Setup up the owners
//    void setup_distribution_maps(const bool nosymm=false);
//    /// Convert an integer to a string
//    std::string to_string(const int &val);
//    /// Ioff lookup function
//    int INDEX(const int &i, const int &j);
//    /// This actually computes the all integrals for the given ik
//    madness::Void compute_shell(const int &i, const int &k);
//    /// Zero out the mat_ array
//    void zero_mat();


//public:
//    distributed_container() :
//        madness::WorldObject<distributed_container>(*Communicator::world->get_madworld())
//    {
//        me_ = Communicator::world->me();
//        nproc_ = Communicator::world->nproc();
//        nthread_ = Communicator::world->nthread();
//        comm_ = Communicator::world->communicator();

//        process_pending();
//    }

//    /// Initialize an n dimensional distributed array
//    void init(const std::vector<int> &index,
//              const std::vector<int> &dist,
//              boost::shared_ptr<Wavefunction> wfn,
//              boost::shared_ptr<BasisSet> basisset,
//              boost::shared_ptr<Molecule> molecule,
//              std::vector< boost::shared_ptr<TwoElectronInt> > &eri,
//              const double bound, const bool nosymm=false);

//    /// Allocate a matrix for the given indices
//    void allocate_shell_block(const int &i, const int &j,
//                              const int &k, const int &l);

//    /** This distributes the integrals according to who owns the
//      * given pair, and submitts a task to compute the integrals
//      */
//    void compute_integrals(const int &i, const int &k);

//    /// Returns the vector that contains the distribution indices
//    std::vector<int> distribution() {return distribution_;}


//};

#endif // have_madness

}



#endif // DIST_CONTAINER_H
