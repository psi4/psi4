#ifndef MAD_MP2_PROC_H
#define MAD_MP2_PROC_H

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <libparallel/parallel.h>
#include <map>

#if HAVE_MADNESS

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class Denominator;
class Matrix; 
class Vector; 
class IntVector; 
class BasisSet; 
class Wavefunction;
class Options;
class PSIO;
class Chkpt;

namespace  mad_mp2 {

class MAD_MP2 : public Wavefunction , public madness::WorldObject<MAD_MP2>
{

protected:

    // => Parallel Variables <= //

#ifdef HAVE_MADNESS
     SharedMadWorld madworld_;
     SharedMutex mutex_;
#endif
    
    // Number of processors
    int nproc_;
    // Number of threads
    int mad_nthread_;
    // My rank
    int rank_; 
    // Communicator Type
    std::string comm_; 
    // Total number of ia pairs   
    ULI nia_;
    // Local number of ia pairs
    ULI nia_local_;
    // ia pair assignments
    std::vector<int> ia_owner_;
    // ablock owner for given i
    std::vector<std::vector<int> > ablock_owner_;
    // global a ablock starts for given i
    std::vector<std::vector<int> > ablock_start_;
    // global a ablock sizes for given i
    std::vector<std::vector<int> > ablock_size_;
    // Local ia pair assignments
    std::vector<ULI> ia_local_to_global_;
    // Global ia pair assignments 
    std::map<ULI, int> ia_global_to_local_;
    // Global i assignments
    std::map<int,int> i_global_to_local_;
    // Local number of i
    int naocc_local_;
    // Local number of a
    int navir_local_;
    // Local i values
    std::vector<int> aocc_local_;
    // Local a values
    std::vector<int> avir_local_;

    /// Energies table
    std::map<std::string, double> energies_;

    /// Total MP2J terms
    double E_MP2J_;
    /// Total MP2K terms
    double E_MP2K_;

    /// Same-spin scale
    double scale_ss_;
    /// Opposite-spin scale
    double scale_os_;
   
    /// Print flag
    int print_;
    /// Debug flag
    int debug_;   
    /// Number of OMP threads
    int omp_nthread_;

    /// Size of auxiliary basis set
    int naux_;
    /// Auxiliary basis set
    boost::shared_ptr<BasisSet> auxiliary_;
    /// Zero basis set
    boost::shared_ptr<BasisSet> zero_;   
    /// Auxiliary basis automagical?
    bool auxiliary_automatic_;

    /// Reference wavefunction pointer
    boost::shared_ptr<Wavefunction> reference_;
    /// SCF energy
    double Eref_;

    // => Orbital sizing/data <= //

    /// Total number of frozen occupied orbitals
    int nfocc_;
    /// Total number of active occupied orbitals
    int naocc_;   
    /// Total number of active virtual orbitals
    int navir_;   
    /// Total number of frozen virtual orbitals
    int nfvir_;

    /// Number of active occupieds per irrep
    int naoccpi_[8];
    /// Number of active virtuals per irrep
    int navirpi_[8];
    /// Cumsum of naoccpi_
    int offset_aocc_[8];   
    /// Cumsum of navirpi_
    int offset_avir_[8];   
 
    /// AO2USO transform matrix
    SharedMatrix AO2USO_; 
    /// AO2USO transform matrix (auxiliary)
    SharedMatrix AO2USO_aux_; 
    /// SOs per irrep in the auxiliary basis
    int nauxpi_[8];
   
    /// Max active occupieds per irrep
    int max_naoccpi_;
    /// Max active occupieds per irrep
    int max_navirpi_;
    /// Max active occupieds per irrep
    int max_nauxpi_;
 
    /// C1 copy of active occupied orbitals
    SharedMatrix Caocc_;
    /// C1 copy of active virtual orbitals
    SharedMatrix Cavir_;
   
    /// C1 copy of active occupied evals
    boost::shared_ptr<Vector> eps_aocc_;
    /// C1 copy of active virtual evals
    boost::shared_ptr<Vector> eps_avir_;

    /// C1 copy of active occupied irreps
    boost::shared_ptr<IntVector> irrep_aocc_;
    /// C1 copy of active virtual irreps
    boost::shared_ptr<IntVector> irrep_avir_;

    // => Key Tensors <= //

    /// J^-1/2
    SharedMatrix Jm12_;
    /// (A|ia) 
    SharedMatrix Aia_; 
    /// \tau_ia^Q
    boost::shared_ptr<Denominator> denom_;   
 
    // => Compute routines <= //

    /// Handle sizing, orbital evals, and C matrix
    virtual void common_init();
    /// Handle parallel initialization/sizing
    virtual void parallel_init();
    /// Print the header
    virtual void print_header();
    /// Build the MP2J Cholesky or Laplace denominator
    virtual void denominator();
    /// Detemine memory, and striping information
    virtual void check_memory();
    /// Build J
    virtual void J();
    /// Build J^-1/2
    virtual void Jm12();
    /// Build Aia, Ami, Amn on the go
    virtual void Aia();
    /// Build I, find contributions to energy
    virtual void I();
    /// Build I_MP2J, find contributions to energy
    virtual void IJ();
    /// Print the energies
    void print_energy();

    madness::Future<std::vector<double> > fetch_Qia_block(const int& i, const int& ablock);
    madness::Void unpack_Qia_block(const std::vector<double>& block, SharedMatrix Q,
                                   const int& astart, const int& asize,
                                   const int &i);
    madness::Future<SharedMatrix> build_Qa(const int &i);
    madness::Future<SharedMatrix> build_Qb(const int &j);
    madness::Future<SharedMatrix> build_I(SharedMatrix Qa, SharedMatrix Qb);
    madness::Future<double> energy_j(const SharedMatrix I,const int &i, const int &j);
    madness::Future<double> energy_k(const SharedMatrix I,const int &i, const int &j);

    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }

public:
    MAD_MP2(Options& options, boost::shared_ptr<PSIO> psio);
    virtual ~MAD_MP2();

    /// Pure virtual from Wavefunction
    virtual double compute_energy();

};

}} // Namespaces
#endif //HAVE_MADNESS

#endif

