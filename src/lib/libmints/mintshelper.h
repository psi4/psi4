#ifndef _psi_src_lib_libmints_helper_h
#define _psi_src_lib_libmints_helper_h

#include <vector>
#include "wavefunction.h"
#include "multipolesymmetry.h"

namespace psi {

class Options;
class CdSalcList;

/**
* The MintsHelper object, places molecular integrals
* (and later derivative integrals) on disk
**/
class MintsHelper {

private:
    /// The Options reference for basis sets and things
    Options& options_;
    boost::shared_ptr<PSIO> psio_;
    boost::shared_ptr<MatrixFactory> factory_;
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<IntegralFactory> integral_;
    boost::shared_ptr<BasisSet> basisset_;
    boost::shared_ptr<SOBasisSet> sobasis_;
    int print_;

    // In-core O(N^5) transqt
    SharedMatrix mo_eri_helper(SharedMatrix Iso, SharedMatrix Co, SharedMatrix Cv);
    // In-core O(N^5) transqt
    SharedMatrix mo_eri_helper(SharedMatrix Iso, SharedMatrix C1, SharedMatrix C2,
                                                                           SharedMatrix C3, SharedMatrix C4);

public:

    void init_helper(boost::shared_ptr<Wavefunction> wavefunction = boost::shared_ptr<Wavefunction>());
    void init_helper_2(boost::shared_ptr<BasisSet> basis);

    /// Constructor, just lines references up
    MintsHelper(Options&, int print = 1);
    /// Constructor, uses a basis set
    MintsHelper(boost::shared_ptr<BasisSet> basis);
    /// Constructor, uses globals
    MintsHelper();
    /// Constructor, using wavefunction
    MintsHelper(boost::shared_ptr<Wavefunction> wavefunction);
    /// Destructor, does nothing
    ~MintsHelper();

    OperatorSymmetry operator_symmetry(int order) {
        return OperatorSymmetry(order, molecule_, integral_, factory_);
    }

    /// Returns the number of basis functions
    int nbf() const;

    /// Returns petite list that is capable of transforming basis functions (nbf) to SO's.
    boost::shared_ptr<PetiteList> petite_list() const;

    enum {
        kFromCartesianAO = true,
        kFromBF = false
    };
    /** Returns petite list that is capable of transforming AO basis functions (nbf) to SO's.
     *  \param include_pure_transform Is either kFromCartesianAO or kFromBF.
     */
    boost::shared_ptr<PetiteList> petite_list(bool include_pure_transform) const;
    /// Basis set being used.
    boost::shared_ptr<BasisSet> basisset() const;
    /// SO basis set being used.
    boost::shared_ptr<SOBasisSet> sobasisset() const;
    /// Matrix factory being used
    boost::shared_ptr<MatrixFactory> factory() const;
    /// Integral factory being used
    boost::shared_ptr<IntegralFactory> integral() const;

    /// Molecular integrals (just like cints used to do)
    void integrals();
    /// Standard one electron integrals (just like oeints used to do)
    void one_electron_integrals();
    /// Derivative integrals (not implemented)
    void integral_gradients();
    /// Hessian integrals (not implemented)
    void integral_hessians();

    /// AO ERI Integrals (Full matrix, not recommended for large systems)
    SharedMatrix ao_eri();
    /// Symmetric MO ERI Integrals, (ov|ov) type  (Full matrix, N^5, not recommended for large systems)
    /// Pass C_ C_ for (aa|aa) type, Cocc_, Cocc_ for (oo|oo) type, or Cvir_, Cvir_ for (vv|vv) type
    SharedMatrix mo_eri(SharedMatrix Cocc, SharedMatrix Cvir);
    /// Non Symmetric MO ERI Omega Integrals, (12|34) type  (Full matrix, N^5, not recommended for large systems)
    SharedMatrix mo_eri(SharedMatrix C1, SharedMatrix C2,
                                     SharedMatrix C3, SharedMatrix C4);
    /// AO ERI Omega Integrals (Full matrix, not recommended for large systems)
    SharedMatrix ao_erf_eri(double omega);
    /// Symmetric MO ERI Omega Integrals, (ov|ov) type  (Full matrix, N^5, not recommended for large systems)
    /// Pass C_ C_ for (aa|aa) type, Cocc_, Cocc_ for (oo|oo) type, or Cvir_, Cvir_ for (vv|vv) type
    SharedMatrix mo_erf_eri(double omega, SharedMatrix Cocc, SharedMatrix Cvir);
    /// Non Symmetric MO ERI Omega Integrals, (12|34) type  (Full matrix, N^5, not recommended for large systems)
    SharedMatrix mo_erf_eri(double omega, SharedMatrix C1, SharedMatrix C2,
                                                       SharedMatrix C3, SharedMatrix C4);

    /// AO Overlap Integrals
    SharedMatrix ao_overlap();
    /// AO Kinetic Integrals
    SharedMatrix ao_kinetic();
    /// AO Potential Integrals
    SharedMatrix ao_potential();
    /// Vector AO Dipole Integrals
    std::vector<SharedMatrix> ao_dipole();
    /// Vector AO Angular Momentum Integrals
    std::vector<SharedMatrix > ao_angular_momentum();
    /// Vector AO Nabla Integrals
    std::vector<SharedMatrix > ao_nabla();
    /// AO Alchemical Potential
    SharedMatrix ao_alchemical_potential();
    /// SO Alchemical Potential
    SharedMatrix so_alchemical_potential();
    /// AO Overlap Integrals
    SharedMatrix so_overlap();
    /// AO Kinetic Integrals
    SharedMatrix so_kinetic();
    /// AO Potential Integrals
    SharedMatrix so_potential();
    /// Vector SO Dipole Integrals
    std::vector<SharedMatrix > so_dipole();
    /// Vector SO Nabla Integrals
    std::vector<SharedMatrix > so_nabla();
    /// Vector SO Angular Momentum Integrals
    std::vector<SharedMatrix > so_angular_momentum();
    /// Vector SO Quadrupole Integrals
    std::vector<SharedMatrix > so_quadrupole();
    /// Vector SO Traceless Quadrupole Integrals
    std::vector<SharedMatrix > so_traceless_quadrupole();

    /// Returns a CdSalcList object
    boost::shared_ptr<CdSalcList> cdsalcs(int needed_irreps=0xF,
                                          bool project_out_translations=true,
                                          bool project_out_rotations=true);

    /// Play function
    void play();
};

}

#endif
