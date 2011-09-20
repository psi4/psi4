#ifndef _psi_src_lib_libmints_helper_h
#define _psi_src_lib_libmints_helper_h

#include <libiwl/iwl.hpp>
#include <vector>
#include "wavefunction.h"
#include "multipolesymmetry.h"

namespace psi {

class Options;
class CdSalcList;

/**
* IWLWriter functor for use with SO TEIs
**/
class IWLWriter {
    IWL& writeto_;
    size_t count_;
    int& current_buffer_count_;

    Label *plabel_;
    Value *pvalue_;
public:

    IWLWriter(IWL& writeto) : writeto_(writeto), count_(0),
        current_buffer_count_(writeto_.index())
    {
        plabel_ = writeto_.labels();
        pvalue_ = writeto_.values();
    }

    void operator()(int i, int j, int k, int l, int , int , int , int , int , int , int , int , double value)
    {
        int current_label_position = 4*current_buffer_count_;

        // Save the labels
        plabel_[current_label_position++] = i;
        plabel_[current_label_position++] = j;
        plabel_[current_label_position++] = k;
        plabel_[current_label_position]   = l;

        // Save the value
        pvalue_[current_buffer_count_++] = value;

        // Increment overall counter
        count_++;

        // If our IWL buffer is full dump to disk.
        if (current_buffer_count_ == writeto_.ints_per_buffer()) {
            writeto_.last_buffer() = 0;
            writeto_.buffer_count() = current_buffer_count_;
            writeto_.put();
            current_buffer_count_ = 0;
        }
    }

    size_t count() const { return count_; }
};

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
    boost::shared_ptr<Matrix> mo_eri_helper(boost::shared_ptr<Matrix> Iso, boost::shared_ptr<Matrix> Co, boost::shared_ptr<Matrix> Cv);
    // In-core O(N^5) transqt
    boost::shared_ptr<Matrix> mo_eri_helper(boost::shared_ptr<Matrix> Iso, boost::shared_ptr<Matrix> C1, boost::shared_ptr<Matrix> C2,
                                                                           boost::shared_ptr<Matrix> C3, boost::shared_ptr<Matrix> C4);

public:

    void init_helper(boost::shared_ptr<Wavefunction> wavefunction = boost::shared_ptr<Wavefunction>());
    /// Constructor, just lines references up
    MintsHelper(Options&, int print = 1);
    /// Constructor, uses globals
    MintsHelper();
    /// Constructor, using wavefunction
    MintsHelper(boost::shared_ptr<Wavefunction> wavefunction);
    /// Destructor, does nothing
    ~MintsHelper();

    OperatorSymmetry operator_symmetry(int order) {
        return OperatorSymmetry(order, molecule_, integral_, factory_);
    }

    /// Basis set being used.
    boost::shared_ptr<BasisSet> basisset() const;
    /// SO basis set being used.
    boost::shared_ptr<SOBasisSet> sobasisset() const;
    /// Matrix factory being used
    boost::shared_ptr<MatrixFactory> factory() const;

    /// Molecular integrals (just like cints used to do)
    void integrals();
    /// Standard one electron integrals (just like oeints used to do)
    void one_electron_integrals();
    /// Derivative integrals (not implemented)
    void integral_gradients();
    /// Hessian integrals (not implemented)
    void integral_hessians();

    /// AO ERI Integrals (Full matrix, not recommended for large systems)
    boost::shared_ptr<Matrix> ao_eri();
    /// Symmetric MO ERI Integrals, (ov|ov) type  (Full matrix, N^5, not recommended for large systems)
    /// Pass C_ C_ for (aa|aa) type, Cocc_, Cocc_ for (oo|oo) type, or Cvir_, Cvir_ for (vv|vv) type
    boost::shared_ptr<Matrix> mo_eri(boost::shared_ptr<Matrix> Cocc, boost::shared_ptr<Matrix> Cvir);
    /// Non Symmetric MO ERI Omega Integrals, (12|34) type  (Full matrix, N^5, not recommended for large systems)
    boost::shared_ptr<Matrix> mo_eri(boost::shared_ptr<Matrix> C1, boost::shared_ptr<Matrix> C2,
                                     boost::shared_ptr<Matrix> C3, boost::shared_ptr<Matrix> C4);
    /// AO ERI Omega Integrals (Full matrix, not recommended for large systems)
    boost::shared_ptr<Matrix> ao_erf_eri(double omega);
    /// Symmetric MO ERI Omega Integrals, (ov|ov) type  (Full matrix, N^5, not recommended for large systems)
    /// Pass C_ C_ for (aa|aa) type, Cocc_, Cocc_ for (oo|oo) type, or Cvir_, Cvir_ for (vv|vv) type
    boost::shared_ptr<Matrix> mo_erf_eri(double omega, boost::shared_ptr<Matrix> Cocc, boost::shared_ptr<Matrix> Cvir);
    /// Non Symmetric MO ERI Omega Integrals, (12|34) type  (Full matrix, N^5, not recommended for large systems)
    boost::shared_ptr<Matrix> mo_erf_eri(double omega, boost::shared_ptr<Matrix> C1, boost::shared_ptr<Matrix> C2,
                                                       boost::shared_ptr<Matrix> C3, boost::shared_ptr<Matrix> C4);

    /// AO Overlap Integrals
    boost::shared_ptr<Matrix> ao_overlap();
    /// AO Kinetic Integrals
    boost::shared_ptr<Matrix> ao_kinetic();
    /// AO Potential Integrals
    boost::shared_ptr<Matrix> ao_potential();
    /// Vector AO Angular Momentum Integrals
    std::vector<boost::shared_ptr<Matrix> > ao_angular_momentum();
    /// Vector AO Nabla Integrals
    std::vector<boost::shared_ptr<Matrix> > ao_nabla();
    /// AO Overlap Integrals
    boost::shared_ptr<Matrix> so_overlap();
    /// AO Kinetic Integrals
    boost::shared_ptr<Matrix> so_kinetic();
    /// AO Potential Integrals
    boost::shared_ptr<Matrix> so_potential();
    /// Vector SO Dipole Integrals
    std::vector<boost::shared_ptr<Matrix> > so_dipole();
    /// Vector SO Nabla Integrals
    std::vector<boost::shared_ptr<Matrix> > so_nabla();
    /// Vector SO Angular Momentum Integrals
    std::vector<boost::shared_ptr<Matrix> > so_angular_momentum();
    /// Vector SO Quadrupole Integrals
    std::vector<boost::shared_ptr<Matrix> > so_quadrupole();
    /// Vector SO Traceless Quadrupole Integrals
    std::vector<boost::shared_ptr<Matrix> > so_traceless_quadrupole();

    /// Returns a CdSalcList object
    boost::shared_ptr<CdSalcList> cdsalcs(int needed_irreps=0xF,
                                          bool project_out_translations=true,
                                          bool project_out_rotations=true);

    /// Play function
    void play();
};

}

#endif
