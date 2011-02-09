#ifndef _psi_src_lib_libmints_helper_h
#define _psi_src_lib_libmints_helper_h

#include <libiwl/iwl.hpp>

namespace psi {

class Options;
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

public:
    /// Constructor, just lines references up
    MintsHelper(Options&);
    /// Constructor, uses globals
    MintsHelper();
    /// Destructor, does nothing
    ~MintsHelper();
    /// Molecular integrals (just like cints used to do)
    void integrals();
    /// One electron integrals (just like oeints used to do)
    void one_electron_integrals();
    /// Derivative integrals (not implemented)
    void integral_gradients();
    /// Hessian integrals (not implemented)
    void integral_hessians();
    /// AO Overlap Integrals
    shared_ptr<Matrix> ao_overlap();
    /// AO Kinetic Integrals
    shared_ptr<Matrix> ao_kinetic();
    /// AO Potential Integrals
    shared_ptr<Matrix> ao_potential();
    /// AO ERI Integrals (Full matrix, not recommended for large systems)
    shared_ptr<Matrix> ao_eri();
    /// AO ERI Omega Integrals (Full matrix, not recommended for large systems)
    shared_ptr<Matrix> ao_erf_eri(double omega, double alpha, double beta);

    /// AO Overlap Integrals
    shared_ptr<Matrix> so_overlap();
    /// AO Kinetic Integrals
    shared_ptr<Matrix> so_kinetic();
    /// AO Potential Integrals
    shared_ptr<Matrix> so_potential();
};

}

#endif
