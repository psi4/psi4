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

    public:

    
    IWLWriter(IWL& writeto) : writeto_(writeto), count_(0)
      { }

    void operator()(int i, int j, int k, int l, int , int , int , int , int , int , int , int , double value)
    {
        writeto_.write_value(i, j, k, l, value, 0, NULL, 0);
        count_++;
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
        /// Uniqueness forcer until SOShellCombinationsIterator is complete
        static int determine_unique_shell_quartets(int usii, int usjj, int uskk, int usll,
                                           int* usi_arr,
                                           int* usj_arr,
                                           int* usk_arr,
                                           int* usl_arr);
    public:
        /// Constructor, just lines references up
        MintsHelper(Options&);
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

};

}

#endif
