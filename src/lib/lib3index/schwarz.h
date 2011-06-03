#ifndef three_index_schwarz_H
#define three_index_schwarz_H

#include <psi4-dec.h>
#include <psiconfig.h>

namespace psi {

class BasisSet;

class SchwarzSieve {

protected:

    // The schwarz cutoff
    double schwarz_;

    // Basis set for this schwarz
    boost::shared_ptr<BasisSet> basis_;

    // Is the sieve initialized
    bool initialized_;

    // number of significant shell pairs
    unsigned long int nshell_pairs_;
    // number of significant function pairs
    unsigned long int nfun_pairs_;

    double max_global_val_;
    int* schwarz_shells_;
    int* schwarz_funs_;

    long int* schwarz_shells_reverse_;
    long int* schwarz_funs_reverse_;

    double* schwarz_shell_vals_;
    double* schwarz_fun_vals_;

    void form_schwarz_ints();

public:
    SchwarzSieve(boost::shared_ptr<BasisSet>, double cutoff);
    virtual ~SchwarzSieve();

    void form_schwarz_sieve(double cutoff);
    // Sizes of the significant bra/ket pairs
    unsigned long int get_nshell_pairs() const { return nshell_pairs_; }
    unsigned long int get_nfun_pairs() const { return nfun_pairs_; }
    // I_global = arr[2*I_local], J_global = arr[2*I_local + 1]
    // These are only defined up to nshell_pairs_ and nfun_pairs_, respectively
    int* get_schwarz_shells() const { return schwarz_shells_; }
    int* get_schwarz_funs() const { return schwarz_funs_; }
    // Canonical compound indexing, -1 if not present
    long int* get_schwarz_shells_reverse() const { return schwarz_shells_reverse_; }
    long int* get_schwarz_funs_reverse() const { return schwarz_funs_reverse_; }

};


}
#endif
