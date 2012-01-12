#ifndef SIEVE_H 
#define SIEVE_H

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class BasisSet;
class TwoBodyAOInt;

/**
 * ERISieve
 *
 * Class to perform ERI sieving
 *
 * At the moment, this class uses Cauchy-Schwarz sieving,
 * and is designed to be used in concert with density-based
 * sieving. In the future, it is hoped that MBIE sieving may 
 * be added to this class. 
 *
 * Use of this class is as follows:
 *
 *     // Initialize the sieve object
 *     boost::shared_ptr<ERISieve> sieve(basisset, sieve_cutoff);
 *
 *     // Reset the sieve cutoff (you can do this wherever)
 *     sieve->set_sieve(new_cutoff);
 *
 *     ...
 *      
 *     // Investigate stuff. This example is for shells, all methods are
 *     // also provided for functions  
 *
 *     // Compute a shell quartet (MN|RS), if (MN|RS) >= sieve_cutoff
 *     if (sieve->shell_significant(M,N,R,S)) eri->compute(M,N,R,S); 
 *
 *     // Compute a shell quartet (MN|RS), if (MN|RS) * D_RS >= sieve_cutoff 
 *     // Squares are used to avoid sqrt()
 *     if (sieve->shell_ceiling2(M,N,R,S) * D_RS * D_RS >= sieve_cutoff * sieve_cutoff)
 *         eri->compute(M,N,R,S); 
 *
 *     // Index the significant MN shell pairs (triangular M,N)
 *     const std::vector<std::pair<int,int> >& MN = sieve->shell_pairs();
 *     for (long int index = 0L; index < MN.size(); ++index) {
 *         int M = MN[index].first;
 *         int N = MN[index].second;
 *     } 
 *
 *     // Check if a triangular index MNindex (M * (M + 1) / 2) + N exists,
 *     // and if so, where it starts in reduced triangular MN
 *     int MNindex = (M * (M + 1) >> 1) + N;
 *     const std::vector<long int> >& MN_reverse = sieve->shell_pairs_reverse();
 *     int MNreduced = MN_reverse[MNindex];
 *     if (MNreduced < 0) {
 *         // The shell pair is not signficant
 *     } else {
 *         // The shell pair is the MNreduced significant shell pair
 *     }  
 *       
 *
 */
class ERISieve {

protected:

    /// Debug flag (defaults to 0)
    int debug_;

    /// Basis set reference
    boost::shared_ptr<BasisSet> primary_;

    /// Number of basis functions
    int nbf_;
    /// Number of shells
    int nshell_;    

    /// Cutoff values
    double sieve_;
    /// Maximum |(mn|ls)|
    double max_;
    /// sieve_ / max_
    double sieve_over_max_;
    /// sieve_ * sieve_
    double sieve2_;
    /// sieve_ * sieve_ / max_
    double sieve2_over_max_;

    /// |(mn|mn)| values (nbf * nbf)
    double* function_pair_values_;
    /// max |(MN|MN)| values (nshell * nshell)
    double* shell_pair_values_;

    /// Significant unique bra- function pairs, in reduced triangular indexing
    std::vector<std::pair<int,int> > function_pairs_;
    /// Significant unique bra- shell pairs, in reduced triangular indexing
    std::vector<std::pair<int,int> > shell_pairs_;
    /// Unique bra- function pair indexing, accessed in triangular order, or -1 for non-significant pair 
    std::vector<long int> function_pairs_reverse_;
    /// Unique bra- shell pair indexing, accessed in triangular order, or -1 for non-significant pair 
    std::vector<long int> shell_pairs_reverse_;
    /// Significant function pairs, indexes by function
    std::vector<std::vector<int> > shell_to_shell_;
    /// Significant shell pairs, indexes by shell
    std::vector<std::vector<int> > function_to_function_;
     
    /// Set initial indexing
    void common_init();
    /// Compute sieve integrals (only done once)
    void integrals();

public:

    /// Constructor, basis set and first sieve cutoff
    ERISieve(boost::shared_ptr<BasisSet> primary, double sieve = 0.0);
    /// Destructor, frees memory
    virtual ~ERISieve();

    /// Set sieve value and redo indexing
    void set_sieve(double sieve);
    /// Get sieve cutoff value 
    double sieve() const { return sieve_; }
    /// Global maximum |(mn|rs)|
    double max() const { return max_; }
    
    // => Significance Checks <= //

    /// Square of ceiling of shell quartet (MN|RS)
    inline double shell_ceiling2(int M, int N, int R, int S) { 
        return shell_pair_values_[N * (unsigned long int) nshell_ + M] * 
               shell_pair_values_[R * (unsigned long int) nshell_ + S]; } 

    /// Square of ceiling of integral (mn|rs)
    inline double function_ceiling2(int m, int n, int r, int s) { 
        return function_pair_values_[m * (unsigned long int) nbf_ + n] * 
               function_pair_values_[r * (unsigned long int) nbf_ + s]; } 

    /// Is the shell quartet (MN|RS) significant according to sieve? (no restriction on MNRS order)
    inline bool shell_significant(int M, int N, int R, int S) { 
        return shell_pair_values_[N * (unsigned long int) nshell_ + M] * 
               shell_pair_values_[R * (unsigned long int) nshell_ + S] >= sieve2_; } 

    /// Is the integral (mn|rs) significant according to sieve? (no restriction on mnrs order)
    inline bool function_significant(int m, int n, int r, int s) { 
        return function_pair_values_[m * (unsigned long int) nbf_ + n] * 
               function_pair_values_[r * (unsigned long int) nbf_ + s] >= sieve2_; } 

    // => Indexing [these change after a call to sieve()] <= //

    /// Significant unique bra- function pairs, in reduced triangular indexing
    const std::vector<std::pair<int,int> >& function_pairs() const { return function_pairs_; }
    /// Significant unique bra- shell pairs, in reduced triangular indexing
    const std::vector<std::pair<int,int> >& shell_pairs() const { return shell_pairs_; }
    /// Unique bra- function pair indexing, accessed in triangular order, or -1 for non-significant pair 
    const std::vector<long int> function_pairs_reverse() const { return function_pairs_reverse_; }
    /// Unique bra- shell pair indexing, accessed in triangular order, or -1 for non-significant pair 
    const std::vector<long int> shell_pairs_reverse() const { return shell_pairs_reverse_; }
    /// Significant function pairs, indexes by function
    const std::vector<std::vector<int> >& function_to_function() const { return function_to_function_; }
    /// Significant shell pairs, indexes by shell
    const std::vector<std::vector<int> >& shell_to_shell() const { return shell_to_shell_; }

    /// Set debug flag (defaults to 0)
    void set_debug(int debug) { debug_ = debug; }

};


}
#endif
