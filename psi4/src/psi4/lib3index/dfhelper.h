/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef three_index_dfhelper
#define three_index_dfhelper

#include "psi4/psi4-dec.h"
#include <psi4/libmints/typedefs.h>
#include "psi4/libpsi4util/exception.h"

#include <map>
#include <list>
#include <vector>
#include <tuple>
#include <string>

namespace psi {

class BasisSet;
class Matrix;
class TwoBodyAOInt;

// COMMON VARIABLE NAMES
// lr_symmetric: Are the two C matrices in our J/K-esque contraction equal?

class PSI_API DFHelper {
   public:
    DFHelper(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> aux);
    ~DFHelper();

    ///
    /// Specify workflow for transforming and contracting integrals
    /// @param method (STORE or DIRECT or DIRECT_iaQ) to indicate workflow
    /// STORE: contract and save AO integrals before transforming
    /// DIRECT: pre-transform integrals before metric contraction
    /// DIRECT_iaQ: special workflow when using integrals of the iaQ form
    /// (defaults to STORE)
    ///
    void set_method(std::string method) { method_ = method; }
    std::string get_method() { return method_; }

    ///
    /// Tells DFHelper how many threads to spawn in parallel regions
    /// @param nthreads specifies number of threads to use
    /// oversubcription is possible
    ///
    void set_nthreads(size_t nthreads) { nthreads_ = nthreads; }
    size_t get_nthreads() { return nthreads_; }

    ///
    /// Indicates the memory (in doubles) DFHelper gets to control
    /// @doubles specifies number of doubles given for memory
    /// defaults to 256,000,000 (2.04GB)
    ///
    void set_memory(size_t doubles) { memory_ = doubles; }
    size_t get_memory() { return memory_; }

    ///
    /// Indicates whether to use the in-core or out-of-core
    //  subalgorithm available in MemDFJK, or to let
    //  Psi4 select automatically based on allocated memory
    /// defaults to AUTO (Psi4 selects by default)
    ///
    void set_subalgo(std::string subalgo) { subalgo_ = subalgo; }

    /// Returns the number of doubles in the *screened* AO integrals
    size_t get_AO_size() { return big_skips_[nbf_]; }

    /// Returns the size of the in-core version in doubles
    size_t get_core_size() {
        AO_core(false);
        return required_core_size_;
    }

    /// Returns the amount of sparsity in the AO integrals
    double ao_sparsity() { return (1.0 - (double)small_skips_[nbf_] / (double)(nbf_ * nbf_)); }

    ///
    /// Sets the AO integrals to in-core. (Defaults to TRUE)
    /// @param core True to indicate in-core
    /// DFHelper will keep track of this memory.
    /// NOTE: DFHelper will automatically revert to on-disk if the
    /// sizes of the AO integrals is greater than 90% of the memory
    /// it controls.
    ///
    void set_AO_core(bool core) { AO_core_ = core; }
    bool get_AO_core() { return AO_core_; }

    /// Sets the flag to release early the memory storing in-core AO integrals.
    /// Methods where the AO integrals are not needed in the
    /// metric contraction step, like SAPT(DFT), can call this to free
    /// up memory for storing intermediates during metric contraction.
    void set_release_core_AO_before_metric(bool release) { release_core_AO_before_metric_ = release; }

    ///
    /// Sets the MO integrals to in-core. (Defaults to FALSE)
    /// @param core True to indicate in-core
    /// DFHelper will not keep track of this memory.
    /// If a seg fault occurs, the MOs were bigger than you thought!
    ///
    void set_MO_core(bool core) { MO_core_ = core; }
    bool get_MO_core() { return MO_core_; }

    /// schwarz screening cutoff (defaults to 1e-12)
    void set_schwarz_cutoff(double cutoff) { cutoff_ = cutoff; }
    double get_schwarz_cutoff() { return cutoff_; }

    /// fitting metric power (defaults to -0.5) to use in
    /// K_{m n} = C_{l a}(m l|Q)(Q|R)^{-1/2}(R|P)^{-1/2}(P|n s)C_{s a}
    void set_metric_pow(double m_pow) { mpower_ = m_pow; }
    double get_metric_pow() { return mpower_; }

    /// fitting metric power for w integrals (defaults to -1.0) to use in
    /// wK_{m n} = C_{l a}(m l|Q)(Q|P)^-1(P|w|n s)C_{s a}
    void set_wmetric_pow(double wm_pow) { wmpower_ = wm_pow; }
    double get_wmetric_pow() { return wmpower_; }

    ///
    /// Sets the fitting metric to be held in core (defaults to FALSE)
    /// @param hold TRUE if metrics are to be held in core.
    /// Memory contraints adapt accordingly.
    ///
    void hold_met(bool hold) { hold_met_ = hold; }
    bool get_hold_met() { return hold_met_; }

    ///
    /// Sets the fitting metric condition
    /// @param condition: tolerrence for metric^pow
    ///
    void set_fitting_condition(double condition) { condition_ = condition; }
    bool get_fitting_condition() { return condition_; }


    ///
    /// Do we calculate omega exchange and regular hf exchange together?
    /// @param wcombine boolean: all exchange in one matrix
    /// TODO: re-enable after all bugs are fixed.
    // void set_wcombine(bool wcombine) {wcombine_ = wcombine;}
    void set_wcombine(bool wcombine) {
        if (wcombine) {
            throw PSIEXCEPTION("JK: wcombine option is currently not available.");
        }
        wcombine_ = wcombine;
    }
    bool get_wcombine() { return wcombine_; }

    ///
    /// Lets me know whether to compute those other type of integrals
    /// @param do_wK boolean indicating to compute other integrals
    ///
    void set_do_wK(bool do_wK) { do_wK_ = do_wK; }
    size_t get_do_wK() { return do_wK_; }

    ///
    /// sets the parameter for the other type of integrals
    /// @param omega double indicating parameter for other type
    ///
    void set_omega(double omega) { omega_ = omega; }
    double get_omega() { return omega_; }

    ///
    /// sets the coefficient for (pq|rs) integrals
    /// @param omega double indicating coefficient for eri
    ///
    void set_omega_alpha(double alpha) { omega_alpha_ = alpha; }
    double get_omega_alpha() { return omega_alpha_; }
    
    ///
    /// sets the parameter for the other type of integrals
    /// @param omega double indicating parameter for other type
    ///
    void set_omega_beta(double beta) { omega_beta_ = beta; }
    double get_omega_beta() { return omega_beta_; }

    ///
    /// set the printing verbosity parameter
    /// @param print_lvl indicating verbosity
    ///
    void set_print_lvl(int print_lvl) { print_lvl_ = print_lvl; }
    int get_print_lvl() { return print_lvl_; }

    /// Initialize the object
    void initialize();

    /// Prepare screening and indexing metadata used for Schwarz screening.
    void prepare_sparsity();

    /// print tons of useful info
    void print_header();

    ///
    /// Add transformation space with key
    /// @param key used to access space orbitals
    /// @param space orbital space matrix
    ///
    void add_space(std::string key, SharedMatrix space);

    ///
    /// Add transformation with name using two space keys
    /// @param name used to access transformed integrals
    /// @param key1 left oribtal space
    /// @param key2 right oribtal space
    /// @param order for direct builds of "Qpq" and "pqQ" forms only
    ///
    void add_transformation(std::string name, std::string key1, std::string key2, std::string order = "Qpq");

    /// invoke transformations
    void transform();

    // => Tensor IO <=
    // many ways to access the 3-index tensors.

    ///
    /// Fill a SharedMatrix with three index pairs.  Slice the same way you do in python.
    /// @param name name of transformation to be accessed
    /// @param M SharedMatrix you want to be filled
    /// Recursive signitures were added if you want a full 3rd index, 2nd and 3rd index, etc.
    /// For example, fill_tensor("ia", M, (0, 15)) will get you ia[0:15, :, :]
    /// I will check to make sure your slice sizes are not larger than the matrix bounds,
    /// but be prepared for a runtime throw.
    ///
    void fill_tensor(std::string name, SharedMatrix M);
    void fill_tensor(std::string name, SharedMatrix M, std::vector<size_t> a1);
    void fill_tensor(std::string name, SharedMatrix M, std::vector<size_t> a1, std::vector<size_t> a2);
    void fill_tensor(std::string name, SharedMatrix M, std::vector<size_t> a1, std::vector<size_t> a2,
                     std::vector<size_t> a3);

    ///
    /// Fill a buffer with three index pairs. Same concept as fill_tensor.
    /// @param name name of transformation to be accessed
    /// @param b pointer to allocated memory
    /// Be cautious, this function does not bound check your buffer against the tuples provied.
    ///
    void fill_tensor(std::string name, double* b, std::vector<size_t> a1, std::vector<size_t> a2,
                     std::vector<size_t> a3);
    void fill_tensor(std::string name, double* b, std::vector<size_t> a1, std::vector<size_t> a2);
    void fill_tensor(std::string name, double* b, std::vector<size_t> a1);
    void fill_tensor(std::string name, double* b);

    ///
    /// return a SharedMatrix, I take care of sizing for you.
    /// @param name name of transformation to be accessed
    /// I always compound the 2nd and 3rd indices.
    /// For example, get_tensor("ia", (0:15), (0:5), (0:5)) will return a
    /// SharedMatrix of size (15, 25), so be careful if you plan to use Matrix::gemm
    ///
    SharedMatrix get_tensor(std::string name);
    SharedMatrix get_tensor(std::string name, std::vector<size_t> a1);
    SharedMatrix get_tensor(std::string name, std::vector<size_t> a1, std::vector<size_t> a2);
    SharedMatrix get_tensor(std::string name, std::vector<size_t> a1, std::vector<size_t> a2, std::vector<size_t> a3);

    ///
    /// Add a 3-index disk tensor (that is not a transformation)
    /// @param name name of tensor - used to be accessed later
    /// @param dimensions shape of tensor
    ///
    void add_disk_tensor(std::string key, std::tuple<size_t, size_t, size_t> dimensions);

    ///
    /// Write to a 3-index disk tensor from a SharedMatrix
    /// Can be a transformation if you want to overwrite one.
    /// @param name name of tensor - used to be accessed later
    /// @param M SharedMatrix with contents to write to disk tensor
    /// Bound checks are performed
    ///
    void write_disk_tensor(std::string name, SharedMatrix M);
    void write_disk_tensor(std::string name, SharedMatrix M, std::vector<size_t> a1);
    void write_disk_tensor(std::string name, SharedMatrix M, std::vector<size_t> a1, std::vector<size_t> a2);
    void write_disk_tensor(std::string name, SharedMatrix M, std::vector<size_t> a1, std::vector<size_t> a2,
                           std::vector<size_t> a3);

    ///
    /// Write to a 3-index disk tensor from a buffer
    /// Can be a transformation if you want to overwrite one.
    /// @param name name of tensor - used to be accessed later
    /// @param b buffer to write to disk tensor
    /// No bound checking on buffer is performed, use caution!
    ///
    void write_disk_tensor(std::string name, double* b, std::vector<size_t> a1, std::vector<size_t> a2,
                           std::vector<size_t> a3);
    void write_disk_tensor(std::string name, double* b, std::vector<size_t> a1, std::vector<size_t> a2);
    void write_disk_tensor(std::string name, double* b, std::vector<size_t> a1);
    void write_disk_tensor(std::string name, double* b);

    /// tranpose a tensor *after* it has been written
    void transpose(std::string name, std::tuple<size_t, size_t, size_t> order);

    /// clear spaces
    void clear_spaces();

    /// clear transformations
    void clear_transformations(); 

    /// clears spaces and transformations
    void clear_all();

    /// get sizes, shapes of tensors
    size_t get_space_size(std::string key);
    size_t get_tensor_size(std::string key);
    std::tuple<size_t, size_t, size_t> get_tensor_shape(std::string key);
    size_t get_naux() { return naux_; }

    /// builds J/K
    void build_JK(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright, std::vector<SharedMatrix> D,
                  std::vector<SharedMatrix> J, std::vector<SharedMatrix> K, std::vector<SharedMatrix> wK,
                  size_t max_nocc, bool do_J, bool do_K, bool do_wK, bool lr_symmetric);

   protected:
    // => basis sets <=
    std::shared_ptr<BasisSet> primary_;
    std::shared_ptr<BasisSet> aux_;
    size_t nbf_;
    size_t naux_;

    // => memory in doubles <=
    size_t memory_ = 256000000;
    size_t required_core_size_;

    // => internal holders <=
    std::string method_ = "STORE";
    bool direct_ = false;
    bool direct_iaQ_ = false;
    bool symm_compute_;
    // Use in-core subalgorithm?
    bool AO_core_ = true;
    bool MO_core_ = false;
    bool release_core_AO_before_metric_ = false;
    size_t nthreads_ = 1;
    double cutoff_ = 1e-12;
    double condition_ = 1e-12;
    double mpower_ = -0.5;
    double wmpower_ = -1.0;
    bool hold_met_ = false;
    bool built_ = false;
    bool transformed_ = false;
    std::pair<size_t, size_t> info_;
    bool ordered_ = false;
    bool do_wK_ = false;
    bool wcombine_ = false;
    double omega_;
    double omega_alpha_;
    double omega_beta_;
    bool debug_ = false;
    // Can we early-exit prepare_sparsity?
    bool sparsity_prepared_ = false;
    int print_lvl_ = 1;
    // Which subalgorithm (in-core or out-of-core) do we use?
    // AUTO lets DFHelper decide based on memory allocated
    // INCORE forces the in-core subalgorithm (AO_core_ = true)
    // OUT_OF_CORE forces the out-of-core subalgorithm (AO_core_ = false)
    std::string subalgo_ = "AUTO";

    // => in-core machinery <=
    void AO_core(bool set_AO_core);
    std::unique_ptr<double[]> Ppq_;
    // Maps x -> (P|Q) ^ x.
    std::map<double, SharedMatrix> metrics_;

    // => in-core wK machinery <=
    std::unique_ptr<double[]> wPpq_;  // if do_wK_ holds (A|w|mn)
    std::unique_ptr<double[]> m1Ppq_;

    // => AO building machinery <=
    void prepare_AO();
    void prepare_AO_core();
    void compute_dense_Qpq_blocking_Q(const size_t start, const size_t stop, double* Mp,
                                      std::vector<std::shared_ptr<TwoBodyAOInt>> eri);
    void compute_sparse_pQq_blocking_Q(const size_t start, const size_t stop, double* Mp,
                                       std::vector<std::shared_ptr<TwoBodyAOInt>> eri);
    void compute_sparse_pQq_blocking_p(const size_t start, const size_t stop, double* Mp,
                                       std::vector<std::shared_ptr<TwoBodyAOInt>> eri);
    void compute_sparse_pQq_blocking_p_symm(const size_t start, const size_t stop, double* Mp,
                                            std::vector<std::shared_ptr<TwoBodyAOInt>> eri);
    void compute_sparse_pQq_blocking_p_symm_abw(const size_t start, const size_t stop, double* just_Mp, double* param_Mp,
                                        std::vector<std::shared_ptr<TwoBodyAOInt>> eri,
                                        std::vector<std::shared_ptr<TwoBodyAOInt>> weri);



    void contract_metric_AO_core_symm(double* Qpq, double* Ppq, double* metp, size_t begin, size_t end);
    void grab_AO(const size_t start, const size_t stop, double* Mp);

    // => wK AO building machinery <=
    void prepare_AO_wK_core();
    void prepare_AO_wK();

    void copy_upper_lower_wAO_core_symm(double* Qpq, double* Ppq, size_t begin, size_t end);

    // first integral transforms
    void first_transform_pQq(size_t bsize, size_t bcount, size_t block_size, double* Mp, double* Tp, double* Bp,
                             std::vector<std::vector<double>>& C_buffers);

    // => index vectors for screened AOs <=
    // ==> Skips for non-symmetric "densities" <==
    // For a given basis function, how many basis functions have significant interactions?
    // Used for "q" indexing.
    std::vector<size_t> small_skips_;
    // For a given basis function, how many significant (basis function, naux) pairs have we already passed?
    // This is a running total, unlike small_skips. Used for "p" indexing.
    std::vector<size_t> big_skips_;
    // ==> Skips for symmetric "densities" <==
    // For a given basis function, how many basis functions have significant interaction AND have lesser index?
    // This counts columns we ignore due to the density of the symmetry. Used for "q" indexing.
    std::vector<size_t> symm_ignored_columns_;
    // For a given basis function, how many basis functions have significant interaction AND have non-lesser index?
    // Should be small_skips[p] - symm_ginored_columns[p]. Used for "q" indexing.
    std::vector<size_t> symm_small_skips_;
    // For a given basis function, how many significant (non-lesser index basis function, naux)
    // pairs have we already passed? This is a running total, unlike symm_small_skips.
    // Used for "p" indexing.
    std::vector<size_t> symm_big_skips_;

    // => Shell info and blocking : no screening involved <=
    // Number of primary shells
    size_t pshells_;
    // Number of auxiliary shells
    size_t Qshells_;
    // greatest number of functions in an aux_ shell
    double Qshell_max_;
    // When we start primary shell i, how many functions have we passed? Length pshells_ + 1.
    std::vector<size_t> pshell_aggs_;
    // When we start auxilliary shell i, how many functions have we passed? Length Qshells_ + 1.
    std::vector<size_t> Qshell_aggs_;
    // Populate all variables in this section. Should only ever be called by the constructor.
    void prepare_blocking();

    // => generalized blocking <=
    std::pair<size_t, size_t> pshell_blocks_for_AO_build(const size_t mem, size_t symm,
                                                         std::vector<std::pair<size_t, size_t>>& b);
    // returns pair(largest buffer size, largest block size)
    std::pair<size_t, size_t> Qshell_blocks_for_transform(const size_t mem, size_t wtmp, size_t wfinal,
                                                          std::vector<std::pair<size_t, size_t>>& b);
    void metric_contraction_blocking(std::vector<std::pair<size_t, size_t>>& steps, size_t blocking_index,
                                     size_t block_sizes, size_t total_mem, size_t memory_factor, size_t memory_bump);

    // => Schwarz Screening <=
    // Does shell pair (m, n) have any significant (mn|mn)-type integrals? Size pshells_ ** 2
    std::vector<bool> schwarz_shell_mask_;
    // What is the index of ij in significant basis function interactions for i? 0 = not significant. 2 means 1 before, 3 means 2 before, ect. Size nbf_ ** 2.
    // Ordering is based on the natural ordering of basis functions, not ordering of significance.
    std::vector<size_t> schwarz_fun_index_;

    // => Coulomb metric handling <=
    std::vector<std::pair<double, std::string>> metric_keys_;
    // Create J and write it to disk.
    void prepare_metric();
    // Create J and cache it in metrics_.
    void prepare_metric_core();
    double* metric_prep_core(double m_pow);
    std::string return_metfile(double m_pow);
    std::string compute_metric(double m_pow);

    double* metric_inverse_prep_core();

    // => metric operations <=
    void contract_metric_Qpq(std::string file, double* metp, double* Mp, double* Fp, const size_t tots);
    void contract_metric(std::string file, double* metp, double* Mp, double* Fp, const size_t tots);
    void contract_metric_core(std::string file);
    void contract_metric_AO(double* Mp);
    void contract_metric_AO_core(double* Qpq, double* metp);

    // => spaces and transformation maps <=
    std::map<std::string, std::tuple<SharedMatrix, size_t>> spaces_;
    std::map<std::string, std::tuple<std::string, std::string, size_t>> transf_;
    std::map<std::string, std::unique_ptr<double[]>> transf_core_;

    // => transformation machinery <=
    std::pair<size_t, size_t> identify_order();
    void print_order();
    void put_transformations_Qpq(int begin, int end, int wsize, int bsize, double* Fp, int ind, bool bleft);
    void put_transformations_pQq(int begin, int end, int block_size, int bcount, int wsize, int bsize, double* Np,
                                 double* Fp, int ind, bool bleft);
    std::vector<std::pair<std::string, size_t>> sorted_spaces_;
    std::vector<std::string> order_;
    std::vector<std::string> bspace_;
    std::vector<size_t> strides_;

    // => FILE IO maintenence <=
    typedef struct StreamStruct {
        StreamStruct();
        StreamStruct(std::string filename, std::string op, bool activate = true);
        ~StreamStruct();

        FILE* get_stream(std::string op);
        void change_stream(std::string op);
        void close_stream();

        FILE* fp_;
        std::string op_;
        bool open_ = false;
        std::string filename_;

    } Stream;

    std::map<std::string, std::shared_ptr<Stream>> file_streams_;
    FILE* stream_check(std::string filename, std::string op);

    // => FILE IO machinery <=
    void put_tensor(std::string file, double* b, std::pair<size_t, size_t> a1, std::pair<size_t, size_t> a2,
                    std::pair<size_t, size_t> a3, std::string op);

    void put_tensor(std::string file, double* b, const size_t start1, const size_t stop1, const size_t start2,
                    const size_t stop2, std::string op);
    void get_tensor_(std::string file, double* b, std::pair<size_t, size_t> a1, std::pair<size_t, size_t> a2,
                     std::pair<size_t, size_t> a3);
    void get_tensor_(std::string file, double* b, const size_t start1, const size_t stop1, const size_t start2,
                     const size_t stop2);
    // Write to `file` from `Mp`, starting at position `start` in `file` and reading length `size`.
    // The file is opened in mode `op`.
    void put_tensor_AO(std::string file, double* Mp, size_t size, size_t start, std::string op);
    // Read from `file` into `Mp`, starting at position `start` in `file` and reading length `size`
    void get_tensor_AO(std::string file, double* Mp, size_t size, size_t start);

    // => internal handlers for FILE IO <=
    // Map a base filename to the fullname of the p variant and then the original tensor.
    std::map<std::string, std::tuple<std::string, std::string>> files_;
    // Map a full filename to the size of its three indices.
    std::map<std::string, std::tuple<size_t, size_t, size_t>> sizes_;
    std::map<std::string, std::tuple<size_t, size_t, size_t>> tsizes_;
    std::vector<size_t> AO_file_sizes_;
    std::vector<std::string> AO_names_;
    // Given the base of a filename, return the full filename.
    // "full" filename adds parts to the base to prevent files from overwriting each other.
    std::string start_filename(std::string start);
    // Given a base filename and the indices of its three dimensions, register the
    // filename in files_ and the dimensions in sizes_. Creates an entry for both
    // the tensor and the p variant. (The distinction is important for direct_iaQ.)
    // Q = aux, p = left primary, q = right primary
    // op = 0 if Qpq, 1 if pQq, 2 if pqQ
    void filename_maker(std::string name, size_t Q, size_t p, size_t q, size_t op = 0);
    // Store a filename for the i'th AO quantity in AO_names_.
    void AO_filename_maker(size_t i);
    void check_file_key(std::string);
    void check_file_tuple(std::string name, std::pair<size_t, size_t> t0, std::pair<size_t, size_t> t1,
                          std::pair<size_t, size_t> t2);
    void check_matrix_size(std::string name, SharedMatrix M, std::pair<size_t, size_t> t0, std::pair<size_t, size_t> t1,
                           std::pair<size_t, size_t> t2);

    // => transpose a tensor <=
    void transpose_core(std::string name, std::tuple<size_t, size_t, size_t> order);
    void transpose_disk(std::string name, std::tuple<size_t, size_t, size_t> order);

    // => JK <=
    void compute_JK(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright, std::vector<SharedMatrix> D,
                    std::vector<SharedMatrix> J, std::vector<SharedMatrix> K, size_t max_nocc, bool do_J, bool do_K,
                    bool do_wK, bool lr_symmetric);
    void compute_D(std::vector<SharedMatrix> D, std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright);
    // Compute the J matrices (mn|rs) D_rs and store them into J, for an array of D.
    // @param D: vector of densities to be contracted against AO basis integrals.
    // @param J: vector of Coulomb integrals, one for each density
    // @param M1p : Intermediate populated now for wK later.
    // @param T1p : Temporary matrix
    // @param T2p : Temporary matrix
    void compute_J(const std::vector<SharedMatrix> D, std::vector<SharedMatrix> J, double* Mp, double* T1p, double* T2p,
                   std::vector<std::vector<double>>& D_buffers, size_t bcount, size_t block_size);
    void compute_J_symm(std::vector<SharedMatrix> D, std::vector<SharedMatrix> J, double* Mp, double* T1p, double* T2p,
                        std::vector<std::vector<double>>& D_buffers, size_t bcount, size_t block_size);
    void compute_J_combined(std::vector<SharedMatrix> D, std::vector<SharedMatrix> J, double* Mp, double* T1p, double* T2p, std::vector<std::vector<double>>& D_buffers, size_t bcount, size_t block_size);
    void compute_K(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright, std::vector<SharedMatrix> K,
                   double* Tp, double* Jtmp, double* Mp, size_t bcount, size_t block_size,
                   std::vector<std::vector<double>>& C_buffers, bool lr_symmetric);
    // returns tuple(largest AO buffer size, largest Q block size)
    std::tuple<size_t, size_t> Qshell_blocks_for_JK_build(std::vector<std::pair<size_t, size_t>>& b, size_t max_nocc,
                                                          bool lr_symmetric);
    void compute_wK(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright, std::vector<SharedMatrix> wK,
                    size_t max_nocc, bool do_J, bool do_K, bool do_wK);

    // => misc <=
    // Utility function to fill double* with zero in parallel
    void fill(double* b, size_t count, double value);

};  // End DF Helper class
}  // psi4 namespace
#endif
