/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

#ifndef three_index_df_helper
#define three_index_df_helper

#include "psi4/psi4-dec.h"
#include <psi4/libmints/typedefs.h>

#include <map>
#include <list>
#include <vector>
#include <tuple>
#include <string>

namespace psi {

class BasisSet;
class Options;
class Matrix;
class ERISieve;
class TwoBodyAOInt;

namespace df_helper{

class DF_Helper {

public:

    DF_Helper(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> aux);
    ~DF_Helper();

    // => user options - set these before calling initialize() <=
    // workflow method (store or direct) (defaults to STORE)
    void set_method(std::string met) { method_ = met;}
    std::string get_method() { return method_; }
    
    // number of threads (defaults to 1)
    void set_nthreads(size_t threads) { nthreads_ = threads;}
    size_t get_nthreads() { return nthreads_; }
    
    // memory in doubles (defaults to 256,000,000) (2.048GB)
    void set_memory(size_t mem) { memory_ = mem;}
    size_t get_memory() { return memory_; }

    // I can tell you how many doubles the *screened* AOs will require
    size_t get_AO_size() { return big_skips_[nao_]; } 

    // want the AO integrals in core? (defaults to True) 
    // (I will adapt blocking  and turn off if necessary to keep memory satisfied)
    void set_AO_core(bool on) {AO_core_ = on;}
    bool get_AO_core() { return AO_core_; }
    
    // Do you want the MO integrals in core? (defaults to FALSE) (not my responsiblity to keep track of)
    void set_MO_core(bool on) {MO_core_ = on;}
    bool get_MO_core() { return MO_core_; }
    
    // schwarz cutoff (defaults to 1e-12)
    void set_schwarz_cutoff(double cutoff) { cutoff_ = cutoff;}
    double get_schwarz_cutoff() { return cutoff_; }
    
    // metric power (defaults to -0.5)  
    void set_metric_pow(double pow) { mpower_ = pow;}
    double get_metric_pow() { return mpower_;}

    // Want the metric to be held in core? (defaults to FALSE) (I will adapt blocking to keep memory satisfied)
    void hold_met(bool hold) {hold_met_ = hold;}
    bool get_hold_met() { return hold_met_;}

    // Tell me the worst MO size to improve blocking. (defaults to 0.5*AO_SIZE, must specify if greater!)
    void set_MO_hint(size_t wMO) {wMO_ = wMO;}
    size_t get_MO_hint() { return wMO_;}
    
    // Enhanced memory use if (for every Cleft = Cright) in JK builds. (defaults to FALSE!)
    void set_JK_hint(bool hint) {JK_hint_ = hint;}
    size_t get_JK_hint() { return JK_hint_;}
    // => end user options <=

    // Initialize the object
    void initialize();
    void print_header();

    // Add transformation space with key
    void add_space(std::string key, SharedMatrix M);

    // add transformation with name using two space keys
    void add_transformation(std::string name, std::string key1, std::string key2, std::string order = "Qpq");

    // invoke transformations
    void transform();

    // => Tensor IO <=
    // Fill a SharedMatrix with three index pairs.  Slice the same way you do in python. 
    // Recursive signitures were added if you want a full 3rd index, 2nd and 3rd index, etc.
    // For example, fill_tensor("ia", M, (0, 15)) will get you ia[0:15, :, :]
    // I will check to make sure your slice sizes are not larger than the matrix bounds, 
    // but be prepared for a runtime throw.
    void fill_tensor(std::string name, SharedMatrix M);
    void fill_tensor(std::string name, SharedMatrix M, std::pair<size_t, size_t> a1);
    void fill_tensor(std::string name, SharedMatrix M, std::pair<size_t, size_t> a1,
      std::pair<size_t, size_t> a2);
    void fill_tensor(std::string name, SharedMatrix M, std::pair<size_t, size_t> a1,
      std::pair<size_t, size_t> a2, std::pair<size_t, size_t> a3);

    // return a SharedMatrix, I take care of sizing for you.
    // I always compound the 2nd and 3rd indices.
    // For example, get_tensor("ia", (0:15), (0:5), (0:5)) will return a 
    // SharedMatrix of size (15, 25), so be careful if you plan to use Matrix::gemm
    SharedMatrix get_tensor(std::string name);
    SharedMatrix get_tensor(std::string name, std::pair<size_t, size_t> a1);
    SharedMatrix get_tensor(std::string name, std::pair<size_t, size_t> a1,
      std::pair<size_t, size_t> a2);
    SharedMatrix get_tensor(std::string name, std::pair<size_t, size_t> a1,
      std::pair<size_t, size_t> a2, std::pair<size_t, size_t> a3);

    // Add a 3-index disk tensor, write to it, or write over a transformed tensor
    void add_disk_tensor(std::string key, std::tuple<size_t, size_t, size_t> dimensions);
    void write_disk_tensor(std::string key, SharedMatrix M);    
    void write_disk_tensor(std::string key, SharedMatrix M, std::pair<size_t, size_t> a0); 
    void write_disk_tensor(std::string key, SharedMatrix M, std::pair<size_t, size_t> a0, 
    std::pair<size_t, size_t> a1);
    void write_disk_tensor(std::string key, SharedMatrix M, std::pair<size_t, size_t> a0, 
        std::pair<size_t, size_t> a1, std::pair<size_t, size_t> a2);    

    // only use this one if you now what you're doing -- I do not bound check!
    void fill_tensor(std::string name, double * b, std::pair<size_t, size_t> a1,
      std::pair<size_t, size_t> a2, std::pair<size_t, size_t> a3);
    void fill_tensor(std::string name, double * b, std::pair<size_t, size_t> a1,
      std::pair<size_t, size_t> a2);
    void fill_tensor(std::string name, double * b, std::pair<size_t, size_t> a1);
    void fill_tensor(std::string name, double * b); 
    
    // tranpose a tensor *after* it has been written
    void transpose(std::string name, std::tuple<size_t, size_t, size_t> order);
    
    // clear spaces or spaces and transformations
    void clear_spaces();
    void clear_all();

    // get sizes, shapes
    size_t get_space_size(std::string key);
    size_t get_tensor_size(std::string key);
    std::tuple<size_t, size_t, size_t> get_tensor_shape(std::string key);
    size_t get_naux() { return naux_; }

    // => Build JK <=
    void build_JK(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright, 
        std::vector<SharedMatrix> J, std::vector<SharedMatrix> K); 

protected:

    // => basis sets <=
    std::shared_ptr<BasisSet> primary_;
    std::shared_ptr<BasisSet> aux_;
    size_t nao_;
    size_t naux_;

    // => memory in doubles <=
    size_t memory_ = 256000000;

    // => internal holders <=
    std::string method_ = "STORE";
    bool direct_;
    bool AO_core_ = 1;
    bool MO_core_ = 0;
    size_t nthreads_ = 1;
    double cutoff_ = 1e-12;
    double tolerance_ = 0.0;
    double condition_ = 1e-12;
    double mpower_ = -0.5;
    bool hold_met_ = false;
    bool JK_hint_ = false;
    bool built = false;    
    bool transformed_ = false;
    std::pair<size_t, size_t> info_;
    bool ordered_=0;
    std::pair<size_t, size_t> identify_order();
    void print_order();
    size_t wMO_;

    // => in-core machinery <=
    std::vector<double> Ppq_;
    std::map<double, SharedMatrix> metrics_;

    // => AO building machinery <=
    void prepare_AO();
    void prepare_AO_core();
    void compute_AO_Q(const size_t start, const size_t stop, double* Mp, std::vector<std::shared_ptr<TwoBodyAOInt>> eri);
    void compute_AO_p(const size_t start, const size_t stop, double* Mp, std::vector<std::shared_ptr<TwoBodyAOInt>> eri);
    void compute_AO_p_symm(const size_t start, const size_t stop, double* Mp, std::vector<std::shared_ptr<TwoBodyAOInt>> eri);
    void contract_metric_AO_core_symm(double* Qpq, double* metp, size_t begin, size_t end);
    void grab_AO(const size_t start, const size_t stop, double* Mp);

    // => index vectors for screened AOs <=
    std::vector<size_t> small_skips_;
    std::vector<size_t> big_skips_;
    std::vector<size_t> symm_skips_;
    std::vector<size_t> symm_sizes_;
    std::vector<size_t> symm_agg_sizes_;

    // => shell info and blocking <=
    size_t pshells_;
    size_t Qshells_;
    std::vector<size_t> pshell_aggs_;
    std::vector<size_t> Qshell_aggs_;
    void prepare_blocking();

    // => generalized blocking <=
    std::pair<size_t, size_t> pshell_blocks_for_AO_build(const size_t mem, size_t symm, std::vector<std::pair<size_t, size_t>>& b);
    std::pair<size_t, size_t> Qshell_blocks_for_transform(const size_t mem, size_t wtmp, size_t wfinal, std::vector<std::pair<size_t, size_t>>& b);

    // => Schwarz Screening <=
    std::vector<size_t> schwarz_fun_mask_;
    std::vector<size_t> schwarz_shell_mask_;
    std::vector<size_t> schwarz_fun_count_;
    void prepare_sparsity();

    // => Coulomb metric handling <=
    std::vector<std::pair<double, std::string>> metric_keys_;
    void prepare_metric();
    void prepare_metric_core();
    double* metric_prep_core(double pow);
    std::string return_metfile(double pow);
    std::string compute_metric(double pow);

    // => metric operations <=
    void contract_metric(std::string file, double* metp, double* Mp, double* Fp, const size_t tots);
    void contract_metric_core(std::string file);
    void contract_metric_AO(double* Mp);
    void contract_metric_AO_core(double* Qpq, double* metp);

    // => spaces and transformation maps <=
    std::map<std::string, std::tuple<SharedMatrix, size_t>> spaces_;
    std::map<std::string, std::tuple<std::string, std::string, size_t>> transf_;
    std::map<std::string, std::vector<double>> transf_core_;

    // => transformation machinery <=
    void transform_core();
    void transform_disk();
    std::vector<std::pair<std::string, size_t>> sorted_spaces_;
    std::vector<std::string> order_;
    std::vector<std::string> bspace_;
    std::vector<size_t> strides_;
    
    // => FILE IO maintenence <= 
    struct stream{
        FILE* fp;
        std::string op;};
    std::map<std::string, stream> file_status_;
    FILE* stream_check(std::string filename, std::string op);

    // => FILE IO machinery <= 
    void put_tensor(std::string file, double* b, std::pair<size_t, size_t> a1,
      std::pair<size_t, size_t> a2, std::pair<size_t, size_t> a3, std::string op);
    void put_tensor(std::string file, double* b, const size_t start1,
      const size_t stop1, const size_t start2, const size_t stop2, std::string op);
    void get_tensor_(std::string file, double* b, std::pair<size_t, size_t> a1,
      std::pair<size_t, size_t> a2, std::pair<size_t, size_t> a3);
    void get_tensor_(std::string file, double* b,  const size_t start1,
      const size_t stop1, const size_t start2, const size_t stop2);
    void put_tensor_AO(std::string file, double* Mp, size_t size, size_t start, std::string op);
    void get_tensor_AO(std::string file, double* Mp, size_t size, size_t start);

    // => internal handlers for FILE IO <= 
    std::map<std::string, std::tuple<std::string, std::string>> files_;
    std::map<std::string, std::tuple<size_t, size_t, size_t>> sizes_;
    std::map<std::string, std::tuple<size_t, size_t, size_t>> tsizes_;
    std::map<std::string, std::string> AO_files_;
    std::vector<size_t> AO_file_sizes_;
    std::vector<std::string> AO_names_;
    void filename_maker(std::string name, size_t a0, size_t a1, size_t a2, size_t op = 0);
    void AO_filename_maker(size_t i);
    void check_file_key(std::string);
    void check_file_tuple(std::string name, std::pair<size_t, size_t> t0, 
        std::pair<size_t, size_t> t1, std::pair<size_t, size_t> t2);
    void check_matrix_size(std::string name, SharedMatrix M, std::pair<size_t, size_t> t0,
        std::pair<size_t, size_t> t1, std::pair<size_t, size_t> t2);

    // => transpose a tensor <= 
    void transpose_core(std::string name, std::tuple<size_t, size_t, size_t> order);
    void transpose_disk(std::string name, std::tuple<size_t, size_t, size_t> order);

    // => JK <=
    void compute_JK(std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright, 
        std::vector<SharedMatrix> J, std::vector<SharedMatrix> K); 
    void compute_D(std::vector<SharedMatrix>& D, std::vector<SharedMatrix> Cleft, std::vector<SharedMatrix> Cright);
    void compute_J(std::vector<SharedMatrix> D, std::vector<SharedMatrix> J, double* Mp, 
        double* T1p, double* T2p, std::vector<std::vector<double>> D_buffers, size_t bcount, size_t block_size);
    void compute_J_symm(std::vector<SharedMatrix> D, std::vector<SharedMatrix> J, double* Mp, 
        double* T1p, double* T2p, std::vector<std::vector<double>> D_buffers, size_t bcount, size_t block_size);
    void compute_K(std::vector<SharedMatrix> Cleft, 
        std::vector<SharedMatrix> Cright, std::vector<SharedMatrix> K, double* Tp, double* Jtmp,    
        double* Mp, size_t bcount, size_t block_size, std::vector<std::vector<double>> C_buffers, std::vector<SharedMatrix> D, std::vector<SharedMatrix> J);
    std::tuple<size_t,size_t,size_t,size_t> Qshell_blocks_for_JK_build(std::vector<std::pair<size_t, size_t>>& b, std::vector<SharedMatrix> Cleft, 
        std::vector<SharedMatrix> Cright);

}; // End DF Helper class

}}//end df_helper/psi4 namespace
#endif
