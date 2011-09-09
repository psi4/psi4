#ifndef _psi_src_lib_libmints_matrix_block_h_
#define _psi_src_lib_libmints_matrix_block_h_

#include <cstdio>
#include <string>
#include <cstring>

#include <libparallel/parallel.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>

namespace boost {
template<class T> class shared_ptr;
// Forward declarations for boost.python used in the extract_subsets

namespace python{
       class tuple;
}}


namespace psi {


class Block
{
protected:
    /// This is a block that is further divided into subblocks
    std::vector<double> block_;

    /// The size of the block;
    int nelements_;

    /// The number of sub_blocks;
    int nsb_;

    /// The number of rows in the entire distributed matrix
    int g_nrows_;
    /// The number of columns in the entire distributed matrix
    int g_ncols_;

    /// The number of rows in the block
    int b_nrows_;
    /// The number of columns in the block
    int b_ncols_;

    /// The offset to the sub_blocks;
    std::vector<int> sb_offset_;

    /// The number of rows for each sub_block;
    std::vector<int> sb_nrows_;

    /// The number of columns for each sub_block;
    std::vector<int> sb_ncols_;

    /// The size of each sub_block;
    std::vector<int> sb_size_;

    /// Maps the global sub_block values to the local sub_block value
    std::map<int, int> global_local_sblock_;
    /// Maps the local sub_block value to the global sub_block value
    std::vector<int> local_global_sblock_;

    /// The block name used for printing
    std::string name_;

    /// Process ID
    int me_;
    /// Number of MPI Processes
    int nprocs_;
    /// Number of threads
    int nthreads_;
    /// Communicator Type
    std::string comm_;

    template<typename T>
    void free_vector(std::vector<T> &vec)
    {
        std::vector<T> tmp;
        vec.clear();
        vec.swap(tmp);
    }
    template<typename T1, typename T2>
    void free_map(std::map<T1,T2> &map)
    {
        std::map<T1,T2> tmp;
        map.clear();
        map.swap(tmp);
    }

public:

    Block() : nelements_(0), nsb_(0), b_nrows_(0), b_ncols_(0),
        g_nrows_(0), g_ncols_(0)
    {
        me_ = Communicator::world->me();
        nprocs_ = Communicator::world->nproc();
        nthreads_ = Communicator::world->nthread();
        comm_ = Communicator::world->communicator();
    }

    ~Block() { Communicator::world->sync(); }

    void clear_block()
    {
        nelements_ = 0;
        nsb_ = 0;
        b_nrows_ = 0;
        b_ncols_ = 0;
        g_nrows_ = 0;
        g_ncols_ = 0;
        name_.clear();
        free_vector(block_);
        free_vector(sb_offset_);
        free_vector(sb_nrows_);
        free_vector(sb_ncols_);
        free_vector(sb_size_);
        free_vector(local_global_sblock_);
        free_map(global_local_sblock_);
    }

    /// Returns the number of elements in the block
    int nelements() const { return nelements_; }
    /// Return the number sub_blocks in the block
    int nsb() const { return nsb_; }
    /// Return the global number of rows
    int g_nrows() const { return g_nrows_; }
    /// Return the global number of columns
    int g_ncols() const { return g_ncols_; }
    /// Return the number of rows in the block
    int b_nrows() const { return b_nrows_; }
    /// Return the number of columns in the block
    int b_ncols() const { return b_ncols_; }
    /// Returns the sub_block offset vector
    std::vector<int> block_sb_offset() const { return sb_offset_; }
    /// Returns the sub_block nrows vector;
    std::vector<int> block_sb_nrows() const { return sb_nrows_; }
    /// Returns the sub_block ncols vector;
    std::vector<int> block_sb_ncols() const { return sb_ncols_; }
    /// Returns the sub_block size vector;
    std::vector<int> block_sb_size() const { return sb_size_; }
    /// Returns the sub_block global_local_sblock map;
    std::map<int, int> global_local_sblock() const { return global_local_sblock_; }
    /// Returns the sub_block local_global_sblock vector;
    std::vector<int> local_global_sblock() const { return local_global_sblock_; }

    /// Returns the sub_block offset in the block
    int sb_offset(const int &sb) { return sb_offset_[global_local_sblock_[sb]]; }
    /// Return the number of rows in a given sub_block
    int sb_nrows(const int &sb) { return sb_nrows_[global_local_sblock_[sb]]; }
    /// Return the number of columns in a give sub_block
    int sb_ncols(const int &sb) { return sb_ncols_[global_local_sblock_[sb]]; }
    /// Return the size of a given sub_block
    int sb_size(const int &sb) { return sb_size_[global_local_sblock_[sb]]; }

#if HAVE_MADNESS == 1
    /// Initialize the block
    madness::Void initialize_block(const int &nrow, const int &ncol,
                                   const std::vector<int> &rsize,
                                   const std::vector<int> &csize,
                                   const std::vector<int> &roffset)
    {
        me_ = Communicator::world->me();
        nprocs_ = Communicator::world->nproc();
        nthreads_ = Communicator::world->nthread();
        comm_ = Communicator::world->communicator();
        if (block_.size()) clear_block();

        int tsb = nprocs_ * nprocs_;

        g_nrows_ = nrow;
        g_ncols_ = ncol;
        name_ = "Block " + to_string(me_);
        nsb_ = nprocs_;

        sb_nrows_ = rsize;
        for (int i=0; i < nsb_; i++) {
            sb_ncols_.push_back(csize[me_]);
        }

        b_nrows_ = nrow;
        b_ncols_ = sb_ncols_[0];

        nelements_ = b_nrows_*b_ncols_;

        sb_offset_.clear();
        for (int i=0; i < roffset.size(); i++) {
            sb_offset_.push_back(roffset[i]*sb_ncols_[i]);
        }


        for (int i=0; i < sb_offset_.size()-1; i++) {
            sb_size_.push_back(sb_offset_[i+1]-sb_offset_[i]);
        }
        sb_size_.push_back(nelements_ - sb_offset_[sb_offset_.size()-1]);

        block_.resize(nelements_, 0.0);

        for (int i=0, local=0; i < tsb; i++) {
            if (me_ == owner(i)) {
                global_local_sblock_.insert(std::pair<int,int>(i,local));
                local_global_sblock_.push_back(i);
                local++;
            }
        }

        return madness::None;
    }

    inline int owner(const int &sb) { return sb%nprocs_; }

    std::vector<double> get_block() { return block_; }

    std::vector<double> get_sblock(const int &sb) {
        int local = global_local_sblock_[sb];
        std::vector<double> tmp(sb_size_[local], 0.0);
        memcpy(&tmp[0], &block_[sb_offset_[local]], sb_size_[local]*sizeof(double));
        return tmp;
    }



    /// Set the appropriate sub_block of the full block to identity
    madness::Void set_block_identity()
    {
        memset(&block_[0], 0, block_.size()*sizeof(double));

        std::vector<double> iden(sb_ncols_[me_], 1.0);
        C_DCOPY(sb_ncols_[me_], &iden[0], 1, &block_[sb_offset_[me_]], sb_ncols_[me_]+1);

        return madness::None;
    }

    double val(const int &i, const int &j)
    {
        return block_[i*b_ncols_ + j];
    }

    void set(const int &i, const int &j, const double &val)
    {
        block_[i*b_ncols_ + j] = val;
    }

    /// set the block
    madness::Void copy_block(const std::vector<double> &copy)
    {
        free_vector(block_);
        block_ = copy;
        return madness::None;
    }
    /// set the sblock
    madness::Void copy_sblock(const int &sb, const std::vector<double> &copy)
    {
        int local = global_local_sblock_[sb];
        memcpy(&block_[local], &copy[0], copy.size()*sizeof(double));
        return madness::None;
    }


    madness::Void sum_block(const Block *blk)
    {
        for (int i=0; i < block_.size(); i++) {
            block_[i] += blk->block_[i];
        }
        return madness::None;
    }

    madness::Void sum_sblock(const int &sb, const std::vector<double> &blk)
    {
        int local = global_local_sblock_[sb];
        int offset = sb_offset_[local];
        int max = sb_size_[local] + offset;

        for (int i=offset; i < max; i++) {
            block_[i] += blk[i];
        }
        return madness::None;
    }

    madness::Void fill_block(const double &val)
    {
        for (int i=0; i < nelements_; i++) {
            block_[i] = val;
        }
        return madness::None;
    }

    madness::Void fill_sblock(const int &sb, const double &val)
    {
        int local = global_local_sblock_[sb];
        int offset = sb_offset_[local];
        int max = offset + sb_size_[local];

        for (int i=offset; i < max; i++) {
            block_[i] = val;
        }
        return madness::None;
    }

    madness::Void mxm(const int &sb, const std::vector<double> a,
                      const std::vector<double> &b,
                      const int &a_row, const int &a_col,
                      const int &b_row, const int &b_col)
    {
        int local = global_local_sblock_[sb];
        int offset = sb_offset_[local];

//        if (a_row == b_col) {
        for (int i=0, ij=0; i < a_row; i++) {
            for (int j=0; j < b_col; j++, ij++) {
                block_[offset + ij] += C_DDOT(a_col, const_cast<double*>(&a[i*a_col]), 1,
                                              const_cast<double*>(&b[j]), b_col);
            }
        }
        //        }
//        else throw PSIEXCEPTION("Rows and columns do not line up for matrix multiplication.\n");

        return madness::None;
    }

    madness::Void multiply_block(const int &sb,
                                 bool transa, bool transb, double alpha,
                                 const std::vector<double> &a,
                                 const std::vector<double> &b, double beta)
    {
        int offset = sb_offset_[global_local_sblock_[sb]];

        char ta = transa ? 't' : 'n';
        char tb = transb ? 't' : 'n';

        int m, n, k, lda, ldb, ldc;
        m = b_nrows_;
        n = b_ncols_;
        k = transa ? b_nrows_ : b_ncols_;
        lda = transa ? m : k;
        ldb = transb ? k : n;
        ldc = n;

        C_DGEMM(ta, tb, m, n, k, alpha, const_cast<double*>(&a[0]), lda,
                const_cast<double*>(&b[0]), ldb, beta, &block_[offset], ldc);

        return madness::None;
    }

#endif // End of HAVE_MADNESS

};


}


#endif // MATRIX_BLOCK_H
