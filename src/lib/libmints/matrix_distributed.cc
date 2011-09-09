/*
 *  matrix_distributed.cc
 *  distributed_matrix
 *
 *  Created by Ben Mintz on 8/24/11.
 *  Copyright 2011 __MintzInc__. All rights reserved.
 *
 */

#include <exception.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <libmints/integral.h>
#include <libdpd/dpd.h>
#include "factory.h"
#include "wavefunction.h"
#include "dimension.h"
#include "matrix_distributed.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <ctype.h>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>

using namespace boost;
using namespace psi;
using namespace std;

// In molecule.cc
namespace psi {
extern int str_to_int(const std::string& s);
extern double str_to_double(const std::string& s);
}

#if HAVE_MADNESS == 1

Distributed_Matrix::Distributed_Matrix()
    : Block(), nblocks_(0), nrows_(0), ncols_(0), tsb_(0), telements_(0)
    , madworld_( Communicator::world->get_madworld() ),
      madness::WorldObject<Distributed_Matrix>(*Communicator::world->get_madworld())
{
    me_ = Communicator::world->me();
    nprocs_ = Communicator::world->nproc();
    nthreads_ = Communicator::world->nthread();
    comm_ = Communicator::world->communicator();
    print_mutex_ = Communicator::world->get_mutex();
    block_roffset_.clear();
    block_coffset_.clear();
    block_csize_.clear();
    block_rsize_.clear();

    process_pending();
}


Distributed_Matrix::Distributed_Matrix(const int &rows, const int &cols)
    : Block()
    , madworld_( Communicator::world->get_madworld() ),
      madness::WorldObject<Distributed_Matrix>(*Communicator::world->get_madworld())
{
    me_ = Communicator::world->me();
    nprocs_ = Communicator::world->nproc();
    nthreads_ = Communicator::world->nthread();
    comm_ = Communicator::world->communicator();
    print_mutex_ = Communicator::world->get_mutex();

    common_init(rows, cols);

    process_pending();
}

Distributed_Matrix::Distributed_Matrix(const std::string &name, const int &rows, const int &cols)
    : Block()
    , madworld_( Communicator::world->get_madworld() ),
      madness::WorldObject<Distributed_Matrix>(*Communicator::world->get_madworld())
{
    me_ = Communicator::world->me();
    nprocs_ = Communicator::world->nproc();
    nthreads_ = Communicator::world->nthread();
    comm_ = Communicator::world->communicator();
    print_mutex_ = Communicator::world->get_mutex();
    name_ = name;
    common_init(rows, cols);

    process_pending();
}

madness::Void Distributed_Matrix::clear()
{
    nblocks_ = 0;
    nrows_ = 0;
    ncols_ = 0;
    tsb_ = 0;
    telements_ = 0;
    block_roffset_.clear();
    block_coffset_.clear();
    block_csize_.clear();
    block_rsize_.clear();
    clear_block();
    return madness::None;
}

madness::Void Distributed_Matrix::common_init(const int &rows, const int &cols)
{    
    nblocks_ = nprocs_;
    nrows_ = rows;
    ncols_ = cols;
    tsb_ = nprocs_*nprocs_;
    telements_ = nrows_*ncols_;

    block_roffset_.clear();
    block_coffset_.clear();
    block_csize_.clear();
    block_rsize_.clear();

    // Everyone will get at least this many rows and columns
    for (int i=0; i < nblocks_; i++) {
        block_rsize_.push_back(nrows_/nblocks_);
        block_csize_.push_back(ncols_/nblocks_);
    }
    // This adds in the left over rows (i.e. rows%nblocks)
    for (int i=0; i < nrows_%nblocks_; i++) {
        block_rsize_[i]++;
    }
    // This adds in the left over columns (i.e. cols%nblocks)
    for (int i=0; i < ncols_%nblocks_; i++) {
        block_csize_[i]++;
    }

    // Here we determine the offsets for the rows and columns
    block_roffset_.push_back(0);
    block_coffset_.push_back(0);
    for (int i=1; i < nprocs_; i++) {
        block_roffset_.push_back(block_roffset_[i-1] + block_rsize_[i-1]);
        block_coffset_.push_back(block_coffset_[i-1] + block_csize_[i-1]);
    }

    initialize_block(nrows_, ncols_, block_rsize_,
                     block_csize_, block_roffset_);

    return madness::None;
}


void Distributed_Matrix::print_all_blocks() const
{
    if (me_ == 0) {
        if (nelements() != 0) {
            for (int i=0; i < nprocs_; i++)
                task(me_, &Distributed_Matrix::print_block, i);
        }
        else {
            if (name_.size() != 0)
                fprintf(outfile, "  ## %s (is empty) ##\n", name_.c_str());
            else
                fprintf(outfile, "  ## Distributed Matrix (is empty) ##\n");
        }
    }
    Communicator::world->sync();
}

madness::Void Distributed_Matrix::print_block(const int &i) const
{
    if (me_ == 0) {
        // Get the block and info from the owner
        int own = owner(i);
        madness::Future<std::vector<double> > block = task(own, &Distributed_Matrix::get_block);
        madness::Future<int> rows = task(own, &Distributed_Matrix::b_nrows);
        madness::Future<int> cols = task(own, &Distributed_Matrix::b_ncols);

        std::string fname;
        if (name_.size() != 0) fname = name_ + ": Block " + to_string(i);
        else fname = "Block " + to_string(i);

        // Print the block
        task(me_, &Distributed_Matrix::print_mat, block, rows, cols, fname);
    }

    return madness::None;
}

void Distributed_Matrix::print_all_sblocks() const
{
    if (me_ == 0) {
        if (nelements() != 0) {
            for (int i=0; i < tsb_; i++)
                task(me_, &Distributed_Matrix::print_sblock, i);
        }
        else {
            if (name_.size() != 0)
                fprintf(outfile, "  ## %s (is empty) ##\n", name_.c_str());
            else
                fprintf(outfile, "  ## Distributed Matrix (is empty) ##\n");
        }
    }
    Communicator::world->sync();
}

madness::Void Distributed_Matrix::print_sblock(const int &i) const
{
    if (me_ == 0) {
        // Get the block and info from the owner
        int own = owner(i);
        madness::Future<std::vector<double> > sblock = task(own, &Distributed_Matrix::get_sblock, i);
        madness::Future<int> rows = task(own, &Distributed_Matrix::sb_nrows, i);
        madness::Future<int> cols = task(own, &Distributed_Matrix::sb_ncols, i);

        std::string fname;
        if (name_.size() != 0) fname = name_ + ": Block " + to_string(own);
        else fname = "Block " + to_string(own);

        fname += ": Sub_Block " + to_string(i);

        // Print the block
        task(me_, &Distributed_Matrix::print_mat, sblock, rows, cols, fname);
    }

    return madness::None;
}

madness::Void Distributed_Matrix::print_mat(const std::vector<double> &a, const int &m,
                                            const int &n, const std::string &name) const
{
    print_mutex_->lock();

    fprintf(outfile, "\n  ## %s ##\n", name.c_str());

    if (a.size() == 0) {
        fprintf(outfile, "\n\t## %s ## (empty)\n", name.c_str());
    }
    else {
        // Should only be called by process 0
        const int print_ncol = 5;
        int num_frames = int(n/print_ncol);
        int num_frames_rem = n%print_ncol; //adding one for changing 0->1 start
        int num_frame_counter = 0;
        //for each frame
        for(num_frame_counter=0;num_frame_counter<num_frames;num_frame_counter++){
            fprintf(outfile,"\n");
            for(int j=print_ncol*num_frame_counter+1;j<print_ncol*num_frame_counter+print_ncol+1;j++){
                if(j==print_ncol*num_frame_counter+1){ fprintf(outfile,"%18d",j); }
                else{ fprintf(outfile,"        %5d",j); }
            }
            fprintf(outfile,"\n\n");

            for(int k=1; k<=m; ++k){
                for(int j=print_ncol*num_frame_counter+1;j<print_ncol*num_frame_counter+print_ncol+2;j++){
                    if(j==print_ncol*num_frame_counter+1){ fprintf(outfile,"%5d",k);}
                    else{ fprintf(outfile," %12.7f",a[(k-1)*n + (j-2)]); }
                }
                fprintf(outfile,"\n");
            }
        }

        // ALREADY DID THE FULL FRAMES BY THIS POINT
        // NEED TO TAKE CARE OF THE REMAINDER
        if(num_frames_rem != 0){
            fprintf(outfile,"\n");
            for(int j=print_ncol*num_frame_counter+1;j<=n;j++){
                if(j==print_ncol*num_frame_counter+1){ fprintf(outfile,"%18d",j); }
                else{ fprintf(outfile,"        %5d",j); }
            }
            fprintf(outfile,"\n\n");

            for(int k=1; k<=m; ++k){
                for(int j=print_ncol*num_frame_counter+1;j<n+2;j++){
                    if(j==print_ncol*num_frame_counter+1){ fprintf(outfile,"%5d",k); }
                    else{ fprintf(outfile," %12.7f",a[(k-1)*n + (j-2)]); }
                }
                fprintf(outfile,"\n");
            }
        }
        fprintf(outfile,"\n\n");
    }
    print_mutex_->unlock();

    return madness::None;
}

madness::Void Distributed_Matrix::identity()
{
    if (nelements() != 0) {
        if (nrows_ == ncols_) {
            task(me_, &Distributed_Matrix::set_block_identity);
        }
        else throw PSIEXCEPTION("The matrix is not square.\n");
    }
    else {
        throw PSIEXCEPTION("You can not set a empty distributed matrix to the identity.\n");
    }
    Communicator::world->sync();
    return madness::None;
}

madness::Future<double> Distributed_Matrix::get_val(const int &row, const int &col)
{
    if (nelements() != 0) {
        int own = nprocs_ - 1;
        for (int i=0; i < nprocs_-1; i++) {
            if (col < block_coffset_[i+1]) {
                own = i;
                break;
            }
        }

        if (me_ == own) {
            int local_column = col - block_coffset_[own];
            return madness::Future<double> ( val(row, local_column) );
        }
        else
            return task(own, &Distributed_Matrix::get_val, row, col);
    }
    else {
        throw PSIEXCEPTION("There are no values to get in an empty distributed matrix.\n");
    }
}

madness::Void Distributed_Matrix::set_val(const int &row, const int &col, const double &val)
{
    // There is no sync here.
    // Be careful here because this will return
    // immediately after the tasks have been submitted (not run)
    if (nelements() != 0) {
        int own = nprocs_ - 1;
        for (int i=0; i < nprocs_-1; i++) {
            if (col < block_coffset_[i+1]) {
                own = i;
                break;
            }
        }

        if (me_ == own) {
            int local_column = col - block_coffset_[own];

            set(row, local_column, val);
            return madness::None;
        }
        else {
            task(own, &Distributed_Matrix::set_val, row, col, val);
            return madness::None;
        }
    }
    else {
        throw PSIEXCEPTION("There are no values to get in an empty distributed matrix.\n");
    }
}

Distributed_Matrix& Distributed_Matrix::operator= (const Distributed_Matrix &rhs)
{
    this->clear();
    int rows = rhs.nrows();
    int cols = rhs.ncols();
    if (name_.size() == 0) name_ = rhs.name();
    common_init(rows, cols);

    for (int i=0; i < nblocks_; i++) {
        int own = owner(i);
        if (me_ == own) {
            madness::Future<std::vector<double> > blk = rhs.task(own, &Distributed_Matrix::get_block);
            this->task(own, &Distributed_Matrix::copy_block, blk);
        }
    }

    // We need a sync here to make sure that all of the copying is done before moving on
    Communicator::world->sync();
    return *this;
}

Distributed_Matrix& Distributed_Matrix::operator =(const Distributed_Matrix *rhs)
{
    this->clear();
    int rows = rhs->nrows();
    int cols = rhs->ncols();
    if (name_.size() == 0) name_ = rhs->name();
    common_init(rows, cols);

    for (int i=0; i < nblocks_; i++) {
        int own = owner(i);
        if (me_ == own) {
            madness::Future<std::vector<double> > blk = rhs->task(own, &Distributed_Matrix::get_block);
            this->task(own, &Distributed_Matrix::copy_block, blk);
        }
    }

    // We need a sync here to make sure that all of the copying is done before moving on
    Communicator::world->sync();
    return *this;
}

Distributed_Matrix::Distributed_Matrix(const Distributed_Matrix &copy)
    : madworld_( Communicator::world->get_madworld() ),
      madness::WorldObject<Distributed_Matrix>(*Communicator::world->get_madworld())
{
    me_ = Communicator::world->me();
    nprocs_ = Communicator::world->nproc();
    nthreads_ = Communicator::world->nthread();
    comm_ = Communicator::world->communicator();
    print_mutex_ = Communicator::world->get_mutex();

    *this = copy;
}

madness::Void Distributed_Matrix::operator +=(const Distributed_Matrix &rhs)
{
    if (*this == rhs) {
        madness::Future<std::vector<double> > block = rhs.task(me_, &Distributed_Matrix::get_block);
        for (int sb=0; sb < tsb_; sb++) {
            if (me_ == owner(sb))
                this->task(me_, &Distributed_Matrix::sum_sblock, sb, block);
        }
    }
    else throw PSIEXCEPTION("The matrices being added are not the same.\n");

    Communicator::world->sync();
    return madness::None;
}

madness::Void Distributed_Matrix::fill(const double &val)
{
    if (this->nelements() != 0) {
        for (int sb=0; sb < tsb_; sb++) {
            if (me_ == owner(sb))
                this->task(me_, &Distributed_Matrix::fill_sblock, sb, val);
        }
    }
    else throw PSIEXCEPTION("The matrix you are tryin to fill is empty.\n");

    Communicator::world->sync();
    return madness::None;
}

madness::Void Distributed_Matrix::operator =(const double &val)
{
    fill(val);
    return madness::None;
}


Distributed_Matrix Distributed_Matrix::operator +(const Distributed_Matrix &rhs)
{
    if (*this == rhs) {
        Distributed_Matrix tmp = *this;
        tmp += rhs;
        // The += operator has a sync, so we don't need one here
        return tmp;
    }
    else throw PSIEXCEPTION("The matrices that are being added are not the same.\n");
}

bool Distributed_Matrix::operator ==(const Distributed_Matrix &rhs)
{
    if ( this->g_nrows_ != rhs.g_nrows() ) return false;
    else if ( this->g_ncols_ != rhs.g_ncols() ) return false;
    else if ( this->nelements_ != rhs.nelements() ) return false;
    else if ( this->nsb_ != rhs.nsb() ) return false;
    else if ( this->sb_offset_ != rhs.block_sb_offset() ) return false;
    else if ( this->sb_nrows_ != rhs.block_sb_nrows() ) return false;
    else if ( this->sb_ncols_ != rhs.block_sb_ncols() ) return false;
    else if ( this->sb_size_ != rhs.block_sb_size() ) return false;
    else return true;
}


bool Distributed_Matrix::operator !=(const Distributed_Matrix &rhs)
{
    if (*this == rhs) return false;
    else return true;
}

Distributed_Matrix Distributed_Matrix::operator* (const Distributed_Matrix &rhs)
{
    bool symmetric = true;

    // I have assumed a symmetric matrice to start with
    if (symmetric) {
        Distributed_Matrix result("multiply", this->nrows(), this->ncols());


//        for (int j=0; j < nblocks_; j++) {
//            if (me_ == owner(j)) {
//                madness::Future<std::vector<double> > b_block = rhs.task(owner(j), &Distributed_Matrix::get_block);
//                for (int i=0; i <= j; i++) {
//                    int sb = i*nblocks_ + j;
//                    madness::Future<std::vector<double> > a_block = this->task(owner(i), &Distributed_Matrix::get_block);

//                    result.task(owner(sb), &Distributed_Matrix::multiply_block,
//                                sb, true, false, 1.0, a_block, b_block, 0.0);
//                }
//            }
//        }

//        Communicator::world->sync();

//        this->print_all_blocks();
//        rhs.print_all_blocks();
//        result.print_all_blocks();

//        Communicator::world->sync();

//        for (int i=0; i < nblocks_; i++) {
//            for (int j=0; j > i; j++) {
//                int ij = i*ncols_ + j;
//                int ji = j*ncols_ + i;
//                if ( owner(ji) )
//                    result.task(owner(ji), &Distributed_Matrix::copy_nonlocal_sblock, ij, ji);
//            }
//        }
//        Communicator::world->sync();



        for (int i=0; i < nblocks_; i++) {
            for (int k=0; k < nblocks_; k++) {
                int ik = i*nblocks_ + k;
                if (me_ == owner(ik)) {
                    for (int j=0; j < nblocks_; j++) {
                        int ij = i*nblocks_ + j;
                        int jk = j*nblocks_ + k;
                        madness::Future<std::vector<double> > A = this->task(owner(ij), &Distributed_Matrix::get_sblock, ij);
                        madness::Future<std::vector<double> > B = rhs.task(owner(jk), &Distributed_Matrix::get_sblock, jk);
                        madness::Future<int> a_row = this->task(owner(ij), &Distributed_Matrix::sb_nrows, ij);
                        madness::Future<int> a_col = this->task(owner(ij), &Distributed_Matrix::sb_ncols, ij);
                        madness::Future<int> b_row = rhs.task(owner(jk), &Distributed_Matrix::sb_nrows, jk);
                        madness::Future<int> b_col = rhs.task(owner(jk), &Distributed_Matrix::sb_ncols, jk);

                        result.task(owner(ik), &Distributed_Matrix::mxm, ik, A, B, a_row, a_col, b_row, b_col);
                    }
                }
            }
        }


        Communicator::world->sync();

        this->print_all_blocks();
        rhs.print_all_blocks();
        result.print_all_blocks();


        return result;
    }
    else {
        throw PSIEXCEPTION("Only symmetric matrix multiplies currently work.\n");
    }
}

#endif // End of HAVE_MADNESS













