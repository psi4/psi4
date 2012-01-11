/*
 *  dist_mat_set.cc
 *  part of distributed matrix
 *
 *  Created by Ben Mintz on 12/14/11.
 *
 */


#include "dist_mat.h"

#ifdef HAVE_MADNESS

namespace psi {

    // default constructor
    Distributed_Matrix::Distributed_Matrix() :
        madness::WorldObject<Distributed_Matrix>(*Communicator::world->get_madworld())
    {
        parallel_init();
        process_pending();
    }

    // constructor
    Distributed_Matrix::Distributed_Matrix(const process_grid<2> &pgrid,
                                           const int &nrows, const int &ncols,
                                           const int &tile_sz,
                                           const std::string &name) :
        madness::WorldObject<Distributed_Matrix>(*Communicator::world->get_madworld())
    {
        parallel_init();
        initialize(pgrid, nrows, ncols, tile_sz, name);
        process_pending();
    }

    // initialize parallel stuff
    void Distributed_Matrix::parallel_init()
    {
        me_ = Communicator::world->me();
        nprocs_ = Communicator::world->nproc();
        nthreads_ = Communicator::world->nthread();
        comm_ = Communicator::world->communicator();
        print_mutex_ = Communicator::world->get_mutex();
        add_mutex_ = Communicator::world->get_mutex();
        mult_mutex_ = Communicator::world->get_mutex();
        set_mutex_ = Communicator::world->get_mutex();
        madworld_ = Communicator::world->get_madworld();
    }

    // initialize the distributed matrix
    void Distributed_Matrix::initialize(const process_grid<2> &pgrid,
                                         const int &nrows, const int &ncols,
                                         const int &tile_sz,
                                         const std::string &name)
    {
        clear_matrix();

        pgrid_ = pgrid;
        name_ = name;

        for (int i=0; i < pgrid_.ndims(); i++) {
            pgrid_dims_.push_back(pgrid_.dim_size(i));
        }

        // initialize the tile size
        tile_sz_ = tile_sz;

        // initialize the global number of rows and columns in the matrix
        nrows_ = nrows;
        ncols_ = ncols;
        // determine the number of elements in the matrix
        nelements_ = nrows_*ncols_;

        // Compute the total number of tile rows and columns
        tile_nrows_ = nrows_/tile_sz_;
        tile_ncols_ = ncols_/tile_sz_;
        if (nrows_%tile_sz_) tile_nrows_++;
        if (ncols_%tile_sz_) tile_ncols_++;

        // compute the total number of tiles in the matrix
        ntiles_ = tile_nrows_*tile_ncols_;

        // set up the (global to local) and (local to global) maps and allocate the tiles
        local_ntiles_ = 0;
        for (int tij=0; tij < ntiles_; tij++) {
            if (me_ == owner(tij)) {
                global_to_local_map_.insert(std::pair<int,int>(tij,local_ntiles_));
                local_to_global_map_.insert(std::pair<int,int>(local_ntiles_,tij));
                local_ntiles_++;
            }
        }

        // allocate the tiles
        for (int tij=0; tij < local_ntiles_; tij++)
            data_.push_back(madness::Tensor<double>());

//        Communicator::world->sync();
//        Communicator::world->sync();

        for (int ti=0, tij=0; ti < tile_nrows_; ti++) {
            for (int tj=0; tj < tile_ncols_; tj++, tij++) {
                if (me_ == owner(tij)) {
                    if (nrows_%tile_sz_ && ncols_%tile_sz_) {
                        if ((ti+1)%tile_nrows_ && (tj+1)%tile_ncols_)
                            data_[local(tij)] = madness::Tensor<double>(tile_sz_, tile_sz_);
                        if ((ti+1)%tile_nrows_ == 0 && (tj+1)%tile_ncols_)
                            data_[local(tij)] = madness::Tensor<double>(nrows_%tile_sz_, tile_sz_);
                        if ((ti+1)%tile_nrows_ && (tj+1)%tile_ncols_ == 0)
                            data_[local(tij)] = madness::Tensor<double>(tile_sz_, ncols_%tile_sz_);
                        else if ((ti+1)%tile_nrows_ == 0 && (tj+1)%tile_ncols_ == 0)
                            data_[local(tij)] = madness::Tensor<double>(nrows_%tile_sz_, ncols_%tile_sz_);
                    }
                    else if (nrows_%tile_sz_ == 0 && ncols_%tile_sz_) {
                        if ((tj+1)%tile_ncols_)
                            data_[local(tij)] = madness::Tensor<double>(tile_sz_, tile_sz_);
                        else
                            data_[local(tij)] = madness::Tensor<double>(tile_sz_, ncols_%tile_sz_);
                    }
                    else if (nrows_%tile_sz_ && ncols_%tile_sz_ == 0) {
                        if ((ti+1)%tile_nrows_)
                            data_[local(tij)] = madness::Tensor<double>(tile_sz_, tile_sz_);
                        else
                            data_[local(tij)] = madness::Tensor<double>(nrows_%tile_sz_, tile_sz_);
                    }
                    else
                        data_[local(tij)] = madness::Tensor<double>(tile_sz_, tile_sz_);
                }
            }
        }


//        print_matrix_info();
    }

    // clear the distributed matrix
    void Distributed_Matrix::clear_matrix()
    {

        nrows_ = 0;
        ncols_ = 0;
        nelements_ = 0;

        nrows_ = 0;
        ncols_ = 0;
        nelements_ = 0;

        tile_sz_ = 0;
        tile_nrows_ = 0;
        tile_ncols_ = 0;
        ntiles_ = 0;
        local_ntiles_ = 0;

        free_map(global_to_local_map_);
        free_map(local_to_global_map_);
        free_vector(pgrid_dims_);
        free_vector(data_);

//        free_vector(data_);
    }

    Distributed_Matrix::~Distributed_Matrix()
    {
        data_.clear();
        pgrid_dims_.clear();
        global_to_local_map_.clear();
        local_to_global_map_.clear();
//        free_vector(data_);
//        if (data_ != NULL)
//            delete(data_);
    }

    Distributed_Matrix& Distributed_Matrix::operator= (const Distributed_Matrix &copy)
    {
        this->clear_matrix();
        Communicator::world->sync();
        if (this->name_.size()) initialize(copy.pgrid_, copy.nrows_, copy.ncols_,
                                            copy.tile_sz_, this->name_);
        else initialize(copy.pgrid_, copy.nrows_, copy.ncols_,
                         copy.tile_sz_, copy.name_);

        for (int tij=0; tij < this->ntiles_; tij++) {
            if (me_ == owner(tij)) {
                madness::Future<madness::Tensor<double> > tile = copy.task(me_, &Distributed_Matrix::get_tile_tij, tij);
                task(me_, &Distributed_Matrix::copy_tile_tij, tij, tile);
            }
        }
        Communicator::world->sync();
        return *this;
    }

    Distributed_Matrix& Distributed_Matrix::operator= (const Distributed_Matrix *copy)
    {
        this->clear_matrix();
        Communicator::world->sync();
        if (this->name_.size()) initialize(copy->pgrid_, copy->nrows_, copy->ncols_,
                                            copy->tile_sz_, this->name_);
        else initialize(copy->pgrid_, copy->nrows_, copy->ncols_,
                         copy->tile_sz_, copy->name_);

        for (int tij=0; tij < this->ntiles_; tij++) {
            if (me_ == owner(tij)) {
                madness::Future<madness::Tensor<double> > tile = copy->task(me_, &Distributed_Matrix::get_tile_tij, tij);
                task(me_, &Distributed_Matrix::copy_tile_tij, tij, tile);
            }
        }
        Communicator::world->sync();
        return *this;
    }

    Distributed_Matrix::Distributed_Matrix(const Distributed_Matrix &copy)
        : madness::WorldObject<Distributed_Matrix>(*Communicator::world->get_madworld()),
          data_(NULL)
    {
        Communicator::world->sync();
        parallel_init();
        process_pending();
        *this = copy;
        Communicator::world->sync();
    }
    Distributed_Matrix::Distributed_Matrix(const Distributed_Matrix *copy)
        : madness::WorldObject<Distributed_Matrix>(*Communicator::world->get_madworld()),
          data_(NULL)
    {
        Communicator::world->sync();
        parallel_init();
        process_pending();
        *this = copy;
        Communicator::world->sync();
    }





} // End of namespace psi


#endif
