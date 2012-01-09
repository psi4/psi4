/*
 *  dist_mat_set.cc
 *  part of distributed matrix
 *
 *  Created by Ben Mintz on 12/14/11.
 *
 */

#if HAVE_MADNESS

#include "dist_mat.h"
#include "libqt/qt.h"


namespace psi {

madness::Void Distributed_Matrix::copy_tile_tij(const int &tij, const madness::Tensor<double> &tile) //const double *tile)
{
    memcpy(data_[local(tij)].ptr(), const_cast<double*>(tile.ptr()), data_[local(tij)].size()*sizeof(double));
    return madness::None;
}

madness::Void Distributed_Matrix::copy_tile(const int &ti, const int &tj, const madness::Tensor<double> &tile)// const double *tile)//
{
    int tij = ti*tile_ncols_+tj;
    memcpy(data_[local(tij)].ptr(), const_cast<double*>(tile.ptr()), data_[local(tij)].size()*sizeof(double));
    return madness::None;
}


madness::Void Distributed_Matrix::copy_row(const int &yi, Distributed_Matrix &X, const int &xi)
{
    Communicator::world->sync();
    if (nelements_) {
        if (*this == X) {
            if (xi < X.nrows_) {
                if (yi < nrows_) {

                    int yti = this->convert_i_to_ti(yi);// yi/tile_sz_;
                    int xti = X.convert_i_to_ti(xi); // xi/X.tile_sz_;
                    for (int tj=0; tj < tile_ncols_; tj++) {
                        int y_tij = yti*tile_ncols_ + tj;
                        int x_tij = xti*X.tile_ncols_ + tj;
                        if (me_ == owner(y_tij)) {
                            bool local_copy = false;
                            if (me_ == X.owner(x_tij)) local_copy = true;

                            int ya = this->convert_i_to_a(yi);// yi%tile_sz_;
                            int xa = X.convert_i_to_a(xi);// xi%X.tile_sz_;


                            if (local_copy) {
                                this->task(me_, &Distributed_Matrix::local_copy_tile_row,
                                           ya, y_tij, const_cast<double*>(&X.data_[X.local(x_tij)](xa,0)));
                            }
                            else {
                                madness::Future<std::vector<double> > xrow = X.task(owner(x_tij), &Distributed_Matrix::get_tile_row, xa, x_tij);
                                this->task(me_, &Distributed_Matrix::non_local_copy_tile_row, xrow, ya, y_tij);
                            }

                        }
                    }

                }
                else throw PSIEXCEPTION("The row being copied into is out of bounds.\n");
            }
            else throw PSIEXCEPTION("The row being copied is out of bounds.\n");
        }
        else throw PSIEXCEPTION("The distributed matrices are not the same.\n");
    }
    Communicator::world->sync();
    return madness::None;
}


madness::Void Distributed_Matrix::local_copy_tile_row(const int &ya, const int &y_tij,
                                                      const double *x)
{
    C_DCOPY(data_[local(y_tij)].dim(1), const_cast<double*>(x), 1,
            &data_[local(y_tij)](ya,0), 1);
    return madness::None;
}

madness::Void Distributed_Matrix::non_local_copy_tile_row(const std::vector<double> &X,
                                                          const int &ya,
                                                          const int &y_tij)
{

//    this->data_[local_ytile](row,madness::_) = X(row, madness::_);

    C_DCOPY(data_[local(y_tij)].dim(1), const_cast<double*>(&X[0]), 1,
            &data_[local(y_tij)](ya,0), 1);
    return madness::None;
}



std::vector<double> Distributed_Matrix::get_tile_row(const int &a, const int &tij)
{
    int loc = local(tij);
    int ncols = data_[loc].dim(1);

    std::vector<double> temp(ncols);

    C_DCOPY(ncols, &data_[loc](a,0), 1, &temp[0], 1);

    return temp;
}

std::vector<double> Distributed_Matrix::get_tile_col(const int &b, const int &tij)
{
    int loc = local(tij);
    int nrows = data_[loc].dim(0);
    int ncols = data_[loc].dim(1);

    std::vector<double> temp(nrows);

    C_DCOPY(nrows, &data_[loc](0,b), ncols, &temp[0], 1);

    return temp;
}



madness::Void Distributed_Matrix::copy_col(const int &yj, Distributed_Matrix &X, const int &xj)
{
    Communicator::world->sync();
    if (nelements_) {
        if (*this == X) {
            if (xj < X.ncols_) {
                if (yj < ncols_) {

                    int y_tj = this->convert_j_to_tj(yj); // yj/tile_sz_;
                    int x_tj = X.convert_j_to_tj(xj); // xj/X.tile_sz_;
                    for (int ti=0; ti < tile_nrows_; ti++) {
                        int y_tij = ti*tile_ncols_ + y_tj;
                        int x_tij = ti*X.tile_ncols_ + x_tj;
                        if (me_ == owner(y_tij)) {
                            bool local_copy = false;
                            if (me_ == X.owner(x_tij)) local_copy = true;

                            int yb = this->convert_j_to_b(yj);// y_j%tile_sz_;
                            int xb = X.convert_j_to_b(xj);//  x_j%X.tile_sz_;


                            if (local_copy) {
                                this->task(me_, &Distributed_Matrix::local_copy_tile_col,
                                           yb, y_tij, const_cast<double*>(&X.data_[X.local(x_tij)](0,xb)),
                                           X.data_[X.local(x_tij)].dim(1));
                            }
                            else {
                                madness::Future<std::vector<double> > x_tile = X.task(owner(x_tij), &Distributed_Matrix::get_tile_col, xb, x_tij);
                                task(me_, &Distributed_Matrix::non_local_copy_tile_col, x_tile, y_tij, yb);
                            }
                        }
                    }

                }
                else throw PSIEXCEPTION("The column being copied into is out of bounds.\n");
            }
            else throw PSIEXCEPTION("The column being copied is out of bounds.\n");
        }
        else throw PSIEXCEPTION("The distributed matrices are not the same.\n");
    }
    Communicator::world->sync();
    return madness::None;
}


madness::Void Distributed_Matrix::local_copy_tile_col(const int &yb, const int &y_tij,
                                                      const double *X, const int &stride)
{
    C_DCOPY(data_[local(y_tij)].dim(0), const_cast<double*>(X), stride,
            &data_[local(y_tij)](0,yb), data_[local(y_tij)].dim(1));
    return madness::None;
}


madness::Void Distributed_Matrix::non_local_copy_tile_col(const std::vector<double> &X,
                                                          const int &y_tij,
                                                          const int &y_b)
{

    C_DCOPY(data_[local(y_tij)].dim(0), const_cast<double*>(&X[0]), 1,
            &data_[local(y_tij)](0,y_b), data_[local(y_tij)].dim(1));


    return madness::None;
}

madness::Void Distributed_Matrix::copy(const int &length, Distributed_Matrix &X,
                                       const int xi, const int xj, const int &x_inc,
                                       const int yi, const int yj, const int &y_inc)
{

//    int xi = x[0];
//    int xj = x[1];
//    int yi = y[0];
//    int yj = y[1];

    Communicator::world->sync();
    bool copy_all = false;
    bool same_dist = false;
    bool copy_all_row = false;
    bool copy_all_col = false;
    bool copy_tile_row = false;
    bool copy_tile_col = false;

    // check to see if the matrices are distributed the same
    if (*this == X) same_dist = true;
    if (same_dist) {
        // check to see if we are copying the entire distributed matrix
        if (length == nelements_)
            copy_all = true;
        // check to see if we are only copying a single row
        else if (length == ncols_ && xj == 0 && yj == 0 && x_inc == 1 && y_inc == 1)
            copy_all_row = true;
        // check to see if we are only copying a single col
        else if (length == nrows_ && xi == 0 && yi == 0 && x_inc == ncols_ && y_inc == ncols_)
            copy_all_col = true;
    }


    // Case: the matrices are distributed the same, copy the entire matrix, no communication
    if (copy_all)
        *this = X;
    // Case: Only copying a single row. There may be communication if xrow and yrow are not distributed the same.
    else if (copy_all_row)
        (*this)[yi]["*"] = X[xi]["*"];
    // Case: Only copying a single column.  There may be communication if xcol and ycol are not distributed the same.
    else if (copy_all_col)
        this->copy_col(yj, X, xj);
    // Case: this copies element by element.
    else {
        for (int i=0, y_index = yi * this->ncols_ + yj,
             x_index = xi * X.ncols_ + xj;
             i < length; i++, y_index += y_inc, x_index = x_inc) {

            int yrow = this->convert_ij_to_i(y_index);
            int ycol = this->convert_ij_to_j(y_index);
            int xrow = X.convert_ij_to_i(x_index);
            int xcol = X.convert_ij_to_j(x_index);

            if (this->me_ == this->owner(yrow/this->tile_sz_,ycol/this->tile_sz_))
                (*this)[yrow][ycol] = X[xrow][xcol].get_future();

        }
    }

    Communicator::world->sync();
    return madness::None;

}


Distributed_Matrix& Distributed_Matrix::transpose()
{

    Communicator::world->sync();
    Distributed_Matrix result(this->pgrid_, this->ncols_, this->nrows_,
                              this->tile_sz_, this->name_);

    for (int ti=0; ti < result.tile_nrows_; ti++) {
        for (int tj=0; tj < result.tile_ncols_; tj++) {
            if (me_ == result.owner(ti,tj) && me_ == result.owner(tj,ti)) {
                int tji = tj*this->ncols() + ti;

                result.task(me_, &Distributed_Matrix::local_copy_invert_tile, ti, tj, const_cast<double*>(this->data_[local(tji)].ptr()),
                            this->data_[local(tji)].dim(1));
            }
            else if (me_ == result.owner(ti,tj) && me_ != result.owner(tj,ti)) {
                madness::Future<madness::Tensor<double> > tile_ji = this->task(this->owner(tj,ti), &Distributed_Matrix::get_tile, tj, ti);
                result.task(me_, &Distributed_Matrix::non_local_copy_invert_tile, ti, tj, tile_ji);
            }
        }
    }

    Communicator::world->sync();
    std::swap(*this, result);
    return *this;
}

madness::Void Distributed_Matrix::local_copy_invert_tile(const int &ti, const int &tj,
                                                         const double *tptr,
                                                         const int &stride)
{
    int loc = local(ti,tj);
    int nrows = data_[loc].dim(0);
    double *ptr, *end;

    for (int t=0; t < nrows; t++) {
        ptr = &data_[loc](t,0);
        end = ptr + data_[loc].dim(1);
        while (ptr < end) {
            *ptr++ = *tptr;
            tptr+=stride;
        }
    }
}


madness::Void Distributed_Matrix::non_local_copy_invert_tile(const int &ti, const int &tj, const madness::Tensor<double> &tile)
{
    int loc = local(ti,tj);
    int stride = tile.dim(1);
    int nrows = data_[loc].dim(0);
    double *ptr, *end, *tptr;

    for (int t=0; t < nrows; t++) {
        ptr = &data_[loc](t,0);
        end = ptr + data_[loc].dim(1);
        tptr = const_cast<double*>(&tile(0,t));
        while (ptr < end) {
            *ptr++ = *tptr;
            tptr+=stride;
        }
    }
}

} // End of namespace psi

#endif
