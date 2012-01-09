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

using namespace psi;
using namespace std;

namespace psi {

madness::Void Distributed_Matrix::zero()
{
    Communicator::world->sync();
    for (int tij=0; tij < ntiles_; tij++) {
        if (me_ == owner(tij)) {
            task(me_, &Distributed_Matrix::zero_tile_tij, tij);
        }
    }
    Communicator::world->sync();
    return madness::None;
}

madness::Void Distributed_Matrix::zero_tile(const int &ti, const int &tj)
{
    int loc = local(ti,tj);
    data_[loc].fill(0.0);
    return madness::None;
}

madness::Void Distributed_Matrix::zero_tile_tij(const int &tij)
{
    int loc = local(tij);
    data_[loc].fill(0.0);
    return madness::None;
}

madness::Void Distributed_Matrix::identity()
{
    Communicator::world->sync();
    zero();
    for (int ti=0; ti < tile_nrows_; ti++) {
        if (me_ == owner(ti,ti)) {
            task(me_, &Distributed_Matrix::set_tile_to_identity, ti, ti);
        }
    }
    Communicator::world->sync();
    return madness::None;
}

madness::Void Distributed_Matrix::set_tile_to_identity(const int &ti, const int &tj)
{
    set_tile_to_identity_tij(ti*tile_ncols_+tj);
    return madness::None;
}

madness::Void Distributed_Matrix::set_tile_to_identity_tij(const int &tij)
{

    int loc = local(tij);
    int stride = data_[loc].dim(1) + 1;
    double* ptr = data_[loc].ptr();
    double* end = ptr + data_[loc].size();


    while (ptr <= end) {
        *ptr = 1.0;
        ptr+=stride;
    }

    return madness::None;
}


madness::Void Distributed_Matrix::zero_diagonal()
{
    Communicator::world->sync();
    for (int ti=0; ti < tile_nrows_; ti++) {
        if (me_ == owner(ti,ti)) {
            task(me_, &Distributed_Matrix::set_tile_diagonal_zero, ti, ti);
        }
    }
    Communicator::world->sync();
    return madness::None;
}

madness::Void Distributed_Matrix::set_tile_diagonal_zero(const int &ti, const int &tj)
{
    set_tile_diagonal_zero_tij(ti*tile_ncols_+tj);
    return madness::None;
}

madness::Void Distributed_Matrix::set_tile_diagonal_zero_tij(const int &tij)
{
    int loc = local(tij);
    int stride = data_[loc].dim(1) + 1;
    double* ptr = data_[loc].ptr();
    double* end = ptr + data_[loc].size();


    while (ptr <= end) {
        *ptr = 0.0;
        ptr+=stride;
    }

    return madness::None;
}


madness::Future<double> Distributed_Matrix::get_val(const int &i, const int &j)
{
    if (nelements_) {
        if (i < nrows_ && j < ncols_) {
            int ti = convert_i_to_ti(i);// i/tile_sz_;
            int tj = convert_j_to_tj(j); // j/tile_sz_;
            int tij = ti*tile_ncols_ + tj;
            if (me_ == owner(tij)) {
                int a = convert_i_to_a(i); // i%tile_sz_;
                int b = convert_j_to_b(j); // j%tile_sz_;
                return task(me_, &Distributed_Matrix::return_tile_val, tij, a, b);
            }
            else return task(owner(tij), &Distributed_Matrix::get_val, i, j);
        }
        else
            throw PSIEXCEPTION("The value is out of bounds.\n");

    }
    else throw PSIEXCEPTION("There are no values to get in an empty distributed matrix.\n");
}

madness::Future<double> Distributed_Matrix::return_tile_val(const int &tij, const int &a,
                                                            const int &b)
{
    return madness::Future<double> ( data_[local(tij)](a,b) );
}


madness::Void Distributed_Matrix::set_val(const int &i, const int &j,
                                          const double &val)
{
    if (nelements_) {
        if (i < nrows_ && j < ncols_) {
            int ti = convert_i_to_ti(i);
            int tj = convert_j_to_tj(j);
            int tij = ti*tile_ncols_ + tj;
            if (me_ == owner(tij)) {
                int a = convert_i_to_a(i);
                int b = convert_j_to_b(j);
                task(me_, &Distributed_Matrix::set_tile_val, tij, a, b, val);
            }
            else task(owner(tij), &Distributed_Matrix::set_val, i, j, val);
        }
        else
            throw PSIEXCEPTION("The value is out of bounds.\n");

    }
    else throw PSIEXCEPTION("There are no values to get in an empty distributed matrix.\n");

    return madness::None;
}

madness::Void Distributed_Matrix::set_tile_val(const int &tij, const int &a,
                                               const int &b, const double &val)
{
    set_mutex_->lock();
    data_[local(tij)](a,b) = val;
    set_mutex_->unlock();
    return madness::None;
}

madness::Void Distributed_Matrix::set_row(const int &i, const double & val)
{
    Communicator::world->sync();
    if (nelements_) {
        if (i < nrows_) {
            int ti = convert_i_to_ti(i);// i/tile_sz_;
            for (int tj=0; tj < tile_ncols_; tj++) {
                int tij = ti*tile_ncols_ + tj;
                if (me_ == owner(tij)) {
                    int a = convert_i_to_a(i);// i%tile_sz_;
                    task(me_, &Distributed_Matrix::set_tile_row, tij, a, val);
                }
            }
        }
    }
    Communicator::world->sync();
    return madness::None;
}

madness::Void Distributed_Matrix::set_tile_row(const int &tij,
                                               const int &a,
                                               const double &val)
{
    int loc = local(tij);
    data_[loc](a,madness::_) = val;

    //    double *ptr = &data_[loc](loc_row,0); //[loc_row].begin();
    //    double *end = ptr + data_[loc].size();

    //    while (ptr != end)
    //        *ptr++ = val;

    return madness::None;
}

madness::Void Distributed_Matrix::set_col(const int &j, const double & val)
{
    Communicator::world->sync();
    if (nelements_) {
        if (j < ncols_) {
            int tj = convert_j_to_tj(j);  // j/tile_sz_;
            for (int ti=0; ti < tile_nrows_; ti++) {
                int tij = ti*tile_ncols_ + tj;
                if (me_ == owner(tij)) {
                    int b = convert_j_to_b(j);  // j%tile_sz_;
                    task(me_, &Distributed_Matrix::set_tile_col, tij, b, val);
                }
            }
        }
    }
    Communicator::world->sync();
    return madness::None;
}


madness::Void Distributed_Matrix::set(const int &i, const int &j, const double &val)
{
    this->set_val(i,j,val);
}

madness::Void Distributed_Matrix::set(const int &i, const std::string &str, const double &val)
{
    if (str == "*") set_row(i, val);
    else throw PSIEXCEPTION("To specify an entire row use set_row(i,*).\n");
}

madness::Void Distributed_Matrix::set(const std::string &str, const int &j, const double &val)
{
    if (str == "*") set_col(j, val);
    else throw PSIEXCEPTION("To specify an entire row use set_row(i,*).\n");
}

madness::Void Distributed_Matrix::set(const int &i, const int &j, madness::Future<double> val)
{
    this->set_val(i,j,val);
}

madness::Void Distributed_Matrix::set(const int &i, const std::string &str, madness::Future<double> val)
{
    if (str == "*") set_row(i, val);
    else throw PSIEXCEPTION("To specify an entire row use set_row(i,*).\n");
}

madness::Void Distributed_Matrix::set(const std::string &str, const int &j, madness::Future<double> val)
{
    if (str == "*") set_col(j, val);
    else throw PSIEXCEPTION("To specify an entire row use set_row(i,*).\n");
}



madness::Void Distributed_Matrix::set_tile_col(const int &tij,
                                               const int &b,
                                               const double &val)
{
    int loc = local(tij);
    //    data_[loc](0,col) = val;
    int stride = data_[loc].dim(1);
    double *ptr = &data_[loc](0,b);
    double *end = ptr + data_[loc].size();

    while (ptr < end) {
        *ptr = val;
        ptr+=stride;
    }

    return madness::None;
}


madness::Void Distributed_Matrix::fill(const double &val)
{
    Communicator::world->sync();
    if (nrows_*ncols_) {
        for (int ti=0; ti < tile_nrows_; ti++) {
            for (int tj=0; tj < tile_ncols_; tj++) {
                if (me_ == owner(ti,tj)) {
                    task(me_, &Distributed_Matrix::fill_tile,
                         ti*tile_ncols_+tj, val);
                }
            }
        }
    }
    else throw PSIEXCEPTION("The matrix you are trying to fill is empty.\n");
    Communicator::world->sync();
    return madness::None;
}

madness::Void Distributed_Matrix::fill_tile(const int &tij, const double &val)
{

    double* ptr = data_[local(tij)].ptr();
    double* end = ptr + data_[local(tij)].size();

    while (ptr != end)
        *ptr++ = val;

    return madness::None;
}

madness::Void Distributed_Matrix::operator +=(const Distributed_Matrix &rhs)
{
    if (*this == rhs) {
        Communicator::world->sync();
        for (int tij=0; tij < ntiles_; tij++) {
            if (me_ == owner(tij)) {
                madness::Future<madness::Tensor<double> > tile = rhs.task(me_, &Distributed_Matrix::get_tile_tij, tij);
                this->task(me_, &Distributed_Matrix::sum_tile_tij, tij, tile);
                //                           const_cast<double*>(&rhs.data_[local(tij)](0,0)));
            }
        }
    }
    else throw PSIEXCEPTION("The matrices being added are not the same.\n");
    Communicator::world->sync();
    return madness::None;
}

madness::Void Distributed_Matrix::operator +=(const Distributed_Matrix *rhs)
{
    Communicator::world->sync();
    if (*this == rhs) {
        for (int tij=0; tij < ntiles_; tij++) {
            if (me_ == owner(tij)) {
                madness::Future<madness::Tensor<double> > tile = rhs->task(me_, &Distributed_Matrix::get_tile_tij, tij);
                this->task(me_, &Distributed_Matrix::sum_tile_tij, tij, tile);
                //                           const_cast<double*>(&rhs->data_[local(tij)](0,0)));
            }
        }
    }
    else throw PSIEXCEPTION("The matrices being added are not the same.\n");
    Communicator::world->sync();
    return madness::None;
}


madness::Void Distributed_Matrix::sum_tile(const int &ti, const int &tj,
                                           //                                           const double *rhs_ptr)
                                           const madness::Tensor<double> &tile)
{
    //    sum_tile_ij(i*global_tile_ncols_+j, tile);
    sum_tile_tij(ti*tile_ncols_+tj, tile);
    return madness::None;
}
madness::Void Distributed_Matrix::sum_tile_tij(//const int &tij, const double *rhs_ptr)
                                               const int &tij, const madness::Tensor<double> &tile)
{
    int loc = local(tij);

    double *ptr = const_cast<double*>(this->data_[loc].ptr());
    double *end = ptr + this->data_[loc].size();
    double *rhs_ptr = const_cast<double*>(tile.ptr());

    add_mutex_->lock();
    while (ptr < end)
        *ptr++ += *rhs_ptr++;
    add_mutex_->unlock();

    return madness::None;
}

Distributed_Matrix Distributed_Matrix::operator +(const Distributed_Matrix &rhs) const
{
    if (*this == rhs) {
        Distributed_Matrix result = *this;
        result += rhs; //const_cast<Distributed_Matrix&>(rhs);
        return result;
    }
    else throw PSIEXCEPTION("The matrices that are being added are not the same.\n");
}

Distributed_Matrix Distributed_Matrix::operator +(const Distributed_Matrix *rhs) const
{
    if (*this == rhs) {
        Distributed_Matrix result = *this;
        result += rhs; //const_cast<Distributed_Matrix*>(rhs);
        return result;
    }
    else throw PSIEXCEPTION("The matrices that are being added are not the same.\n");
}

madness::Void Distributed_Matrix::scale(const double &val)
{
    Communicator::world->sync();
    if (this->nelements_) {
        for (int tij=0; tij < ntiles_; tij++) {
            if (me_ == owner(tij)) {
                this->task(me_, &Distributed_Matrix::scale_tile_tij, tij, val);
            }
        }
    }
    else throw PSIEXCEPTION("The matrix you are trying to scale is empty.\n");

    Communicator::world->sync();
    return madness::None;
}
madness::Void Distributed_Matrix::operator*=(const double &val)
{
    this->scale(val);
    return madness::None;
}

madness::Void Distributed_Matrix::scale_tile(const int &ti, const int &tj,
                                             const double &val)
{
    scale_tile_tij(ti*tile_ncols_+tj, val);
    return madness::None;
}
madness::Void Distributed_Matrix::scale_tile_tij(const int &tij, const double &val)
{
    int loc = local(tij);
    //    this->data_[loc].scale(val);

    double *ptr = this->data_[loc].ptr();
    double *end = ptr + this->data_[loc].size();

    while (ptr < end)
        *ptr++ *= val;

    return madness::None;
}

double Distributed_Matrix::trace() const
{
    Communicator::world->sync();
    double trace_val = 0.0;
    for (int ti=0; ti < tile_nrows_; ti++) {
        if (me_ == owner(ti,ti)) {
            madness::Future<double> val = task(me_, &Distributed_Matrix::trace_tile, ti, ti);
            trace_val += val.get();
        }
    }
    Communicator::world->sync();
    Communicator::world->sum(&trace_val, 1);

    return trace_val;
}

double Distributed_Matrix::trace_tile(const int &ti, const int &tj)
{
    return trace_tile_tij(ti*tile_ncols_ + tj);
}

double Distributed_Matrix::trace_tile_tij(const int &tij)
{
    int loc = local(tij);
    int stride = data_[loc].dim(1)+1;
    double val = 0.0;
    double *ptr = this->data_[loc].ptr();
    double *end = ptr + this->data_[loc].size();

    while (ptr < end) {
        val += *ptr;
        ptr+=stride;
    }

    return val;
}



double Distributed_Matrix::vector_dot(const Distributed_Matrix &rmat)
{
    Communicator::world->sync();
    double dot = 0.0;
    if (*this == rmat) {
        for (int tij=0; tij < ntiles_; tij++) {
            if (me_ == owner(tij)) {

                //                madness::Future<madness::Tensor<double> > rtile = rmat.task(me_, &Distributed_Matrix::get_tile_ij, t);
                //                madness::Future<madness::Tensor<double> > ltile = this->task(me_, &Distributed_Matrix::get_tile_ij, t);
                //                dot += task(me_, &Distributed_Matrix::vector_dot_tile, ltile, rtile);

                int size = data_[local(tij)].size();
                double *lptr = data_[local(tij)].ptr();
                double *rptr = const_cast<double*>(rmat.data_[local(tij)].ptr());
                madness::Future<double> val = task(me_, &Distributed_Matrix::vector_dot_tile, rptr, rptr, size);
                dot += val.get();
            }
        }
    }
    else throw PSIEXCEPTION("The matrices are not the same.\n");

    Communicator::world->sync();
    Communicator::world->sum(&dot, 1);
    return dot;
}

double Distributed_Matrix::vector_dot_tile(const double *lptr, const double *rptr, const int &size)
//double Distributed_Matrix::vector_dot_tile(const madness::Tensor<double> &ltile,
//                                           const madness::Tensor<double> &rtile)
{
    return C_DDOT(size, const_cast<double*>(lptr), 1,
                  const_cast<double*>(rptr), 1);

    //    return C_DDOT(ltile.size(), const_cast<double*>(ltile.ptr()), 1,
    //                  const_cast<double*>(rtile.ptr()), 1);
}




//madness::Void Distributed_Matrix::copy_from_tile(const int &row, const int &col,
//                                                 const int &inc, const madness::Tensor<double> tile)
//{
//    for (int i=0, index = row * this->ncols_ + col;
//         i < tile.size(); i++, index += inc) {
//        int r = this->convert_ij_to_i(index);
//        int c = this->convert_ij_to_j(index);
//        int tr = i/tile.dim(1);
//        int tc = i%tile.dim(1);

//        if (this->me_ == owner(r/this->tile_sz_,c/this->tile_sz_))
//            (*this)[r][c] = tile(tr,tc);

//    }
//    return madness::None;
//}

//madness::Void Distributed_Matrix::dswap(const int &length, Distributed_Matrix &X,
//                                        const int &xrow, const int &xcol, const int &x_inc,
//                                        const int &yrow, const int &ycol, const int &y_inc)
//{

//    Communicator::world->sync();

//    int yr,yc,xr,xc;
//    for (int i=0, y_index = yrow * this->ncols_ + ycol,
//         x_index = xrow * X.ncols_ + xcol;
//         i < length; i++, y_index += y_inc, x_index = x_inc) {

//        int yr = this->convert_ij_to_i(y_index); // y_index/this->ncols_;
//        int yc = this->convert_ij_to_j(y_index); //  y_index%this->ncols_;
//        int xr = X.convert_ij_to_i(x_index); //x_index/X.ncols_;
//        int xc = X.convert_ij_to_j(x_index); //x_index%X.ncols_;

//        madness::Future<double> new_y, new_x;
//        if (this->me_ == this->owner(yr/this->tile_sz_,yc/this->tile_sz_))
//            new_y = X[xr][xc].get_future();
//        if (X.me_ == X.owner(yr/X.tile_sz_,yc/X.tile_sz_))
//            new_x = (*this)[yr][yc].get_future();

//        // make sure that the new values have been sent to new node
//        Communicator::world->sync();

//        if (this->me_ == this->owner(yr/this->tile_sz_,yc/this->tile_sz_))
//            (*this)[yr][yc] = new_y;
//        if (X.me_ == X.owner(yr/X.tile_sz_,yc/X.tile_sz_))
//            X[xr][xc] = new_x;

//    }


////    int y_row_inc = y_inc/nrows_;
////    int y_col_inc = y_inc%ncols_;

////    int x_row_inc = x_inc/nrows_;
////    int x_col_inc = x_inc%ncols_;

////    for (int i=0, xr=xrow, xc=xcol, yr=yrow, yc=ycol;
////         i < length; i++, xr+=x_row_inc, xc+=x_col_inc, yr+=y_row_inc, yc+=y_col_inc) {
////        if (me_ == owner(rc_tile(yr),rc_tile(yc)))
////            (*this)[yr][yc] = X[xr][xc].get_future();
////        if (me_ == owner(rc_tile(xr),rc_tile(xc)))
////            X[xr][xc] = (*this)[yr][yc].get_future();

////    }
//    Communicator::world->sync();

//    return madness::None;

//}

//double Distributed_Matrix::ddot(const int &length, Distributed_Matrix &X,
//                                const int x[2], const int &x_inc,
//                                const int y[2], const int &y_inc)
//{

//    Communicator::world->sync();
//    int y_row_inc = y_inc/nrows_;
//    int y_col_inc = y_inc%ncols_;

//    int x_row_inc = x_inc/nrows_;
//    int x_col_inc = x_inc%ncols_;

//    double val = 0.0;

//    int start = x[0]*X.ncols_ + x[1];
//    int end = start + length;

//    for (int i=start; i < end; i++) {

//    }


//    for (int i=0, xr=x[0], xc=x[1], yr=y[0], yc=y[1];
//         i < length; i++, xr+=x_row_inc, xc+=x_col_inc, yr+=y_row_inc, yc+=y_col_inc) {
//        if (me_ == owner(rc_to_tile(yr),rc_to_tile(yc))) {
//            madness::Future<double> xval = X[xr][xc].get_future();
//            madness::Future<double> yval = (*this)[yr][yc].get_future();
//            val += xval.get()*yval.get();
//        }
//    }
//    Communicator::world->sync();

//    Communicator::world->sum(&val, 1);

//    return val;

//}


//Distributed_Matrix Distributed_Matrix::operator* (Distributed_Matrix &B_mat)
//{
//    if (this->ncols_ == B_mat.nrows_) {
//        if (this->tile_sz_ == B_mat.tile_sz_) {
//            madness::TaskAttributes attr;
//            attr.set_highpriority(true);

//            Distributed_Matrix C_mat(pgrid_, this->nrows_, B_mat.ncols_,
//                                     this->tile_sz_, "Multiply Result");

//            C_mat.zero();

//            int A_tile_nrow = this->tile_nrows_;
//            int A_tile_ncol = this->tile_ncols_;
//            int B_tile_ncol = B_mat.tile_ncols_;
//            int C_tile_ncol = C_mat.tile_ncols_;

//            for (int ti=0; ti < A_tile_nrow; ti++) {

//                std::vector<madness::Future<madness::Tensor<double> > > A =
//                        madness::future_vector_factory<madness::Tensor<double> >(A_tile_ncol);

//                std::vector<int> proc_row = C_mat.pgrid_row(ti);
//                for (int prow=0; prow < proc_row.size(); prow++) {
//                    if (C_mat.me_ == proc_row[prow]) {
//                        for (int tj=0; tj < A_tile_ncol; tj++) {
//                            int tij = ti*A_tile_ncol + tj;
//                            A[tj] = this->task(this->owner(tij), &Distributed_Matrix::get_tile_tij, tij);
//                        }
//                    }
//                }

//                for (int tk=0; tk < B_tile_ncol; tk++) {
//                    int tik = ti*C_tile_ncol + tk;
//                    if (C_mat.me_ == C_mat.owner(tik)) {

//                        for (int tj=0; tj < A_tile_ncol; tj++) {
//                            int tjk = tj*B_tile_ncol + tk;
//                            madness::Future<madness::Tensor<double> > B = B_mat.task(B_mat.owner(tjk), &Distributed_Matrix::get_tile_tij, tjk);
//                            madness::Future<madness::Tensor<double> > temp = C_mat.task(C_mat.owner(tik), &Distributed_Matrix::mxm, A[tj], B);
//                            task(C_mat.owner(tik), &Distributed_Matrix::sum_tile_tij, tik, temp);
//                        }
//                    }
//                }
//            }
//            Communicator::world->sync();

//            C_mat.print("C_mat test mxm.\n");

//            return C_mat;
//        }
//        else
//            throw PSIEXCEPTION("The tile sizes of A and B do not match.\n");
//    }
//    else {
//        throw PSIEXCEPTION("The columns of A do not match the rows of B.\n");
//    }


//}

//madness::Tensor<double> Distributed_Matrix::mxm(const madness::Tensor<double> &A,
//                                                const madness::Tensor<double> &B)
//{

//    madness::Tensor<double> temp(A.dim(0),B.dim(1));

//    C_DGEMM('n', 'n', A.dim(0), B.dim(1), A.dim(1), 1.0,
//            const_cast<double*>(A.ptr()), A.dim(1),
//            const_cast<double*>(B.ptr()), B.dim(1), 0.0,
//            temp.ptr(), B.dim(1));

//    return temp;
//}


Distributed_Matrix Distributed_Matrix::operator* (Distributed_Matrix &b_mat)
{

    if (this->ncols_ == b_mat.nrows_) {
        if (this->tile_sz_ == b_mat.tile_sz_) {
            madness::TaskAttributes attr;
            attr.set_highpriority(true);

            Distributed_Matrix c_mat(pgrid_, this->nrows_, b_mat.ncols_,
                                     this->tile_sz_, "Multiply Result");

            c_mat.zero();

            int a_mat_row = this->tile_nrows_;
            int a_mat_col = this->tile_ncols_;
            int b_mat_col = b_mat.tile_ncols_;
            int c_mat_col = c_mat.tile_ncols_;

            for (int i=0; i < a_mat_row; i++) {
                for (int j=0; j < b_mat_col; j++) {
                    int ij = i*c_mat_col + j;

                    if (me_ == owner(ij)) {
                        for (int k=0; k < a_mat_col; k++) {
                            int ik = i*a_mat_col + k;
                            int kj = k*b_mat_col + j;
                            madness::Future<madness::Tensor<double> > A = this->task(this->owner(kj), &Distributed_Matrix::get_tile_tij, kj);
                            madness::Future<madness::Tensor<double> > B = b_mat.task(b_mat.owner(kj), &Distributed_Matrix::get_tile_tij, kj);
                            c_mat.task(me_, &Distributed_Matrix::mxm, A, B, ik, 1.0);
                        }
                    }
                }
            }

            Communicator::world->sync();

            return c_mat;
        }
        else
            throw PSIEXCEPTION("The tile sizes of A and B do not match.\n");
    }
    else {
        throw PSIEXCEPTION("The columns of A do not match the rows of B.\n");
    }


//    if (this->ncols_ == b_mat.nrows_) {
//        if (this->tile_sz_ == b_mat.tile_sz_) {
//            madness::TaskAttributes attr;
//            attr.set_highpriority(true);

//            Distributed_Matrix c_mat(pgrid_, this->nrows_, b_mat.ncols_,
//                                     this->tile_sz_, "Multiply Result");

//            c_mat.zero();

//            int a_mat_row = this->tile_nrows_;
//            int a_mat_col = this->tile_ncols_;
//            int b_mat_col = b_mat.tile_ncols_;
//            int c_mat_col = c_mat.tile_ncols_;

//            for (int i=0; i < a_mat_row; i++) {

//                std::vector<madness::Future<madness::Tensor<double> > > A =
//                        madness::future_vector_factory<madness::Tensor<double> >(a_mat_col);

//                std::vector<int> proc_row = c_mat.pgrid_row(i);
//                for (int prow=0; prow < proc_row.size(); prow++) {
//                    if (c_mat.me_ == proc_row[prow]) {
//                       for (int j=0; j < a_mat_col; j++) {
//                            int ij = i*a_mat_col + j;
//                            A[j] = this->task(this->owner(ij), &Distributed_Matrix::get_tile_tij, ij);
//                        }
//                    }
//                }


//                for (int k=0; k < b_mat_col; k++) {
//                    int ik = i*c_mat_col + k;

//                    if (c_mat.me_ == c_mat.owner(ik)) {

//                        for (int j=0; j < a_mat_col; j++) {
//                            int jk = j*b_mat_col + k;
//                            madness::Future<madness::Tensor<double> > B = b_mat.task(b_mat.owner(jk), &Distributed_Matrix::get_tile_tij, jk);
//                            c_mat.task(c_mat.owner(ik), &Distributed_Matrix::mxm,
//                                       A[j], B, ik, 1.0);
//                        }
//                    }
//                }
//            }

//            Communicator::world->sync();

//            return c_mat;
//        }
//        else
//            throw PSIEXCEPTION("The tile sizes of A and B do not match.\n");
//    }
//    else {
//        throw PSIEXCEPTION("The columns of A do not match the rows of B.\n");
//    }
}

madness::Void Distributed_Matrix::mxm(const madness::Tensor<double> &a,
                                      const madness::Tensor<double> &b,
                                      const int &c,
                                      const double &c_scale)
{
    mult_mutex_->lock();

    C_DGEMM('n', 'n', a.dim(0), b.dim(1), a.dim(1), 1.0,
            const_cast<double*>(a.ptr()), a.dim(1),
            const_cast<double*>(b.ptr()), b.dim(1), c_scale,
            const_cast<double*>(this->data_[local(c)].ptr()), b.dim(1));
    mult_mutex_->unlock();

    return madness::None;
}




} // End of namespace psi

#endif
