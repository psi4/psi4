/*
 * File:   matrix.h
 * Author: Ben Mintz, and anyone else
 *
 * Created on March 14, 2012, 3:34 PM
 */

#ifndef _psi_src_lib_libmatrix_matrixwrapper_h_
#define	_psi_src_lib_libmatrix_matrixwrapper_h_

#include <psi4-dec.h>
#include <psiconfig.h>

#include <cstring>
#include <cstdio>

#include <libparallel/parallel.h>
#include <libmints/matrix.h>
#include <libdist_matrix/dist_mat.h>

#if HAVE_MADNESS
    #define parallel 1
#else
    #define serial 1
#endif

namespace psi {

extern FILE *outfile;

class MatrixWrapper;
class SerialMatrix;
class PsiDistMatrix;


class MatrixWrapper {
protected:
    SharedMatrix serialmat_;
#if parallel == 1
    boost::shared_ptr<Distributed_Matrix> psidistmat_;
#endif

public:

    int me_;
    int nprocs_;
    int nthreads_;
    std::string comm_;

    MatrixWrapper();
    virtual ~MatrixWrapper();

    virtual void print() = 0;

    virtual void fill(double val) = 0;

    virtual void add(boost::shared_ptr<MatrixWrapper>) = 0;

    friend class SerialMatrix;
    friend class PsiDistMatrix;
};

class SerialMatrix : public MatrixWrapper {

public:
    SerialMatrix() {/* does nothing */ }
    SerialMatrix(int rows, int cols, int tile_sz=0);

    virtual ~SerialMatrix();

    virtual void print();

    virtual void fill(double val);

    virtual void add(boost::shared_ptr<MatrixWrapper> rhs);

};

#if parallel

class PsiDistMatrix : public MatrixWrapper {

public:

    PsiDistMatrix() {/* does nothing */ }
    PsiDistMatrix(int rows, int cols, int tile_sz=100);

    virtual ~PsiDistMatrix();

    virtual void fill(double val);

    virtual void print();

    virtual void add(boost::shared_ptr<MatrixWrapper> rhs);


};

#endif // End of parallel

}

#endif  /* _psi_src_lib_libmatrix_matrixwrapper_h_ */
