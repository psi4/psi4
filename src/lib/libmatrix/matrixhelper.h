/*
 * File:   matrix.h
 * Author: Ben Mintz, and anyone else
 *
 * Created on March 14, 2012, 3:34 PM
 */

#ifndef _psi_src_lib_libmatrix_matrixhelper_h_
#define	_psi_src_lib_libmatrix_matrixhelper_h_

#include <psi4-dec.h>
#include <psiconfig.h>

#include <cstring>
#include <cstdio>

#include <libmatrix/matrixwrapper.h>

namespace psi {

extern FILE *outfile;


class CreateMatrix {
protected:
    typedef boost::shared_ptr<MatrixWrapper> matrixwrapper;

    matrixwrapper matrix;

public:


    CreateMatrix(int rows, int cols, int tile_sz=0)
    {
#if serial == 1
        matrix = matrixwrapper(new SerialMatrix(rows, cols, tile_sz));
#elif parallel == 1
        matrix = matrixwrapper(new PsiDistMatrix(rows, cols, 2));
#endif
    }

    ~CreateMatrix() { }

    matrixwrapper get() { return matrix; }

};

} // End of namespace psi

#endif  /* _psi_src_lib_libmatrix_matrixhelper_h_ */
