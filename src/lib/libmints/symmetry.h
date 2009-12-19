#ifndef _psi_src_lib_libmints_symmetry_h_
#define _psi_src_lib_libmints_symmetry_h_

#include <cstdio>
#include <libchkpt/chkpt.hpp>

#include <libmints/molecule.h>

namespace psi {

/*! \ingroup MINTS */
class SymmOp {
private:
    double rep_[3][3];
public:
    // default constructor
    SymmOp();
    // copy constructor
    SymmOp(const SymmOp &);
    // destructor
    ~SymmOp();

    double trace() const { return rep_[0][0] + rep_[1][1] + rep_[2][2]; }
    
    double*       operator[](int a)       { return rep_[a]; }
    const double* operator[](int a) const { return rep_[a]; }

    double& operator()(int a, int b)       { return rep_[a][b]; }
    double  operator()(int a, int b) const { return rep_[a][b]; }

    void zero()     { memset(rep_, 0, sizeof(double)*9); }
    void E()        { zero(); rep_[0][0] = rep_[1][1] = rep_[2][2] =  1.0; }
    void i()        { zero(); rep_[0][0] = rep_[1][1] = rep_[2][2] = -1.0; }
    void sigma_h()  { E(); rep_[2][2] = -1.0; }
    void sigma_xz() { E(); rep_[1][1] = -1.0; }
    void sigma_yz() { E(); rep_[0][0] = -1.0; }
};

}

#endif // _basis_symmetry_h_
