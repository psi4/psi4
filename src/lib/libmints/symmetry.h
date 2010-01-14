#ifndef _psi_src_lib_libmints_symmetry_h_
#define _psi_src_lib_libmints_symmetry_h_

#include <cstdio>
#include <libchkpt/chkpt.hpp>

#include <libmints/molecule.h>

namespace psi {

    class SymmetryOperation {
    private:
	double d_[3][3];

    public:
	SymmetryOperation();
	SymmetryOperation(const SymmetryOperation &);
	virtual ~SymmetryOperation();

	double trace() const { return d_[0][0] + d_[1][1] + d_[2][2]; }

	double* operator[](int i) { return d_[i]; }

	const double* operator[](int i) const { return d_[i]; }

	double& operator()(int i, int j) { return d_[i][j]; }

	double operator()(int i, int j) const { return d_[i][j]; }

	void zero() { memset(d_, 0, sizeof(double) * 9); }

	SymmetryOperation operate(const SymmetryOperation& r) const;

	SymmetryOperation transform(const SymmetryOperation& r) const;

	void unit() { zero(); d_[0][0] = d_[1][1] = d_[2][2] = 1.0; }

	void E() { unit(); }

	void i() { zero(); d_[0][0] = d_[1][1] = d_[2][2] = -1.0; }

	void sigma_h() { unit(); d_[2][2] = -1.0; }
	void sigma_xz() { unit(); d_[1][1] = -1.0; }
	void sigma_yz() { unit(); d_[0][0] = -1.0; }

	void rotation(int n);
	void rotation(double theta);

	void c2_x() { i(); d_[0][0] = 1.0; }
	void c2_y() { i(); d_[1][1] = 1.0; }

	void transpose();
    };

    class SymRep {
    private:
	int n_;
	double d_[5][5];

    public:
	SymRep(int =0);
	SymRep(const SymmetryOperation&);
	virtual ~SymRep();

	operator SymmetryOperation() const;

	inline double trace() const;

	void set_dim(int i) { n_ = i; }

	double* operator[](int i) { return d_[i]; }
	const double* operator[](int i) const { return d_[i]; }

	double& operator()(int i, int j) { return d_[i][j]; }
	double operator()(int i, int j) const { return d_[i][j]; }

	void zero() memset(d_, 0, sizeof(double) * 25); }

	SymRep operate(const SymRep& r) const;
}

#endif // _basis_symmetry_h_
