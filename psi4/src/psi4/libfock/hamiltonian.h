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
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

namespace psi {

class Matrix;
class Vector;
class JK;
class VBase;

// => BASE CLASSES <= //

class Hamiltonian {

protected:

    /// Print flag, defaults to 1
    int print_;
    /// Debug flag, defaults to 0
    int debug_;
    /// Bench flag, defaults to 0
    int bench_;
    /// Use exact diagonal, if available?
    bool exact_diagonal_;
    /// jk object
    std::shared_ptr<JK> jk_;
    /// v object
    std::shared_ptr<VBase> v_;

    void common_init();

public:
    // => Constructors < = //

    Hamiltonian(std::shared_ptr<JK> jk);
    Hamiltonian(std::shared_ptr<JK> jk, std::shared_ptr<VBase> v);
    /// Destructor
    virtual ~Hamiltonian();

    // => Accessors <= //

    /**
    * Pointer to the JK object
    * @return current JK object
    */
    std::shared_ptr<JK> jk() const { return jk_; }
    /**
    * Pointer to the V object
    * @return current VBase object
    */
    std::shared_ptr<VBase> v() const { return v_; }

    /**
    * Print header information regarding Hamiltonian
    * type on output file
    */
    virtual void print_header() const = 0;

    // => Knobs <= //

    /**
    * Knob to swap out a JK object
    * @param jk new JK object
    */
    void set_JK(std::shared_ptr<JK> jk) { jk_ = jk; }
    /**
    * Knob to swap out a V object
    * @param v new V object
    */
    void set_V(std::shared_ptr<VBase> v) { v_ = v; }
    /// Print flag (defaults to 1)
    void set_print(int print) { print_ = print; }
    /// Debug flag (defaults to 0)
    void set_debug(int debug) { debug_ = debug; }
    /// Bench flag (defaults to 0)
    void set_bench(int bench) { bench_ = bench; }
    /// User the exact diagonal, if available? (defaults to false)
    void set_exact_diagonal(bool diag) { exact_diagonal_ = diag; }
};

class RHamiltonian : public Hamiltonian {

public:
    // => Constructors < = //

    RHamiltonian(std::shared_ptr<JK> jk);
    RHamiltonian(std::shared_ptr<JK> jk, std::shared_ptr<VBase> v);
    /// Destructor
    virtual ~RHamiltonian();

    // => Required Methods <= //

    /**
    * Return the approximate diagonal of the Hamiltonian
    * (typically orbital energy differences). This is used
    * for vector guess, eigenvector occupation, and preconditioning.
    * @return diagonal approximation, blocked by symmetry
    */
    virtual std::shared_ptr<Vector> diagonal() = 0;
    /**
    * Form the product Hx for each x found in the argument, placing in second argument
    * @param x vector of state functions, blocked by symmetry
    * @param b vector of product functions, blocked by symmetry (preallocated)
    */
    virtual void  product(const std::vector<std::shared_ptr<Vector> >& x,
                                std::vector<std::shared_ptr<Vector> >& b) = 0;

    /**
    * Form the explicit hamiltonian for debugging purposes
    */
    SharedMatrix explicit_hamiltonian();

};

class UHamiltonian : public Hamiltonian {

public:
    // => Constructors < = //

    UHamiltonian(std::shared_ptr<JK> jk);
    UHamiltonian(std::shared_ptr<JK> jk, std::shared_ptr<VBase> v);
    /// Destructor
    virtual ~UHamiltonian();

    // => Required Methods <= //

    /**
    * Return the approximate diagonal of the Hamiltonian
    * (typically orbital energy differences). This is used
    * for vector guess, eigenvector occupation, and preconditioning.
    * @return diagonal approximation, \alpha and \beta vectors,
    * blocked by symmetry
    */
    virtual std::pair<std::shared_ptr<Vector>,
                      std::shared_ptr<Vector> > diagonal() = 0;
    /**
    * Form the product Hx for each x found in the argument, placing in second argument
    * @param x vector of state functions, \alpha and \beta, blocked by symmetry
    * @param b vector of product functions, \alpha and \beta, blocked by symmetry (preallocated)
    */
    virtual void product(const std::vector<std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > >& x,
                               std::vector<std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > >& b) = 0;
    /**
    * Form the explicit hamiltonian for debugging purposes
    */
    std::pair<SharedMatrix, SharedMatrix > explicit_hamiltonian();
};

// => APPLIED CLASSES <= //

class MatrixRHamiltonian : public RHamiltonian {

protected:

    SharedMatrix M_;

public:

    MatrixRHamiltonian(SharedMatrix M);
    virtual ~MatrixRHamiltonian();

    virtual void print_header() const;
    virtual std::shared_ptr<Vector> diagonal();
    virtual void product(const std::vector<std::shared_ptr<Vector> >& x,
                               std::vector<std::shared_ptr<Vector> >& b);

};

class MatrixUHamiltonian : public UHamiltonian {

protected:

    std::pair<SharedMatrix, SharedMatrix > M_;

public:

    MatrixUHamiltonian(std::pair<SharedMatrix, SharedMatrix > M);
    virtual ~MatrixUHamiltonian();

    virtual void print_header() const;
    virtual std::pair<std::shared_ptr<Vector>,
                      std::shared_ptr<Vector> > diagonal();
    virtual void product(const std::vector<std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > >& x,
                               std::vector<std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > >& b);

};

class CISRHamiltonian : public RHamiltonian {

protected:

    bool singlet_;
    SharedMatrix Caocc_;
    SharedMatrix Cavir_;
    std::shared_ptr<Vector> eps_aocc_;
    std::shared_ptr<Vector> eps_avir_;

public:
    CISRHamiltonian(std::shared_ptr<JK> jk,
                    SharedMatrix Caocc,
                    SharedMatrix Cavir,
                    std::shared_ptr<Vector> eps_aocc,
                    std::shared_ptr<Vector> eps_avir,
                    std::shared_ptr<VBase> v = std::shared_ptr<VBase>());
    virtual ~CISRHamiltonian();

    virtual void print_header() const;
    virtual std::shared_ptr<Vector> diagonal();
    virtual void product(const std::vector<std::shared_ptr<Vector> >& x,
                               std::vector<std::shared_ptr<Vector> >& b);

    virtual std::vector<SharedMatrix > unpack(const std::shared_ptr<Vector>& x);

    void set_singlet(bool singlet) { singlet_ = singlet; }

};

class TDHFRHamiltonian : public RHamiltonian {

protected:

    bool singlet_;
    SharedMatrix Caocc_;
    SharedMatrix Cavir_;
    std::shared_ptr<Vector> eps_aocc_;
    std::shared_ptr<Vector> eps_avir_;

public:
    TDHFRHamiltonian(std::shared_ptr<JK> jk,
                    SharedMatrix Caocc,
                    SharedMatrix Cavir,
                    std::shared_ptr<Vector> eps_aocc,
                    std::shared_ptr<Vector> eps_avir,
                    std::shared_ptr<VBase> v = std::shared_ptr<VBase>());
    virtual ~TDHFRHamiltonian();

    virtual void print_header() const;
    virtual std::shared_ptr<Vector> diagonal();
    virtual void product(const std::vector<std::shared_ptr<Vector> >& x,
                               std::vector<std::shared_ptr<Vector> >& b);

    void set_singlet(bool singlet) { singlet_ = singlet; }
};

class CPHFRHamiltonian : public RHamiltonian {

protected:

    SharedMatrix Caocc_;
    SharedMatrix Cavir_;
    std::shared_ptr<Vector> eps_aocc_;
    std::shared_ptr<Vector> eps_avir_;

public:
    CPHFRHamiltonian(std::shared_ptr<JK> jk,
                     SharedMatrix Caocc,
                     SharedMatrix Cavir,
                     std::shared_ptr<Vector> eps_aocc,
                     std::shared_ptr<Vector> eps_avir,
                     std::shared_ptr<VBase> v = std::shared_ptr<VBase>());
    virtual ~CPHFRHamiltonian();

    virtual void print_header() const;
    virtual std::shared_ptr<Vector> diagonal();
    virtual void product(const std::vector<std::shared_ptr<Vector> >& x,
                               std::vector<std::shared_ptr<Vector> >& b);

    virtual std::map<std::string, SharedVector> pack(const std::map<std::string, std::shared_ptr<Matrix> >& b);
    virtual std::vector<SharedMatrix > unpack(const std::vector<SharedVector >& x);
};

class TDARHamiltonian : public CISRHamiltonian {

protected:

    SharedMatrix Cocc_;

public:
    TDARHamiltonian(std::shared_ptr<JK> jk,
                    std::shared_ptr<VBase> v,
                    SharedMatrix Cocc,
                    SharedMatrix Caocc,
                    SharedMatrix Cavir,
                    std::shared_ptr<Vector> eps_aocc,
                    std::shared_ptr<Vector> eps_avir);
    virtual ~TDARHamiltonian();

    virtual void print_header() const;
    virtual void product(const std::vector<std::shared_ptr<Vector> >& x,
                               std::vector<std::shared_ptr<Vector> >& b);
};

class TDDFTRHamiltonian : public TDHFRHamiltonian {

protected:

    SharedMatrix Cocc_;

public:
    TDDFTRHamiltonian(std::shared_ptr<JK> jk,
                    std::shared_ptr<VBase> v,
                    SharedMatrix Cocc,
                    SharedMatrix Caocc,
                    SharedMatrix Cavir,
                    std::shared_ptr<Vector> eps_aocc,
                    std::shared_ptr<Vector> eps_avir);
    virtual ~TDDFTRHamiltonian();

    virtual void print_header() const;
    virtual void product(const std::vector<std::shared_ptr<Vector> >& x,
                               std::vector<std::shared_ptr<Vector> >& b);
};

class CPKSRHamiltonian : public CPHFRHamiltonian {

protected:

    SharedMatrix Cocc_;

public:
    CPKSRHamiltonian(std::shared_ptr<JK> jk,
                    std::shared_ptr<VBase> v,
                    SharedMatrix Cocc,
                    SharedMatrix Caocc,
                    SharedMatrix Cavir,
                    std::shared_ptr<Vector> eps_aocc,
                    std::shared_ptr<Vector> eps_avir);
    virtual ~CPKSRHamiltonian();

    virtual void print_header() const;
    virtual void product(const std::vector<std::shared_ptr<Vector> >& x,
                               std::vector<std::shared_ptr<Vector> >& b);
};

// "Hamiltonian" for UHF stability analysis.
class USTABHamiltonian : public UHamiltonian {

protected:

    SharedMatrix Cocca_;
    SharedMatrix Cvira_;
    SharedMatrix Coccb_;
    SharedMatrix Cvirb_;
    std::shared_ptr<Vector> eps_occa_;
    std::shared_ptr<Vector> eps_vira_;
    std::shared_ptr<Vector> eps_occb_;
    std::shared_ptr<Vector> eps_virb_;

public:

// Should really use a map of matrices here. But meh.
    USTABHamiltonian(std::shared_ptr<JK> jk,
                     SharedMatrix Cocca,
                     SharedMatrix Cvira,
                     SharedMatrix Coccb,
                     SharedMatrix Cvirb,
                     std::shared_ptr<Vector> eps_occa,
                     std::shared_ptr<Vector> eps_vira,
                     std::shared_ptr<Vector> eps_occb,
                     std::shared_ptr<Vector> eps_virb,
                     std::shared_ptr<VBase> v = std::shared_ptr<VBase>());
    virtual ~USTABHamiltonian();

    virtual void print_header() const;
    virtual std::pair<std::shared_ptr<Vector>,
                      std::shared_ptr<Vector> > diagonal();
    virtual void product(const std::vector<std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > >& x,
                               std::vector<std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> > >& b);

// Working with a pair is annoying, so we define a new function below
    virtual std::vector<std::pair<SharedMatrix,SharedMatrix > > unpack(
            const std::pair<std::shared_ptr<Vector>, std::shared_ptr<Vector> >& x){} ;
    virtual std::vector<std::pair<SharedMatrix,SharedMatrix > > unpack_paired(
            const std::shared_ptr<Vector>& x);

};

}
#endif
