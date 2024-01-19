/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <vector>

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
    ~RHamiltonian() override;

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
    virtual void product(const std::vector<std::shared_ptr<Vector> >& x, std::vector<std::shared_ptr<Vector> >& b) = 0;

    /**
    * Form the explicit hamiltonian for debugging purposes
    */
    SharedMatrix explicit_hamiltonian();
};

// => APPLIED CLASSES <= //

class CPHFRHamiltonian : public RHamiltonian {
   protected:
    SharedMatrix Caocc_;
    SharedMatrix Cavir_;
    std::shared_ptr<Vector> eps_aocc_;
    std::shared_ptr<Vector> eps_avir_;

   public:
    CPHFRHamiltonian(std::shared_ptr<JK> jk, SharedMatrix Caocc, SharedMatrix Cavir, std::shared_ptr<Vector> eps_aocc,
                     std::shared_ptr<Vector> eps_avir, std::shared_ptr<VBase> v = std::shared_ptr<VBase>());
    ~CPHFRHamiltonian() override;

    void print_header() const override;
    std::shared_ptr<Vector> diagonal() override;
    void product(const std::vector<std::shared_ptr<Vector> >& x, std::vector<std::shared_ptr<Vector> >& b) override;

    virtual std::map<std::string, SharedVector> pack(const std::map<std::string, std::shared_ptr<Matrix> >& b);
    virtual std::vector<SharedMatrix> unpack(const std::vector<SharedVector>& x);
};
}
#endif
