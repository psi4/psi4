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

#ifndef LRERI_H
#define LRERI_H

#include <map>
#include <string>
 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/libmints/typedefs.h"

namespace psi {

class Options;
class Tensor;
class CoreTensor;
class DiskTensor;
class BasisSet;
class LRERI {

protected:

    // => Utility <= //

    /// Print flag
    int print_;
    /// Debug flag
    int debug_;
    /// Bench flag
    int bench_;
    /// Memory in doubles
    unsigned long int memory_;

    // => Basis Set <= //

    /// Primary orbital basis (nso)
    std::shared_ptr<BasisSet> primary_;

    // => Orbital Spaces <= //

    /// Occupation matrix coefficients (nso x nmo)
    std::shared_ptr<Matrix> C_;
    /// Orbital spaces, each defined by a keyword and the index range in <start, end+1>.
    std::map<std::string, std::pair<int, int> > spaces_;
    /// Orbital spaces order buffer, to keep the printing nice.
    std::vector<std::string> spaces_order_;

    // => Utility Routines <= //

    /// Set defaults
    void common_init();

    /// Inverse fitting metric
    std::shared_ptr<Matrix> Jm12(std::shared_ptr<BasisSet> auxiliary, double condition);

public:
    // => Constructors <= //

    LRERI(std::shared_ptr<BasisSet> primary);
    virtual ~LRERI();

    // => Defaults <= //

    /// O: Load the usual orbital spaces for an RHF or UHF wavefunction
    virtual void load_wavefunction(std::shared_ptr<Wavefunction> ref);
    /// O: Load the usual options objects
    virtual void load_options(Options& options);

    // => Orbital Space Control <= //

    /// R: Set the overall C matrix (calls clear before starting)
    void set_C(std::shared_ptr<Matrix> C);
    /// R: Add an orbital subspace to the queue. The subspace will be referred to by key, and ranges from start to end-1.
    void add_space(const std::string& key, int start, int end);
    /// Clear the C matrix and orbital spaces list
    virtual void clear();

    // => Computers <= //

    /// Print info
    virtual void print_header(int level = 1) = 0;
    /// R: Compute the desired ERI factorization
    virtual void compute() = 0;

    // => Setters <= //

    /// Set the print flag
    void set_print(int print) { print_ = print; }
    /// Set the debug flag
    void set_debug(int debug) { debug_ = debug; }
    /// Set the bench flag
    void set_bench(int bench) { bench_ = bench; }
    /// Set the allowed memory in doubles
    void set_memory(unsigned long int memory) { memory_ = memory; }

};

/**
 * DFERI
 **/
class DFERI : public LRERI {

protected:

    // => DF Tech <= //

    /// Auxiliary orbital-pair basis
    std::shared_ptr<BasisSet> auxiliary_;
    /// Relative condition number in J^-1/2
    double J_cutoff_;
    /// Schwarz sieve tolerance
    double schwarz_cutoff_;

    // => HACK for LRC-ERIs <= //

    double omega_;

    // => Targets <= //

    /// Three-center integrals, by name, sorted e.g. (ov|Q), DiskTensor
    std::map<std::string, std::shared_ptr<Tensor> > ints_;
    /// Requested pair spaces
    std::map<std::string, std::pair<std::string, std::string> > pair_spaces_;
    /// Requested pair space powers
    std::map<std::string, double> pair_powers_;
    /// Requested pair space transpositions (ab|Q) -> (ba|Q)
    std::map<std::string, bool> pair_transposes_;
    /// Order of pair spaces, to keep printing nice
    std::vector<std::string> pair_spaces_order_;

    /// Keep the raw (Q|ia)-type integrals?
    bool keep_raw_integrals_;

    // => Utility Routines <= //

    /// Set defaults
    void common_init();

    void allocate();
    void transform();
    void fit();

public:
    // => Constructors <= //

    DFERI(std::shared_ptr<BasisSet> primary,
          std::shared_ptr<BasisSet> auxiliary);
    virtual ~DFERI();

    static std::shared_ptr<DFERI> build(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options);
    static std::shared_ptr<DFERI> build(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options, std::shared_ptr<Wavefunction> ref);

    // => Defaults <= //

    /// O: Load the usual options objects
    virtual void load_options(Options& options);

    // => Pair Space Control <= //

    // R: add an orbital pair space to the list of tasks
    void add_pair_space(const std::string& name, const std::string& space1, const std::string& space2, double pow = -1.0/2.0, bool transpose12 = false);
    // Clear this DFERI object of tasks
    virtual void clear();
    // Clear this DFERI object of pair spaces
    virtual void clear_pair_spaces();

    // => Computers <= //

    /// Print info
    virtual void print_header(int level = 1);
    /// The size of the auxiliary basis
    int size_Q();
    /// R: Compute the requested DF 3-index integrals
    virtual void compute();
    /// Handle to computed disk tensors, by name in add_pair above
    std::map<std::string, std::shared_ptr<Tensor> >& ints() { return ints_; }
    /// Return the J matrix raised to the desired power
    std::shared_ptr<Matrix> Jpow(double power = -1.0/2.0);

    // => Setters <= //

    /// Set the relative eigenvalue cutoff in J^{-1/2}
    void set_J_cutoff(double J_cutoff) { J_cutoff_ = J_cutoff; }
    /// Set the schwarz sieve cutoff
    void set_schwarz_cutoff(double schwarz_cutoff) { schwarz_cutoff_ = schwarz_cutoff; }
    /// Set to keep the raw integrals in (Q|ia) striping?
    void set_keep_raw_integrals(bool val) { keep_raw_integrals_ = val; }

    /// Set the omega value. ALL integrals and metrics will use this value if set.
    void set_omega(double omega) { omega_ = omega; }

};

/**
 * LSTHCERI
 **/
class LSTHCERI : public LRERI {


protected:

    // => LS-LSTHC Tech <= //

    /// Raw AO-basis X matrix (nso x nP, for now)
    std::shared_ptr<Matrix> X_;
    /// Auxiliary orbital-pair basis
    std::shared_ptr<BasisSet> auxiliary_;
    /// Relative condition number in J^-1/2
    double J_cutoff_;
    /// Relative condition number in S^-1
    double S_cutoff_;
    /// Schwarz sieve tolerance
    double schwarz_cutoff_;
    /// Balance the X matrices?
    bool balance_;

    // => Targets <= //

    /// ERI factors, CoreTensor, swapped out
    std::map<std::string, std::vector<std::shared_ptr<Tensor> > > ints_;
    /// METH factors, CoreTensor, swapped out
    std::map<std::string, std::vector<std::shared_ptr<Tensor> > > meths_;
    /// Requested ERI spaces
    std::map<std::string, std::vector<std::string> > eri_spaces_;
    /// Order of ERI spaces, to keep printing nice
    std::vector<std::string> eri_spaces_order_;

    // => Utility Routines <= //

    /// Set defaults
    void common_init();

    /// Build all requred X matrices (np x nP, core)
    std::map<std::string, std::shared_ptr<Tensor> > build_X(bool meth = false);
    /// Build all required E matrices (nA x nP, disk) [deleted]
    std::map<std::string, std::shared_ptr<Tensor> > build_E(std::map<std::string, std::shared_ptr<Tensor> >& Xs);
    /// Build all required inverse S matrices (nP x nP, core, swapped)
    std::map<std::string, std::shared_ptr<Tensor> > build_S(std::map<std::string, std::shared_ptr<Tensor> >& Xs, bool meth = false);
    /// Build all requred L matrices (nP x nA, core, swapped)
    std::map<std::string, std::shared_ptr<Tensor> > build_L(std::map<std::string, std::shared_ptr<Tensor> >& Es,
                                                              std::map<std::string, std::shared_ptr<Tensor> >& Ss);
    /// Build all required Z matrices (nP x nP, core, swapped)
    std::map<std::string, std::shared_ptr<Tensor> > build_Z(std::map<std::string, std::shared_ptr<Tensor> >& Ls);
    /// Pack up the integrals
    void pack(std::map<std::string, std::shared_ptr<Tensor> >& Xs,
              std::map<std::string, std::shared_ptr<Tensor> >& Zs,
              std::map<std::string, std::shared_ptr<Tensor> >& Ls,
              std::map<std::string, std::shared_ptr<Tensor> >& Ss);
    /// Pack up the meth intermediates
    void pack_meth(std::map<std::string, std::shared_ptr<Tensor> >& Xs,
                   std::map<std::string, std::shared_ptr<Tensor> >& Ss);

public:
    // => Constructors <= //

    LSTHCERI(std::shared_ptr<BasisSet> primary,
           std::shared_ptr<BasisSet> auxiliary,
           std::shared_ptr<Matrix> X);
    virtual ~LSTHCERI();

    static std::shared_ptr<LSTHCERI> build(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, std::shared_ptr<Matrix> X, Options& options);
    static std::shared_ptr<LSTHCERI> build(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, std::shared_ptr<Matrix> X, Options& options, std::shared_ptr<Wavefunction> ref);

    // => Defaults <= //

    /// O: Load the usual options objects
    virtual void load_options(Options& options);

    // => ERI Space Control <= //

    // R: add an eri space to the list of tasks
    void add_eri_space(const std::string& name, const std::string& space1, const std::string& space2, const std::string& space3, const std::string& space4);
    // Clear this LSTHCERI object of tasks
    virtual void clear();

    // => Computers <= //

    /// Print info
    virtual void print_header(int level = 1);
    /// R: Compute the requested LS-LSTHC factors
    virtual void compute();
    /// LS-LSTHC factors [X1,X2,Z,X3,X4,L12,L34,Sinv12,Sinv34]
    std::map<std::string, std::vector<std::shared_ptr<Tensor> > >& ints() { return ints_; };

    /// O: Compute the METH X and Sinv matrices
    virtual void compute_meth();
    /// METH helper factors [X1,X2,Sinv]
    std::map<std::string, std::vector<std::shared_ptr<Tensor> > >& meths() { return meths_; }

    // => Setters <= //

    /// Set the relative eigenvalue cutoff in J^{-1/2}
    void set_J_cutoff(double J_cutoff) { J_cutoff_ = J_cutoff; }
    /// Set the relative eigenvalue cutoff in S^{-1}
    void set_S_cutoff(double S_cutoff) { S_cutoff_ = S_cutoff; }
    /// Set the schwarz sieve cutoff
    void set_schwarz_cutoff(double schwarz_cutoff) { schwarz_cutoff_ = schwarz_cutoff; }
    /// Set to balance or not?
    void set_balance(bool balance) { balance_ = balance; }
};


} // End namespace

#endif
