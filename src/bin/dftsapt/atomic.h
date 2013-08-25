/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef ASAPT_ATOMIC_H
#define ASAPT_ATOMIC_H

#include <libmints/typedefs.h>
#include <libmints/wavefunction.h>
#include <map>
#include <set>

namespace psi {

class MolecularGrid;

class AtomicDensity {

protected:

    // ==> Utility Metadata <== //
    
    /// Print flag
    int print_;
    /// Debug flag
    int debug_;
    /// Bench flag
    int bench_;

    // ==> Input Specification <== //

    /// Molecule this atomic density is built around
    boost::shared_ptr<Molecule> molecule_;
    /// Basis this atomic density is built around
    boost::shared_ptr<BasisSet> primary_;
    /// Density matrix this atomic density is built around
    boost::shared_ptr<Matrix> D_;
    /// MolecularGrid to collocate on (Becke-style, yo?)
    boost::shared_ptr<MolecularGrid> grid_;

    // ==> Targets <== //

    /// X coordinates of grid
    boost::shared_ptr<Vector> x_;
    /// Y coordinates of grid
    boost::shared_ptr<Vector> y_;
    /// Z coordinates of grid
    boost::shared_ptr<Vector> z_;
    /// Weights of grid
    boost::shared_ptr<Vector> w_;
    /// Total density of grid
    boost::shared_ptr<Vector> rho_;
    /// Matrix of atomic densities (true atoms x grid points)
    boost::shared_ptr<Matrix> Q_;

    // ==> Utility Routines <== //

    /// Ubiqitous header
    virtual void print_header() const = 0;
   
    /// Protected constructor
    AtomicDensity();

public:
    // ==> Master Routines <== //

    /// Master destructor
    virtual ~AtomicDensity();

    /// Master builder
    static boost::shared_ptr<AtomicDensity> build(const std::string& type, boost::shared_ptr<BasisSet> basis, boost::shared_ptr<Matrix> D, Options& options);

    /// Master compute routine
    virtual void compute() = 0;

    // ==> Accessors <== //

    boost::shared_ptr<Vector> x()   const { return x_; }
    boost::shared_ptr<Vector> y()   const { return y_; }
    boost::shared_ptr<Vector> z()   const { return z_; }
    boost::shared_ptr<Vector> w()   const { return w_; }
    boost::shared_ptr<Vector> rho() const { return rho_; }
    boost::shared_ptr<Matrix> Q()   const { return Q_; }

    boost::shared_ptr<Molecule> molecule() const { return molecule_; }
    boost::shared_ptr<BasisSet> primary() const { return primary_; }
    boost::shared_ptr<Matrix> D() const { return D_; }
    boost::shared_ptr<MolecularGrid> grid() const { return grid_; }

    // ==> Setters <== //
    
    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }
    void set_bench(int bench) { bench_ = bench; }

};

class StockholderDensity : public AtomicDensity {

protected:

    /// Ubiqitous header
    virtual void print_header() const;
    /// Protected constructor
    StockholderDensity();
    
public:
    /// Subclass destructor
    virtual ~StockholderDensity();

    /// Master compute routine
    virtual void compute();
};

} // End namespace

#endif

