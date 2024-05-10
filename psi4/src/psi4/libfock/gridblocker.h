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

#ifndef libmints_gridblocker_H
#define libmints_gridblocker_H

#include "psi4/psi4-dec.h"

#include "psi4/libmints/vector3.h"
#include "psi4/libpsi4util/exception.h"

namespace psi {

class BasisExtents;
class BlockOPoints;

/**
 * Class to determine groupings of DFT grid points for
 * efficient sparse evaluation of density
 */
class GridBlocker {
   protected:
    int debug_;
    int print_;
    int bench_;

    // Reference to previous grid layout
    const int npoints_ref_;
    double const* x_ref_;
    double const* y_ref_;
    double const* z_ref_;
    double const* w_ref_;

    const size_t tol_max_points_;
    const size_t tol_min_points_;
    const double tol_max_radius_;
    std::shared_ptr<BasisExtents> extents_;

    // New grid layout (built in blocks -- method specific)
    int npoints_;
    int max_points_;
    int max_functions_;
    // The total collocation size
    size_t collocation_size_;
    double* x_;
    double* y_;
    double* z_;
    double* w_;
    std::vector<std::shared_ptr<BlockOPoints>> blocks_;

   public:
    GridBlocker(const int npoints_ref, double const* x_ref, double const* y_ref, double const* z_ref,
                double const* w_ref, const int max_points, const int min_points, const double max_radius,
                std::shared_ptr<BasisExtents> extents);
    virtual ~GridBlocker();

    virtual void block() = 0;

    int npoints() const { return npoints_; }
    int max_points() const { return max_points_; }
    int max_functions() const { return max_functions_; }
    int collocation_size() const { return collocation_size_; }
    double* x() const { return x_; }
    double* y() const { return y_; }
    double* z() const { return z_; }
    double* w() const { return w_; }
    const std::vector<std::shared_ptr<BlockOPoints>>& blocks() const { return blocks_; }
    virtual const std::vector<std::vector<std::shared_ptr<BlockOPoints>>>& atomic_blocks() {
        throw PSIEXCEPTION("GridBlocker: Atomic blocks not implemented in parent class.");
    }

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }
    void set_bench(int bench) { bench_ = bench; }
};

/**
 * Atomic blocking
 */
class AtomicGridBlocker : public GridBlocker {
   public:
    AtomicGridBlocker(const int npoints_ref, double const* x_ref, double const* y_ref, double const* z_ref,
                      double const* w_ref, const int max_points, const int min_points, const double max_radius,
                      std::shared_ptr<BasisExtents> extents, const std::shared_ptr<Molecule> molecule,
                      const std::vector<std::vector<MassPoint>> atomic_grids);
    ~AtomicGridBlocker() override;

    void block() override;
    std::shared_ptr<Molecule> molecule_;
    std::vector<std::vector<MassPoint>> atomic_grids_;
    std::vector<std::vector<std::shared_ptr<BlockOPoints>>> atomic_blocks_;
    const std::vector<std::vector<std::shared_ptr<BlockOPoints>>>& atomic_blocks() override { return atomic_blocks_; };
};

/**
 * Naive stride-based blocking
 */
class NaiveGridBlocker : public GridBlocker {
   public:
    NaiveGridBlocker(const int npoints_ref, double const* x_ref, double const* y_ref, double const* z_ref,
                     double const* w_ref, const int max_points, const int min_points, const double max_radius,
                     std::shared_ptr<BasisExtents> extents);
    ~NaiveGridBlocker() override;

    void block() override;
};

/**
 * Octree-based blocking
 */
class OctreeGridBlocker : public GridBlocker {
   public:
    OctreeGridBlocker(const int npoints_ref, double const* x_ref, double const* y_ref, double const* z_ref,
                      double const* w_ref, const int max_points, const int min_points, const double max_radius,
                      std::shared_ptr<BasisExtents> extents);
    ~OctreeGridBlocker() override;

    void block() override;
};
}  // namespace psi
#endif
