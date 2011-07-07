/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/
#ifndef MADNESS_MRA_FUNCDEFAULTS_H__INCLUDED
#define MADNESS_MRA_FUNCDEFAULTS_H__INCLUDED

/// \file funcdefaults.h
/// \brief Provides FunctionDefaults and utilities for coordinate transformation
/// \ingroup mrabcext

#include <world/world.h>
#include <world/array.h>
#include <world/worlddc.h>
#include <tensor/tensor.h>
#include <mra/key.h>

namespace madness {
    template <typename T, std::size_t NDIM> class FunctionImpl;

    /// The maximum wavelet order presently supported
    static const int MAXK = 30;

    /// The maximum depth of refinement possible
    static const int MAXLEVEL = 8*sizeof(Translation)-2;

    enum BCType {BC_ZERO, BC_PERIODIC, BC_FREE, BC_DIRICHLET, BC_ZERONEUMANN, BC_NEUMANN};

    /*!
      \brief This class is used to specify boundary conditions for all operators
      \ingroup mrabcext

      Exterior boundary conditions (i.e., on the simulation domain)
      are associated with operators (not functions).  The types of
      boundary conditions available are in the enum BCType.

      The default boundary conditions are obtained from the FunctionDefaults.
      For non-zero Dirichlet and Neumann conditions additional information
      must be provided when derivative operators are constructed. For integral
      operators, only periodic and free space are supported.
    */
    template<std::size_t NDIM>
    class BoundaryConditions {
    private:
        // Used to use STL vector but static data on  a MAC was
        // causing problems.
        BCType bc[NDIM*2];

    public:
        /// Constructor. Default boundary condition set to free space
        BoundaryConditions(BCType code=BC_FREE)
        {
            for (std::size_t i=0; i<NDIM*2; ++i) bc[i] = code;
        }

        /// Copy constructor is deep
        BoundaryConditions(const BoundaryConditions<NDIM>& other)
        {
            *this = other;
        }

        /// Assignment makes deep copy
        BoundaryConditions<NDIM>&
        operator=(const BoundaryConditions<NDIM>& other) {
            if (&other != this) {
                for (std::size_t i=0; i<NDIM*2; ++i) bc[i] = other.bc[i];
            }
            return *this;
        }

        /// Returns value of boundary condition

        /// @param d Dimension (0,...,NDIM-1) for boundary condition
        /// @param i Side (0=left, 1=right) for boundary condition
        /// @return Value of boundary condition
        BCType operator()(std::size_t d, int i) const {
            MADNESS_ASSERT(d<NDIM && i>=0 && i<2);
            return bc[2*d+i];
        }

        /// Returns reference to boundary condition

        /// @param d Dimension (0,...,NDIM-1) for boundary condition
        /// @param i Side (0=left, 1=right) for boundary condition
        /// @return Value of boundary condition
        BCType& operator()(std::size_t d, int i) {
            MADNESS_ASSERT(d<NDIM && i>=0 && i<2);
            return bc[2*d+i];
        }

        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & bc;
        }

        /// Translates code into human readable string

        /// @param code Code for boundary condition
        /// @return String describing boundary condition code
        static const char* code_as_string(BCType code) {
            static const char* codes[] = {"zero","periodic","free","Dirichlet","zero Neumann","Neumann"};
            return codes[code];
        }

        /// Convenience for application of integral operators

        /// @return Returns a vector indicating if each dimension is periodic
        std::vector<bool> is_periodic() const {
            std::vector<bool> v(NDIM);
            for (std::size_t d=0; d<NDIM; ++d)
                v[d] = (bc[2*d]==BC_PERIODIC);
            return v;
        }
    };


    template <std::size_t NDIM>
    static
    inline
    std::ostream& operator << (std::ostream& s, const BoundaryConditions<NDIM>& bc) {
        s << "BoundaryConditions(";
        for (int d=0; d<NDIM; ++d) {
            s << bc.code_as_string(bc(d,0)) << ":" << bc.code_as_string(bc(d,1));
            if (d == NDIM-1)
                s << ")";
            else
                s << ", ";
        }
        return s;
    }

    /// FunctionDefaults holds default paramaters as static class members

    /// Declared and initialized in mra.cc and/or funcimpl::initialize.
    ///
    /// Currently all functions of the same dimension share the same cell dimensions
    /// since they are stored inside FunctionDefaults ... if you change the
    /// cell dimensions *all* functions of that dimension are affected.
    ///
    /// N.B.  Ultimately, we may need to make these defaults specific to each
    /// world, as should be all global state.
    /// \ingroup mra
    template <std::size_t NDIM>
    class FunctionDefaults {
    private:
        static int k;                 ///< Wavelet order
        static double thresh;          ///< Truncation threshold
        static int initial_level;      ///< Initial level for fine scale projection
        static int max_refine_level;   ///< Level at which to stop refinement
        static int truncate_mode;    ///< Truncation method
        static bool refine;            ///< Whether to refine new functions
        static bool autorefine;        ///< Whether to autorefine in multiplication, etc.
        static bool debug;             ///< Controls output of debug info
        static bool truncate_on_project; ///< If true initial projection inserts at n-1 not n
        static bool apply_randomize;   ///< If true use randomization for load balancing in apply integral operator
        static bool project_randomize; ///< If true use randomization for load balancing in project/refine
        static BoundaryConditions<NDIM> bc; ///< Default boundary conditions
        static Tensor<double> cell ;   ///< cell[NDIM][2] Simulation cell, cell(0,0)=xlo, cell(0,1)=xhi, ...
        static Tensor<double> cell_width;///< Width of simulation cell in each dimension
        static Tensor<double> rcell_width; ///< Reciprocal of width
        static double cell_volume;      ///< Volume of simulation cell
        static double cell_min_width;   ///< Size of smallest dimension
        static std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > > pmap; ///< Default mapping of keys to processes

        static void recompute_cell_info() {
            MADNESS_ASSERT(cell.dim(0)==NDIM && cell.dim(1)==2 && cell.ndim()==2);
            cell_width = cell(_,1)-cell(_,0);
            cell_volume = cell_width.product();
            cell_min_width = cell_width.min();
            rcell_width = copy(cell_width);
            for (std::size_t i=0; i<NDIM; ++i) rcell_width(i) = 1.0/rcell_width(i);
        }

    public:

        /// Used to set defaults to k=7, thresh=1-5, for a unit cube [0,1].
        static void set_defaults(World& world);

        /// Returns the default wavelet order
        static int get_k() {
            return k;
        }

        /// Sets the default wavelet order

        /// Existing functions are unaffacted.
        static void set_k(int value) {
            k=value;
            MADNESS_ASSERT(k>0 && k<=MAXK);
        }

        /// Returns the default threshold
        static double get_thresh() {
            return thresh;
        }

        /// Sets the default threshold

        /// Existing functions are unaffected
        static void set_thresh(double value) {
            thresh=value;
        }

        /// Returns the default initial projection level
        static int get_initial_level() {
            return initial_level;
        }

        /// Sets the default initial projection level

        /// Existing functions are unaffected
        static void set_initial_level(int value) {
            initial_level=value;
            MADNESS_ASSERT(value>0 && value<MAXLEVEL);
        }

        /// Gets the default maximum adaptive refinement level
        static int get_max_refine_level() {
            return max_refine_level;
        }

        /// Sets the default maximum adaptive refinement level

        /// Existing functions are unaffected
        static void set_max_refine_level(int value) {
            max_refine_level=value;
            MADNESS_ASSERT(value>0 && value<MAXLEVEL);
        }

        /// Gets the default truncation mode
        static int get_truncate_mode() {
            return truncate_mode;
        }

        /// Sets the default truncation mode

        /// Existing functions are unaffected
        static void set_truncate_mode(int value) {
            truncate_mode=value;
            MADNESS_ASSERT(value>=0 && value<3);
        }

        /// Gets the default adaptive refinement flag
        static bool get_refine() {
            return refine;
        }

        /// Sets the default adaptive refinement flag

        /// Existing functions are unaffected
        static void set_refine(bool value) {
            refine=value;
        }

        /// Gets the default adaptive autorefinement flag
        static bool get_autorefine() {
            return autorefine;
        }

        /// Sets the default adaptive autorefinement flag

        /// Existing functions are unaffected
        static void set_autorefine(bool value) {
            autorefine=value;
        }

        /// Gets the default debug flag (is this used anymore?)
        static bool get_debug() {
            return debug;
        }

        /// Sets the default debug flag (is this used anymore?)

        /// Not sure if this does anything useful
        static void set_debug(bool value) {
            debug=value;
        }

        /// Gets the default truncate on project flag
        static bool get_truncate_on_project() {
            return truncate_on_project;
        }

        /// Sets the default truncate on project flag

        /// Existing functions are unaffected
        static void set_truncate_on_project(bool value) {
            truncate_on_project=value;
        }

        /// Gets the random load balancing for integral operators flag
        static bool get_apply_randomize() {
            return apply_randomize;
        }

        /// Sets the random load balancing for integral operators flag
        static void set_apply_randomize(bool value) {
            apply_randomize=value;
        }


        /// Gets the random load balancing for projection flag
        static bool get_project_randomize() {
            return project_randomize;
        }

        /// Sets the random load balancing for projection flag
        static void set_project_randomize(bool value) {
            project_randomize=value;
        }

        /// Returns the default boundary conditions
        static const BoundaryConditions<NDIM>& get_bc() {
            return bc;
        }

        /// Sets the default boundary conditions
        static void set_bc(const BoundaryConditions<NDIM>& value) {
            bc=value;
        }

        /// Gets the user cell for the simulation
        static const Tensor<double>& get_cell() {
            return cell;
        }

        /// Gets the user cell for the simulation

        /// Existing functions are probably rendered useless
        static void set_cell(const Tensor<double>& value) {
            cell=copy(value);
            recompute_cell_info();
        }

        /// Sets the user cell to be cubic with each dimension having range \c [lo,hi]

        /// Existing functions are probably rendered useless
        static void set_cubic_cell(double lo, double hi) {
            cell(_,0)=lo;
            cell(_,1)=hi;
            recompute_cell_info();
        }

        /// Returns the width of each user cell dimension
        static const Tensor<double>& get_cell_width() {
            return cell_width;
        }

        /// Returns the reciprocal of the width of each user cell dimension
        static const Tensor<double>& get_rcell_width() {
            return rcell_width;
        }

        /// Returns the minimum width of any user cell dimension
        static double get_cell_min_width() {
            return cell_min_width;
        }

        /// Returns the volume of the user cell
        static double get_cell_volume() {
            return cell_volume;
        }

        /// Returns the default process map
        static std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& get_pmap() {
            return pmap;
        }

        /// Sets the default process map (does \em not redistribute existing functions)

        /// Existing functions are probably rendered useless
        static void set_pmap(const std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& value) {
            pmap = value;
        }

        /// Sets the default process map and redistributes all functions using the old map
        static void redistribute(World& world, const std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& newpmap) {
            pmap->redistribute(world,newpmap);
            pmap = newpmap;
        }

    };


    /// Convert user coords (cell[][]) to simulation coords ([0,1]^ndim)
    template <std::size_t NDIM>
    static inline void user_to_sim(const Vector<double,NDIM>& xuser, Vector<double,NDIM>& xsim) {
        for (std::size_t i=0; i<NDIM; ++i)
            xsim[i] = (xuser[i] - FunctionDefaults<NDIM>::get_cell()(i,0)) * FunctionDefaults<NDIM>::get_rcell_width()[i];
    }


    /// Convert simulation coords ([0,1]^ndim) to user coords (FunctionDefaults<NDIM>::get_cell())
    template <std::size_t NDIM>
    static void sim_to_user(const Vector<double,NDIM>& xsim, Vector<double,NDIM>& xuser) {
        for (std::size_t i=0; i<NDIM; ++i)
            xuser[i] = xsim[i]*FunctionDefaults<NDIM>::get_cell_width()[i] + FunctionDefaults<NDIM>::get_cell()(i,0);
    }


}
#endif // MADNESS_MRA_FUNCDEFAULTS_H__INCLUDED
