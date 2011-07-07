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
  
  $Id: sdf_shape_2D.h 1888 2010-06-02 23:03:53Z mgr522 $
*/

/**
  \file mra/sdf_shape_2D.h
  \brief Implements the SignedDFInterface for common 2-D geometric objects.
  \ingroup mrabcint

  This file provides signed distance functions for common 2-D geometric objects:
  - Circle
  - Rectangle
  
  \note The signed distance functions should be the shortest distance between
  a point and \b any point on the surface.  This is hard to calculate in many
  cases, so we use contours here.  The surface layer may not be equally thick
  around all points on the surface.  Some shapes have the exact shortest
  distances, which will be noted on a class-by-class basis.  Serious usage
  of the contour-based signed distance functions is not recommended.
*/  
  
#ifndef MADNESS_MRA_SDF_SHAPE_2D_H__INCLUDED
#define MADNESS_MRA_SDF_SHAPE_2D_H__INCLUDED

#include <mra/sdf_domainmask.h>

namespace madness {

    /// \brief A circle (2 dimensions)
    class SDFCircle : public SignedDFInterface<2> {
    protected:
        const double radius; ///< Radius of circle
        const coord_2d center; ///< Center of circle

    public:
        /** \brief SDF for a sphere

            \param radius The radius of the sphere
            \param center The center of the sphere */
        SDFCircle(const double radius, const coord_2d &center) 
            : radius(radius)
            , center(center) 
        {}

        /** \brief Computes the normal distance

            This SDF is exact, and easy to show.

            \param pt Point at which to compute the distance from the surface
            \return The signed distance from the surface */
        double sdf(const coord_2d& pt) const {
            double temp, r;
            int i;
            
            r = 0.0;
            for(i = 0; i < 2; ++i) {
                temp = pt[i] - center[i];
                r += temp * temp;
            }
            
            return sqrt(r) - radius;
        }

        /** \brief Computes the gradient of the SDF.

            \param pt Point at which to compute the gradient
            \return the gradient */
        coord_2d grad_sdf(const coord_2d& pt) const {
            double x = pt[0] - center[0];
            double y = pt[1] - center[1];
            double r = sqrt(x*x + y*y);
            coord_2d g;
            if(r < 1.0e-6) {
                g[0] = g[1] = 0.0;
            }
            else {
                g[0] = x/r;
                g[1] = y/r;
            }
            return g;
        }
    };

    /** \brief A rectangle (2 dimensions)

        This SDF naively uses contours, and should be improved for serious
        usage.

        \note LIMIT -- the 2 primary axes must be x and y */
    class SDFRectangle : public SignedDFInterface<2> {
    protected:
        const coord_2d lengths;  ///< Half the length of each side
        const coord_2d center; ///< the center

    public:
        /** \brief Constructor for rectangle

            \param length The lengths of the rectangle
            \param center The center of the rectangle */
        SDFRectangle(const coord_2d& length, const coord_2d& center) 
            : lengths(length*0.5), center(center) 
        {}

        /** \brief Computes the normal distance

            \param pt Point at which to compute the distance from the surface
            \return The signed distance from the surface */
        double sdf(const coord_2d& pt) const {
            double diff, max;
            
            max = fabs(pt[0] - center[0]) - lengths[0];
            diff = fabs(pt[1] - center[1]) - lengths[1];
            if(diff > max)
                max = diff;
            
            return max;
        }

        /** Computes the gradient of the SDF.

            \param pt Point at which to compute the gradient
            \return the gradient */
        coord_2d grad_sdf(const coord_2d& pt) const {
            MADNESS_EXCEPTION("gradient method is not yet implemented for this shape",0);
        }
    };

} // end of madness namespace

#endif // MADNESS_MRA_SDF_SHAPE_2D_H__INCLUDED
