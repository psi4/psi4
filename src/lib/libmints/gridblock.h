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

#ifndef libmints_gridblock_H
#define libmints_gridblock_H
/*
* gridblock.h
* Definition of class GridBlock for use with numerical integrators 
* (as in KS-DFT) and the various point property calculators
*
* Created by Robert Parrish on 04/15/2010
*/
#include <boost/shared_ptr.hpp>

namespace psi { 
/*! \ingroup LIBMINTS */
//! Integration Point/Weight container class (blocks, not individual) 
class GridBlock {
public:
    /// Weight vector 
    double* w_;
    /// x vector 
    double* x_;
    /// y vector 
    double* y_;
    /// z vector 
    double* z_;
    /// Maximum number of points in block at the moment
    int max_points_;
    /// Actual number of valid points
    int true_points_; 
    GridBlock() {}
    ~GridBlock() {}

    int getMaxPoints() const {return max_points_; }
    int getTruePoints() const {return true_points_; }
    double* getWeights() const {return w_; }
    double* getX() const {return x_; }
    double* getY() const {return y_; }
    double* getZ() const {return z_; }

    void setGrid(double* x, double* y, double* z, double* w) {
        x_ = x;
        y_ = y;
        z_ = z;
        w_ = w;
    }

    void setTruePoints(int n) { true_points_ = n; }
    void setMaxPoints(int n) { max_points_ = n; }
};
typedef boost::shared_ptr<GridBlock> SharedGridBlock;
}
#endif
