#ifndef libmints_gridblock_H
#define libmints_gridblock_H
/*
* gridblock.h
* Definition of class GridBlock for use with numerical integrators 
* (as in KS-DFT) and the various point property calculators
*
* Created by Robert Parrish on 04/15/2010
*/
#include <psi4-dec.h>
using namespace std;

namespace psi { 
/*! \ingroup LIBMINTS */
//! Integration Point/Weight container class (blocks, not individual) 
class GridBlock {
protected:
    /// Weight vector [max_points_]
    double* w_;
    /// x vector [max_points_]
    double* x_;
    /// y vector [max_points_]
    double* y_;
    /// z vector [max_points_]
    double* z_;
    /// Maximum number of points in block at the moment
    int max_points_;
    /// Actual number of valid points
    int true_points_; 
public:
    /** Constructor, allocates block object with max_points size **/
    GridBlock(int max_points);
    /** Destructor, deallocates arrays **/
    ~GridBlock();
    /** Factory Constructor **/
    static GridBlock * createGridBlock(int max_points) {
        return new GridBlock(max_points);
    }
    int getMaxPoints() const {return max_points_; }
    int getTruePoints() const {return true_points_; }
    double* getWeights() const {return w_; }
    double* getX() const {return x_; }
    double* getY() const {return y_; }
    double* getZ() const {return z_; }
    void setGrid(double* x, double* y, double* z, double* w, int n);
    void setTruePoints(int n) { true_points_ = n; }
    void setMaxPoints(int n) { max_points_ = n; }
};
typedef shared_ptr<GridBlock> SharedGridBlock;
}
#endif
