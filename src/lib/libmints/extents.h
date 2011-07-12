#ifndef _psi_src_lib_libmints_extents_h_
#define _psi_src_lib_libmints_extents_h_

#include <cstdio>
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace psi {

class BasisSet;
class Vector;
class Matrix;

/*! \ingroup MINTS */

//! Basis set extents calaculator
class Extents {

protected:
    /// Center 1 basis functions
    boost::shared_ptr<BasisSet> basis1_;
    /// Center 2 basis functions
    boost::shared_ptr<BasisSet> basis2_;

public:
    Extents(boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2);
    virtual ~Extents();

    /// Radii beyond which basis1 functions are strictly < epsilon
    boost::shared_ptr<Vector> extents1(double epsilon);
    /// Convenience routine for x center of basis1 
    boost::shared_ptr<Vector> x1();    
    /// Convenience routine for y center of basis1 
    boost::shared_ptr<Vector> y1();    
    /// Convenience routine for z center of basis1 
    boost::shared_ptr<Vector> z1();    

    /// Radii beyond which basis2 functions are strictly < epsilon
    boost::shared_ptr<Vector> extents2(double epsilon);
    /// Convenience routine for x center of basis2 
    boost::shared_ptr<Vector> x2();    
    /// Convenience routine for y center of basis2 
    boost::shared_ptr<Vector> y2();    
    /// Convenience routine for z center of basis2 
    boost::shared_ptr<Vector> z2();    

    /// Radii beyond which the product of basis1 x basis2 functions are strictly < epsilon
    boost::shared_ptr<Matrix> extents12(double epsilon);
    /// x center of basis1 x basis2 functions 
    boost::shared_ptr<Matrix> x12();    
    /// y center of basis1 x basis2 functions 
    boost::shared_ptr<Matrix> y12();    
    /// z center of basis1 x basis2 functions 
    boost::shared_ptr<Matrix> z12();    

};

}
#endif
