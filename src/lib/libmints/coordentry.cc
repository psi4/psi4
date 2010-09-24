#include <cmath>
#include <libmints/vector3.h>
#include <libmints/molecule.h>
#include <exception.h>
#include "coordentry.h"

using namespace psi;

CoordEntry::CoordEntry(int entry_number, int Z, double charge, double mass, std::string& label)
    : entry_number_(entry_number), computed_(false), Z_(Z),
      charge_(charge), mass_(mass), label_(label), ghosted_(false)
{
}

ZMatrixEntry::ZMatrixEntry(int entry_number, int Z, double charge, double mass, std::string& label,
                           ZMatrixEntry *rto, CoordValue *rval,
                           ZMatrixEntry *ato, CoordValue *aval,
                           ZMatrixEntry *dto, CoordValue *dval)
    : CoordEntry(entry_number, Z, charge, mass, label),
      rto_(rto), rval_(rval),
      ato_(ato), aval_(aval),
      dto_(dto), dval_(dval)
{

}

ZMatrixEntry::~ZMatrixEntry()
{
    if (rval_)
        delete rval_;
    if (aval_)
        delete aval_;
    if (dval_)
        delete dval_;
}

const Vector3& ZMatrixEntry::compute()
{
    if (computed_)
        return coordinates_;

    if(rto_ == 0 && ato_ == 0 && dto_ == 0){
        coordinates_[0] = 0.0;
        coordinates_[1] = 0.0;
        coordinates_[2] = 0.0;
    }else if(ato_ == 0 && dto_ == 0){
        coordinates_[0] = 0.0;
        coordinates_[1] = 0.0;
        coordinates_[2] = rval_->compute();
    }else if(dto_ == 0){
        double a = aval_->compute() * M_PI / 180.0;
        double r = rval_->compute();
        if(rto_->entry_number_ == 0){
            coordinates_[0] = r*sin(a);
            coordinates_[1] = 0.0;
            coordinates_[2] = r*cos(a);
        }else{
            coordinates_[0] = r*sin(a);
            coordinates_[1] = 0.0;
            coordinates_[2] = rto_->compute()[2] - r*cos(a);
        }
    }else{
        double r = rval_->compute();
        double a = aval_->compute();
        double d = dval_->compute();

        a *= M_PI / 180.0;
        d *= M_PI / 180.0;

        /*
         * The atom specification is
         *      this       rTo   rVal  aTo  aVal   dTo   dVal
         *        D         C           B           A
         * which allows us to define the vector from B->C (eBC) as the +z axis, and eAB
         * lies in the xz plane.  Then eX, eY and eZ (=eBC) are the x, y, and z axes, respecively.
         */
        const Vector3& rTo = rto_->compute();
        const Vector3& aTo = ato_->compute();
        const Vector3& dTo = dto_->compute();

        Vector3 eAB = aTo - dTo;
        Vector3 eBC = rTo - aTo;
        eAB.normalize();
        eBC.normalize();
        double cosBCD = cos(a);
        double sinBCD = sin(a);
        double cosABCD = cos(d);
        double sinABCD = sin(d);
        double cosABC = -eAB.dot(eBC);
        double sinABC = sqrt(1.0 - cosABC * cosABC);
        if(sinABC < LINEAR_A_TOL || sinBCD < LINEAR_A_TOL)
            throw PSIEXCEPTION("Atom defines a dihedral using a linear bond angle");
        Vector3 eY = eAB.perp_unit(eBC);
        Vector3 eX = eY.perp_unit(eBC);

        coordinates_[0] = rTo[0] + r * ( -eBC[0] * cosBCD + eX[0] * sinBCD * cosABCD + eY[0] * sinBCD * sinABCD);
        coordinates_[1] = rTo[1] + r * ( -eBC[1] * cosBCD + eX[1] * sinBCD * cosABCD + eY[1] * sinBCD * sinABCD);
        coordinates_[2] = rTo[2] + r * ( -eBC[2] * cosBCD + eX[2] * sinBCD * cosABCD + eY[2] * sinBCD * sinABCD);
    }

    computed_ = true;
    return coordinates_;
}
