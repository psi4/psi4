#include <cmath>
#include <libmints/vector3.h>
#include <libmints/molecule.h>
#include <exception.h>
#include "coordentry.h"

using namespace psi;

double
VariableValue::compute()
{
    if(geometryVariables_.count(name_) == 0)
        throw PSIEXCEPTION("Variable " + name_ + " used in geometry specification has not been defined");
    return negate_ ? -geometryVariables_[name_] : geometryVariables_[name_];
}

CoordEntry::CoordEntry(int entry_number, int Z, double charge, double mass, const std::string& symbol, const std::string& label)
    : entry_number_(entry_number), computed_(false), Z_(Z),
      charge_(charge), mass_(mass), symbol_(symbol), label_(label), ghosted_(false)
{
}

CoordEntry::CoordEntry(int entry_number, int Z, double charge, double mass, const std::string& symbol,
           const std::string& label, const std::map<std::string, std::string>& basis)
    : entry_number_(entry_number), computed_(false), Z_(Z),
      charge_(charge), mass_(mass), symbol_(symbol), label_(label), ghosted_(false),
      basissets_(basis)
{

}

CoordEntry::~CoordEntry()
{

}

void CoordEntry::set_basisset(const std::string& name, const std::string& type)
{
    basissets_[type] = name;
}

const std::string& CoordEntry::basisset(const std::string& type) const
{
    std::map<std::string, std::string>::const_iterator iter = basissets_.find(type);

    if (iter == basissets_.end())
        throw PSIEXCEPTION("CoordEntry::basisset: Basisset not set for "+label_+" and type of " + type);

    return (*iter).second;
}


CartesianEntry::CartesianEntry(int entry_number, int Z, double charge, double mass, const std::string& symbol, const std::string& label,
                               boost::shared_ptr<CoordValue> x, boost::shared_ptr<CoordValue> y, boost::shared_ptr<CoordValue> z)
    : CoordEntry(entry_number, Z, charge, mass, symbol, label), x_(x), y_(y), z_(z)
{
}

CartesianEntry::CartesianEntry(int entry_number, int Z, double charge, double mass, const std::string& symbol, const std::string& label,
               boost::shared_ptr<CoordValue> x, boost::shared_ptr<CoordValue> y, boost::shared_ptr<CoordValue> z,
               const std::map<std::string, std::string>& basis)
    : CoordEntry(entry_number, Z, charge, mass, symbol, label, basis), x_(x), y_(y), z_(z)
{
}

const Vector3& CartesianEntry::compute()
{
    if(computed_)
        return coordinates_;

    coordinates_[0] = x_->compute();
    coordinates_[1] = y_->compute();
    coordinates_[2] = z_->compute();
    computed_ = true;
    return coordinates_;
}

void
CartesianEntry::set_coordinates(double x, double y, double z)
{
    coordinates_[0] = x;
    coordinates_[1] = y;
    coordinates_[2] = z;

    x_->set(x);
    y_->set(y);
    z_->set(z);

    computed_ = true;
}

ZMatrixEntry::ZMatrixEntry(int entry_number, int Z, double charge, double mass, const std::string& symbol, const std::string& label,
                           boost::shared_ptr<CoordEntry> rto, boost::shared_ptr<CoordValue> rval,
                           boost::shared_ptr<CoordEntry> ato, boost::shared_ptr<CoordValue> aval,
                           boost::shared_ptr<CoordEntry> dto, boost::shared_ptr<CoordValue> dval)
    : CoordEntry(entry_number, Z, charge, mass, symbol, label),
      rto_(rto), rval_(rval),
      ato_(ato), aval_(aval),
      dto_(dto), dval_(dval)
{
}

ZMatrixEntry::ZMatrixEntry(int entry_number, int Z, double charge, double mass, const std::string& symbol, const std::string& label,
             const std::map<std::string, std::string>& basis,
             boost::shared_ptr<CoordEntry> rto,
             boost::shared_ptr<CoordValue> rval,
             boost::shared_ptr<CoordEntry> ato,
             boost::shared_ptr<CoordValue> aval,
             boost::shared_ptr<CoordEntry> dto,
             boost::shared_ptr<CoordValue> dval)
    : CoordEntry(entry_number, Z, charge, mass, symbol, label, basis),
      rto_(rto), rval_(rval),
      ato_(ato), aval_(aval),
      dto_(dto), dval_(dval)
{
}

ZMatrixEntry::~ZMatrixEntry()
{
}

void
ZMatrixEntry::set_coordinates(double x, double y, double z)
{
    coordinates_[0] = x;
    coordinates_[1] = y;
    coordinates_[2] = z;

    if(rto_ != 0){
        if(!rto_->is_computed())
             throw PSIEXCEPTION("Coordinates have been set in the wrong order");
        rval_->set(r(coordinates_, rto_->compute()));
    }
    if(ato_ != 0){
        if(!ato_->is_computed())
             throw PSIEXCEPTION("Coordinates have been set in the wrong order");
        aval_->set(180.0*a(coordinates_, rto_->compute(), ato_->compute())/M_PI);
    }
    if(dto_ != 0){
        if(!dto_->is_computed())
             throw PSIEXCEPTION("Coordinates have been set in the wrong order");
        dval_->set(180.0*d(coordinates_, rto_->compute(), ato_->compute(), dto_->compute())/M_PI);
    }

    computed_ = true;
}


/**
 * Computes the coordinates of the current atom's entry
 * @Return The Cartesian Coordinates, in Bohr
 */
const Vector3& ZMatrixEntry::compute()
{
    if (computed_){
        return coordinates_;
}

    if(rto_ == 0 && ato_ == 0 && dto_ == 0){
        coordinates_[0] = 0.0;
        coordinates_[1] = 0.0;
        coordinates_[2] = 0.0;
    }else if(ato_ == 0 && dto_ == 0){
        coordinates_[0] = 0.0;
        coordinates_[1] = 0.0;
        coordinates_[2] = rval_->compute();
    }else if(dto_ == 0){
        double r = rval_->compute();
        double a = aval_->compute() * M_PI/180.0;
        if(rto_->entry_number() == 0){
            coordinates_[0] = r*sin(a);
            coordinates_[1] = 0.0;
            coordinates_[2] = r*cos(a);
        }else{
            coordinates_[0] = r*sin(a);
            coordinates_[1] = 0.0;
            coordinates_[2] = rto_->compute()[2] - r*cos(a);
        }
//TODO generalize this, for general orientations of the first two atoms - this is all that's missing for any zmat/cart combination in geometries
    }else{
        double r = rval_->compute();
        double a = aval_->compute() * M_PI/180.0;
        double d = dval_->compute() * M_PI/180.0;

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

