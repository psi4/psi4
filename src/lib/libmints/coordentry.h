#ifndef COORDENTRY_H
#define COORDENTRY_H

#include <map>
#include <string>
#include "vector3.h"
#include <boost/shared_ptr.hpp>

namespace psi {

class Vector3;


/**
 * An abstract class to handle storage of Cartesian coordinate values, which
 * may be defined in terms of other variables through this mechanism, greatly
 * simplifying Z-matrix specification, for example.
 */
class CoordValue
{
public:
    /// Computes the current value, and returns it.
    virtual double compute() =0;
    /// Sets the current value.
    virtual void set(double val) =0;
};

/**
 * Specialization of CoordValue that is simply a number to be stored.
 */
class NumberValue : public CoordValue
{
    double value_;
public:
    NumberValue(double value) :value_(value) {}
    double compute() { return value_; }
    void set(double val) { value_ = val; }
};

/**
 * Specialization of CoordValue, where the current value depends on the list of
 * geometry values stored by the molecule.
 */
class VariableValue : public CoordValue
{
    const std::string name_;
    std::map<std::string, double>& geometryVariables_;
    bool negate_;
public:
    VariableValue(const std::string name, std::map<std::string, double>& geometryVariables, bool negate=false)
        : name_(name), geometryVariables_(geometryVariables), negate_(negate) {}
    double compute() { return negate_ ? -geometryVariables_[name_] : geometryVariables_[name_]; }
    void set(double val) { geometryVariables_[name_] = negate_ ? -val : val; }
};

//class EquationValue
//{

//};

class CoordEntry
{
protected:
    int entry_number_;
    bool computed_;
    Vector3 coordinates_;

    int Z_;
    double charge_;
    double mass_;
    std::string label_;
    bool ghosted_;
    /// Computes the distance between two sets of coordinates
    static double r(const Vector3 &a1, const Vector3 &a2) { return a1.distance(a2); }
    /// Computes the angle (in rad.) between three sets of coordinates.
    static double a(const Vector3 &a1, const Vector3 &a2, const Vector3 &a3)
        { Vector3 eBA(a2-a1), eBC(a2-a3); eBA.normalize(); eBC.normalize(); return acos(eBA.dot(eBC)); }
    /// Computes the dihedral (in rad.) between four sets of coordinates.
    static double d(const Vector3 &a1, const Vector3 &a2, const Vector3 &a3, const Vector3 &a4)
        { Vector3 eBA(a1-a2), eBC(a3-a2), eCD(a4-a3); eBA.normalize(); eBC.normalize();
          eCD.normalize(); return acos(eBA.dot(eBC)*eBC.dot(eCD)/(sin(a(a1,a2,a3)*sin(a(a2,a3,a4)))));}


public:
    CoordEntry(int entry_number, int Z, double charge, double mass, std::string& label);
    virtual ~CoordEntry();

    /// Computes the values of the coordinates (in whichever units were inputted), returning them in a Vector.
    virtual const Vector3& compute() = 0;
    /// Given the current set of coordinates, updates the values of this atom's coordinates,
    /// and any variables that may depend on it.
    virtual void set_coordinates(double x, double y, double z) = 0;

    /// Whether the current atom's coordinates are up-to-date.
    bool is_computed() const { return computed_; }
    /// Flags the current coordinates as being outdated.
    void invalidate() { computed_ = false; }

    /// Whether the current atom is ghosted or not.
    bool is_ghosted() const { return ghosted_; }
    /// Flag the atom as either ghost or real.
    void set_ghosted(bool ghosted) { ghosted_ = ghosted; }

    /// The nuclear charge of the current atom (0 if ghosted).
    int Z() const { return ghosted_ ? 0 : Z_; }
    /// The "atomic charge" of the current atom (for SAD purposes).
    double charge() const { return charge_; }
    /// The atomic mass of the current atom.
    double mass() const { return mass_; }
    /// The atomic symbol.
    const std::string& label() const { return label_; }
    /// The order in which this appears in the full atom list.
    int entry_number() const { return entry_number_; }
};

class CartesianEntry : public CoordEntry{
    boost::shared_ptr<CoordValue> x_;
    boost::shared_ptr<CoordValue> y_;
    boost::shared_ptr<CoordValue> z_;
public:
    CartesianEntry(int entry_number, int Z, double charge, double mass, std::string& label,
                   boost::shared_ptr<CoordValue>(x), boost::shared_ptr<CoordValue>(y), boost::shared_ptr<CoordValue>(z));

    const Vector3& compute();
    void set_coordinates(double x, double y, double z);
    
};

class ZMatrixEntry : public CoordEntry
{
    boost::shared_ptr<CoordEntry> rto_;
    boost::shared_ptr<CoordValue> rval_;
    boost::shared_ptr<CoordEntry> ato_;
    boost::shared_ptr<CoordValue> aval_;
    boost::shared_ptr<CoordEntry> dto_;
    boost::shared_ptr<CoordValue> dval_;

public:
    ZMatrixEntry(int entry_number, int Z, double charge, double mass, std::string& label,
                 boost::shared_ptr<CoordEntry> rto=boost::shared_ptr<CoordEntry>(),
                 boost::shared_ptr<CoordValue> rval=boost::shared_ptr<CoordValue>(),
                 boost::shared_ptr<CoordEntry> ato=boost::shared_ptr<CoordEntry>(),
                 boost::shared_ptr<CoordValue> aval=boost::shared_ptr<CoordValue>(),
                 boost::shared_ptr<CoordEntry> dto=boost::shared_ptr<CoordEntry>(),
                 boost::shared_ptr<CoordValue> dval=boost::shared_ptr<CoordValue>());
    virtual ~ZMatrixEntry();

    const Vector3& compute();
    void set_coordinates(double x, double y, double z);
};

}

#endif // COORDENTRY_H
