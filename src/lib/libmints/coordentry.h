#ifndef COORDENTRY_H
#define COORDENTRY_H

#include <map>
#include <string>

// Probably but the best, but I'm lazy.
#include "vector3.h"
#include <boost/shared_ptr.hpp>

namespace psi {

class Vector3;

class CoordValue
{
public:
    virtual double compute() =0;
};

class NumberValue : public CoordValue
{
    double value_;
public:
    NumberValue(double value) :value_(value) {}
    double compute() { return value_; }
};

class VariableValue : public CoordValue
{
    const std::string& name_;
    std::map<std::string, double>& geometryVariables_;
    bool negate_;
public:
    VariableValue(const std::string& name, std::map<std::string, double>& geometryVariables, bool negate=false)
        : name_(name), geometryVariables_(geometryVariables), negate_(negate) {}
    double compute() { return negate_ ? -geometryVariables_[name_] : geometryVariables_[name_]; }
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

public:
    CoordEntry(int entry_number, int Z, double charge, double mass, std::string& label);
    virtual ~CoordEntry();

    virtual const Vector3& compute() = 0;
    virtual void set(const Vector3&) = 0;

    bool is_computed() const { return computed_; }
    void invalidate() { computed_ = false; }

    bool is_ghosted() const { return ghosted_; }
    void set_ghosted(bool ghosted) { ghosted_ = ghosted; }

    int Z() const { return ghosted_ ? 0 : Z_; }
    double charge() const { return charge_; }
    double mass() const { return mass_; }
    const std::string& label() const { return label_; }

    int entry_number() const { return entry_number_; }
};

// This version does not support variables in coordinates
class CartesianEntry : public CoordEntry{
public:
    CartesianEntry(int entry_number, int Z, double x, double y, double z, double charge, double mass, std::string& label);

    const Vector3& compute();
    void set(const Vector3&);
};

class ZMatrixEntry : public CoordEntry
{
    boost::shared_ptr<CoordEntry> rto_;
    boost::shared_ptr<CoordValue>   rval_;
    boost::shared_ptr<CoordEntry> ato_;
    boost::shared_ptr<CoordValue>   aval_;
    boost::shared_ptr<CoordEntry> dto_;
    boost::shared_ptr<CoordValue>   dval_;

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
    void set(const Vector3&);
};

}

#endif // COORDENTRY_H
