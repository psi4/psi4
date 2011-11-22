#include "python.h"

#include <boost/python.hpp>

namespace psi {

PythonDataType::PythonDataType()
{
}

PythonDataType::PythonDataType(const boost::python::object &p)
    : python_object_(p)
{
}

PythonDataType::~PythonDataType()
{ }

std::string PythonDataType::type() const{
    return std::string("python");
}

void PythonDataType::assign(const boost::python::object &p)
{
    python_object_ = p;
    changed();
}

const boost::python::object& PythonDataType::to_python() const
{
    return python_object_;
}

}
