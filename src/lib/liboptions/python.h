#ifndef _psi_src_lib_liboptions_python_h
#define _psi_src_lib_liboptions_python_h

#include <boost/python.hpp>
#include <boost/python/object.hpp>

namespace psi {

class PythonDataType : public DataType
{
    boost::python::object python_object_;
public:
    PythonDataType();
    PythonDataType(const boost::python::object& p);
    virtual ~PythonDataType();

    virtual std::string type() const;

    const boost::python::object& to_python() const;

    void assign(const boost::python::object& p);
};

}

#endif // _psi_src_lib_liboptions_python_h
