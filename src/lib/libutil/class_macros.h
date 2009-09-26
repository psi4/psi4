#ifndef _psi_src_lib_libutil_class_macros_h_
#define _psi_src_lib_libutil_class_macros_h_

// This macro may be used in a class to declare
// a private variable with a set and get function
#define READ_WRITE_VAR(type,name) \
public: \
type get_name() const {return name;} \
void  set_name(type _value_) {name = _value_;}  \
private: \
type name;

// This macro may be used in a class to declare
// a private variable with a get function
#define READ_VAR(type,name) \
public: \
type get_name() const {return name;} \
private: \
type name;

#endif // _psi_src_lib_libutil_class_macros_h_
