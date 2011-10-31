#ifndef libmints_typedefs_h
#define libmints_typedefs_h

// Forward declare boost
namespace boost {
template <class>
class shared_ptr;
}

// Forward declare psi
namespace psi {
class Matrix;
}

typedef boost::shared_ptr<psi::Matrix> SharedMatrix;

#endif // libmints_typedefs_h
