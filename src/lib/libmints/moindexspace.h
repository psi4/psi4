#ifndef _psi_src_lib_libmints_moindexspace_h_
#define _psi_src_lib_libmints_moindexspace_h_

#include <string>

namespace boost {
template<class T>
class shared_ptr;
}

namespace psi {

class Matrix;
class IntegralFactory;
class Vector;
class Dimension;
class BasisSet;

class MOIndexSpace
{
    /// Number of irreps.
    int nirrep_;

    /// Name of the orbital space.
    std::string name_;

    /// AO->MO transformation matrix (ao x mo)
    boost::shared_ptr<Matrix> C_;

    /// MO "eigenvalues"
    boost::shared_ptr<Vector> evals_;

    /// AO basis set
    boost::shared_ptr<BasisSet> basis_;
    /// Integral factory that as
    boost::shared_ptr<IntegralFactory> ints_;

    /// MO Dimensionality
    Dimension dim_; // dim_.n() better equal nirrep_

    MOIndexSpace();
public:
    MOIndexSpace(const std::string& name,
                 const boost::shared_ptr<Matrix>& full_C,   // Should this be 4 C's instead?
                 const boost::shared_ptr<Vector>& evals,
                 const boost::shared_ptr<BasisSet>& basis,
                 const boost::shared_ptr<IntegralFactory>& ints);

    MOIndexSpace(const std::string& name,
                 const boost::shared_ptr<Wavefunction>& wave);

    int nirrep() const;
    const std::string& name() const;
    const boost::shared_ptr<Matrix>& C() const;
    const boost::shared_ptr<Vector>& evals() const;
    const boost::shared_ptr<BasisSet>& basis() const;
    const boost::shared_ptr<IntegralFactory>& integral() const;
    const Dimension& dim() const;

    /** Returns the matrix that transforms space1 to space2. Throws if unable to.
        The transform can only be constructed if the overlap of space1.basis()
        and space2.basis() is nonzero.

        The returned matrix has dimensions of space2.C().coldim() and space1.C().coldim().
      */
    static boost::shared_ptr<Matrix> transform(const MOIndexSpace& space2, const MOIndexSpace& space1);

    /** Returns the overlap matrix between space2 and space1.
        The matrix has dimensions of space2.C().coldim() and
        space1.C().coldim().
        Throws if the overlap cannot be computed.
      */
    static boost::shared_ptr<Matrix> overlap(const MOIndexSpace& space2, const MOIndexSpace& space1);
};

}

#endif // _psi_src_lib_libmints_moindexspace_h_
