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

class OrbitalSpace
{
    /// Number of irreps.
    int nirrep_;

    /// Name of the orbital space.
    std::string name_;

    /// AO->MO transformation matrix (ao x mo)
    SharedMatrix C_;

    /// MO "eigenvalues"
    boost::shared_ptr<Vector> evals_;

    /// AO basis set
    boost::shared_ptr<BasisSet> basis_;
    /// Integral factory that as
    boost::shared_ptr<IntegralFactory> ints_;

    /// MO Dimensionality
    Dimension dim_; // dim_.n() better equal nirrep_

    OrbitalSpace();
public:
    OrbitalSpace(const std::string& name,
                 const SharedMatrix& full_C,   // Should this be 4 C's instead?
                 const boost::shared_ptr<Vector>& evals,
                 const boost::shared_ptr<BasisSet>& basis,
                 const boost::shared_ptr<IntegralFactory>& ints);

    OrbitalSpace(const std::string& name,
                 const boost::shared_ptr<Wavefunction>& wave);

    int nirrep() const;
    const std::string& name() const;

    /// C - transformation matrix (AO x MO)
    const SharedMatrix& C() const;

    /// "Eigenvalues" of the C matrix
    const boost::shared_ptr<Vector>& evals() const;

    /// The AO basis set used to create C
    const boost::shared_ptr<BasisSet>& basis() const;

    /// Integral factory used to create C
    const boost::shared_ptr<IntegralFactory>& integral() const;

    /// MO dimensionality
    const Dimension& dim() const;

    /** Creates an OrbitalSpace from 'from' to the given basis set 'to'
      */
    static OrbitalSpace transform(const OrbitalSpace& from, const boost::shared_ptr<BasisSet>& to);

    /** Returns the overlap matrix between space1 and space2.
        The matrix has dimensions of space2.C().coldim() and
        space1.C().coldim().
        Throws if the overlap cannot be computed.
      */
    static SharedMatrix overlap(const OrbitalSpace& space1, const OrbitalSpace& space2);
    /** Returns the overlap matrix between basis1 and basis2.
        Throws if the overlap cannot be computed.
      */
    static SharedMatrix overlap(const boost::shared_ptr<BasisSet>& basis1,
                                const boost::shared_ptr<BasisSet>& basis2);
};

}

#endif // _psi_src_lib_libmints_moindexspace_h_
