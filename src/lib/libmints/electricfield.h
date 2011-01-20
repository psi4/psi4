#ifndef _psi_src_lib_libmints_electricfield_h_
#define _psi_src_lib_libmints_electricfield_h_

namespace psi {

    class BasisSet;
    class GaussianShell;
    class ObaraSaikaTwoCenterRecursion;
    class OneBodyAOInt;
    class IntegralFactory;
    class SphericalTransform;
    class SimpleMatrix;

/*! \ingroup MINTS
 *  \class ElectricFieldInt
 *  \brief Computes electric field integrals.
 *
 *  Use an IntegralFactory to create this object.
 */
class ElectricFieldInt : public OneBodyAOInt
{
    //! Obara and Saika recursion object to be used.
    ObaraSaikaTwoCenterElectricField efield_recur_;                 // Both ElectricField and VIDeriv should give the same result
//    ObaraSaikaTwoCenterVIDerivRecursion efield_recur_;

    //! Number of atoms.
    int natom_;

    //! Computes the electric field between two gaussian shells.
    void compute_pair(const boost::shared_ptr<GaussianShell>&, const boost::shared_ptr<GaussianShell>&);

public:
    //! Constructor. Do not call directly use an IntegralFactory.
    ElectricFieldInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, int deriv=0);
    //! Virtual destructor
    virtual ~ElectricFieldInt();

    //! Compute dipole derivative between two shells, result stored in buffer_.
    //void compute_shell_deriv1(int, int);

    /** Compute all dipole integrals and store them in an array of matrices.
     *  @param result Contains the dipole moment integrals. Order is [mu_x, mu_y, mu_].
     */
//    void compute(std::vector<boost::shared_ptr<SimpleMatrix> > &result);
    /** Compute all dipole derivatives and store them in an array of matrices.
     *  @param result Contains the dipole moment derivative integrals. Order is [mu_x(Aix,Aiy,Aiz...An), mu_y..., mu_z...]
     */
    //void compute_deriv1(std::vector<boost::shared_ptr<SimpleMatrix> > &result);

    //! Does the method provide first derivatives?
    bool has_deriv1() { return false; }
};

}

#endif

