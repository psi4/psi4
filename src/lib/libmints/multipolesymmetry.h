#ifndef MULTIPOLESYMMETRY_H
#define MULTIPOLESYMMETRY_H

#include <vector>
#include <string>

namespace boost {
template<class T>
class shared_ptr;
}

namespace psi{

class Molecule;
class IntegralFactory;
class MatrixFactory;

class OperatorSymmetry
{
    // The order of the multipole (dipole=1, quadrupole=2, etc...)
    int order_;

    // Variables we need from the user
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<IntegralFactory> integral_;
    boost::shared_ptr<MatrixFactory> matrix_;

    /**
     * The symmetry of each component of the multipole.
     * Length = INT_NCART(order_)
     */
    std::vector<int> component_symmetry_;

    /**
     * Forms a string that describes what the name of the
     * componet is.
     * (x=1, y=0, z=0) => "x"
     * (x=1, y=0, z=1) => "xz"
     * (x=0, y=2, z=0) => "y2"
     */
    std::string form_suffix(int x, int y, int z);

    void common_init();

public:
    enum Operator {
        Dipole = 1,
        Quadrupole = 2,

        L = -1,
        AngularMomentum = -1,
        P = Dipole,
        Nabla = Dipole
    };

    /** Constructor
     * Constructs an object that determines the symmetry of the different
     * components of the order-th multipole. Orders with negative values
     * have a special meaning (see OperatorSymmetry enum)
     *
     * @param order Order of the multipole (1 = dipole, 2 = quadrupole, etc.)
     * @param mol Molecule the the multipole will be computed for. Needed to obtain
     *            point group object.
     * @param ints Integral factory. Needed for creation of ShellRotation objects.
     * @param mats Matrix factory. Used by create_matrices to create matrices of the
     *             proper size and symmetry.
     */
    OperatorSymmetry(int order,
                     boost::shared_ptr<Molecule> mol,
                     boost::shared_ptr<IntegralFactory> ints,
                     boost::shared_ptr<MatrixFactory> mats);
    OperatorSymmetry(int order,
                     boost::shared_ptr<Molecule> mol,
                     boost::shared_ptr<IntegralFactory> ints);
    virtual ~OperatorSymmetry();

    std::string name_of_component(int i);
    int component_symmetry(int i) const { return component_symmetry_[i]; }

    std::vector<SharedMatrix > create_matrices(const std::string& basename);
};

}

#endif // MULTIPOLESYMMETRY_H
