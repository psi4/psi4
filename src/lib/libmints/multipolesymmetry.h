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

class MultipoleSymmetry
{
    // The order of the multipole (dipole=1, quadrupole=2, etc...)
    int order_;

    // Variables we need from the user
    const boost::shared_ptr<Molecule>& molecule_;
    const boost::shared_ptr<IntegralFactory>& integral_;
    const boost::shared_ptr<MatrixFactory>& matrix_;

    /**
     * The symmetry of each component of the multipole.
     * Length = INT_NCART(order_)
     */
    int *component_symmetry_;

    /**
     * Forms a string that describes what the name of the
     * componet is.
     * (x=1, y=0, z=0) => "x"
     * (x=1, y=0, z=1) => "xz"
     * (x=0, y=2, z=0) => "y2"
     */
    std::string form_suffix(int x, int y, int z);

public:
    MultipoleSymmetry(int order,
                      const boost::shared_ptr<Molecule>& mol,
                      const boost::shared_ptr<IntegralFactory>& ints,
                      const boost::shared_ptr<MatrixFactory>& mats);
    virtual ~MultipoleSymmetry();

    std::string name_of_component(int i);
    int component_symmetry(int i) const { return component_symmetry_[i]; }

    std::vector<boost::shared_ptr<Matrix> > create_matrices(const std::string& basename);
};

}

#endif // MULTIPOLESYMMETRY_H
