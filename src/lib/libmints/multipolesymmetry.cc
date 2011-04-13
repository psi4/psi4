#include <psi4-dec.h>
#include "mints.h"

using namespace boost;
using namespace std;

namespace psi{

MultipoleSymmetry::MultipoleSymmetry(int order,
                                     const boost::shared_ptr<Molecule>& mol,
                                     const boost::shared_ptr<IntegralFactory>& ints,
                                     const boost::shared_ptr<MatrixFactory>& mats)
    : order_(order), molecule_(mol), integral_(ints), matrix_(mats)
{
    // Nabla operator has same symmetry as dipole
    if (order_ == P)
        order_ = 1;

    if (order_ > 0) {
        int ncart = INT_NCART(order_);
        component_symmetry_ = new int[ncart];
        memset(component_symmetry_, 0, sizeof(int)*ncart);

        CharacterTable ct = molecule_->point_group()->char_table();
        SymmetryOperation so;
        int nirrep = ct.nirrep();

        double *t = new double[ncart];

        for (int irrep=0; irrep<nirrep; ++irrep) {
            IrreducibleRepresentation gamma = ct.gamma(irrep);

            memset(t, 0, sizeof(double)*ncart);

            // Apply the projection
            for (int G=0; G<nirrep; ++G) {
                SymmetryOperation so = ct.symm_operation(G);
                ShellRotation rr(order_, so, integral_.get(), false);

                // rr(xyz, xyz) tells us how the orbitals transform in this
                // symmetry operation, then we multiply by the character in
                // the irrep
                for (int xyz=0; xyz<ncart; ++xyz) {
                    t[xyz] += rr(xyz, xyz) * gamma.character(G) / nirrep;
                }
            }

            // Print out t
            for (int xyz=0; xyz<ncart; ++xyz) {
                if (t[xyz] != 0) {
                    component_symmetry_[xyz]= irrep;
                }
            }
        }

        delete[] t;
    }
    else if (order_ == L) {
        // Angular momentum operator is a rotation operator
        // so symmetry of Lz is Lx ^ Ly which can
        // come from quadrupole
        MultipoleSymmetry quad(2, mol, ints, mats);

        // Make sure order is 1 for the other routines in this class.
        order_ = 1;

        int n = 3;
        component_symmetry_ = new int[n];

        component_symmetry_[0] = quad.component_symmetry(4);  // Lx === yz
        component_symmetry_[1] = quad.component_symmetry(2);  // Ly === xz
        component_symmetry_[2] = quad.component_symmetry(1);  // Lz === xy
    }
    else {
        throw PSIEXCEPTION("MultipoleSymmetry: Don't understand the multipole order given.");
    }
}

MultipoleSymmetry::~MultipoleSymmetry()
{
    delete[] component_symmetry_;
}

string MultipoleSymmetry::form_suffix(int x, int y, int z)
{
    string suffix;

    if (x) {
        suffix += "x";
        if (x>1)
            suffix += to_string(x);
    }

    if (y) {
        suffix += "y";
        if (y > 1)
            suffix += to_string(y);
    }

    if (z) {
        suffix += "z";
        if (z > 1)
            suffix += to_string(z);
    }

    return suffix;
}

string MultipoleSymmetry::name_of_component(int i)
{
    Vector3 components = BasisSet::exp_ao[order_][i];
    return form_suffix(components[0], components[1], components[2]);
}

vector<boost::shared_ptr<Matrix> > MultipoleSymmetry::create_matrices(const std::string &basename)
{
    vector<boost::shared_ptr<Matrix> > matrices;
    string name;

    for (int i=0; i<INT_NCART(order_); ++i) {
        name = basename + " " + name_of_component(i);

        matrices.push_back(matrix_->create_shared_matrix(name, component_symmetry_[i]));
    }

    return matrices;
}

}
