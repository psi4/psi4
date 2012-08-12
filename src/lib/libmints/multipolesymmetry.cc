#include <psi4-dec.h>
#include "mints.h"

using namespace boost;
using namespace std;

namespace psi{

OperatorSymmetry::OperatorSymmetry(int order,
                                   boost::shared_ptr<Molecule> mol,
                                   boost::shared_ptr<IntegralFactory> ints)
    : order_(order), molecule_(mol), integral_(ints)
{
    common_init();
}

OperatorSymmetry::OperatorSymmetry(int order,
                                   boost::shared_ptr<Molecule> mol,
                                   boost::shared_ptr<IntegralFactory> ints,
                                   boost::shared_ptr<MatrixFactory> mats)
    : order_(order), molecule_(mol), integral_(ints), matrix_(mats)
{
    common_init();
}

void OperatorSymmetry::common_init()
{
    // Nabla operator has same symmetry as dipole
    if (order_ == P)
        order_ = 1;

    if (order_ > 0) {
        int ncart = INT_NCART(order_);
        component_symmetry_.resize(ncart, 0);

        CharacterTable ct = molecule_->point_group()->char_table();
        SymmetryOperation so;
        int nirrep = ct.nirrep();

        double *t = new double[ncart];

        for (int irrep=0; irrep<nirrep; ++irrep) {
            IrreducibleRepresentation gamma = ct.gamma(irrep);

            ::memset(t, 0, sizeof(double)*ncart);

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
        OperatorSymmetry quad(2, molecule_, integral_, matrix_);

        // Make sure order is 1 for the other routines in this class.
        order_ = 1;

        int n = 3;
        component_symmetry_.resize(n, 0);

        component_symmetry_[0] = quad.component_symmetry(4);  // Lx === yz
        component_symmetry_[1] = quad.component_symmetry(2);  // Ly === xz
        component_symmetry_[2] = quad.component_symmetry(1);  // Lz === xy
    }
    else {
        throw PSIEXCEPTION("MultipoleSymmetry: Don't understand the multipole order given.");
    }
}

OperatorSymmetry::~OperatorSymmetry()
{
}

string OperatorSymmetry::form_suffix(int x, int y, int z)
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

string OperatorSymmetry::name_of_component(int i)
{
    Vector3 components = BasisSet::exp_ao[order_][i];
    return form_suffix(components[0], components[1], components[2]);
}

vector<SharedMatrix > OperatorSymmetry::create_matrices(const std::string &basename)
{
    vector<SharedMatrix > matrices;
    if (matrix_) {
        string name;
    
        for (int i=0; i<INT_NCART(order_); ++i) {
            name = basename + " " + name_of_component(i);
    
            matrices.push_back(matrix_->create_shared_matrix(name, component_symmetry_[i]));
        }
    }

    return matrices;
}



MultipoleSymmetry::MultipoleSymmetry(int order,
                                   boost::shared_ptr<Molecule> mol,
                                   boost::shared_ptr<IntegralFactory> ints)
    : order_(order), molecule_(mol), integral_(ints)
{
    common_init();
}

MultipoleSymmetry::MultipoleSymmetry(int order,
                                   boost::shared_ptr<Molecule> mol,
                                   boost::shared_ptr<IntegralFactory> ints,
                                   boost::shared_ptr<MatrixFactory> mats)
    : order_(order), molecule_(mol), integral_(ints), matrix_(mats)
{
    common_init();
}

void MultipoleSymmetry::common_init()
{
    int n_components = (order_+1)*(order_+2)*(order_+3)/6 - 1;
    component_symmetry_.resize(n_components, 0);

    int l_offset = 0;
    for(int l = 1; l <= order_; ++l){
        int ncart = INT_NCART(l);

        CharacterTable ct = molecule_->point_group()->char_table();
        int nirrep = ct.nirrep();

        double *t = new double[ncart];

        for (int irrep=0; irrep<nirrep; ++irrep) {
            IrreducibleRepresentation gamma = ct.gamma(irrep);

            ::memset(t, 0, sizeof(double)*ncart);

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

            for (int xyz=0; xyz<ncart; ++xyz) {
                if (t[xyz] != 0) {
                    component_symmetry_[l_offset + xyz]= irrep;
//                    fprintf(outfile, "Component %d %d = %d\n", l, l_offset+xyz, irrep);
                }
            }
        }
        l_offset+= ncart;

        delete[] t;
    }
}

MultipoleSymmetry::~MultipoleSymmetry()
{
}



vector<SharedMatrix > MultipoleSymmetry::create_matrices(const std::string &basename)
{
    vector<SharedMatrix > matrices;

    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();

    int component = 0;
    for(int l = 1; l <= order_; ++l){
        for(int ii = 0; ii <= l; ii++) {
            int lx = l - ii;
            for(int lz = 0; lz <= ii; lz++) {
                int ly = ii - lz;
                std::stringstream sstream;
                if(l == 1){
                    sstream << "Dipole ";
                }else if(l == 2){
                    sstream << "Quadrupole ";
                }else if (l == 3){
                    sstream << "Octupole ";
                }else if (l == 4){
                    sstream << "Hexadecapole ";
                }else{
                    int n = (1 << l);
                    sstream << n << "-pole ";
                }
                std::string name = sstream.str();
                for(int xval = 0; xval < lx; ++xval)
                    name += "X";
                for(int yval = 0; yval < ly; ++yval)
                    name += "Y";
                for(int zval = 0; zval < lz; ++zval)
                    name += "Z";

                name = basename + name;
                int sym = component_symmetry_[component];
                if (matrix_) {
                    matrices.push_back(matrix_->create_shared_matrix(name, component_symmetry_[component]));
                }else{
                    SharedMatrix mat(new Matrix(name, wfn->nsopi(), wfn->nsopi(), sym));
                    matrices.push_back(mat);
                }
                ++component;
            }
        }
    }

    return matrices;
}

}
