#ifndef _psi_src_lib_libmints_molecule_h_
#define _psi_src_lib_libmints_molecule_h_

#include <vector>
#include <string>
#include <cstdio>

#include <libutil/ref.h>
#include <libmints/vector3.h>
#include <libmints/vector.h>
#include <libmints/matrix.h>

#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>


namespace psi {

extern FILE *outfile;

/*! \ingroup MINTS
 *  \class Molecule
 *  \brief Molecule information class.
 */
class Molecule
{
public:
    typedef struct atom_info {
        double x, y, z;
        int Z;
        double charge;
        double mass;
        std::string label;
    };

protected:
    /// Number of atoms.
    int natoms_;
    /// Atom info vector
    std::vector<atom_info> atoms_;
    /// Symmetry information about the molecule
    int nirreps_;
    /// Zero it out
    void clear();

public:
    Molecule();
    virtual ~Molecule();

    /// Pull information from a chkpt object created from psio
    void init_with_chkpt(shared_ptr<PSIO> psio);
    /// Pull information from the chkpt object passed
    void init_with_chkpt(shared_ptr<Chkpt> chkpt);

    /// Add an atom to the molecule
    void add_atom(int Z, double x, double y, double z,
                  const char * = 0, double mass = 0.0,
                  int have_charge = 0, double charge = 0.0);

    /// Number of atoms
    int natom() const { return natoms_; }
    /// Nuclear charge of atom
    int Z(int atom) const { return atoms_[atom].Z; }
    // x position of atom
    double x(int atom) const { return atoms_[atom].x; }
    // y position of atom
    double y(int atom) const { return atoms_[atom].y; }
    // z position of atom
    double z(int atom) const { return atoms_[atom].z; }
    /// Return reference to atom_info struct for atom
    const atom_info &r(int atom) const { return atoms_[atom]; }
    /// Return copy of atom_info for atom
    atom_info r(int atom) { return atoms_[atom]; }
    /// Returns a Vector3 with x, y, z position of atom
    const Vector3 xyz(int atom) const { return Vector3(atoms_[atom].x, atoms_[atom].y, atoms_[atom].z); }
    /// Returns mass atom atom
    double mass(int atom) const;
    /// Returns label of atom
    const std::string label(int atom) const;
    /// Returns charge of atom
    double charge(int atom) const { return atoms_[atom].charge; }

    /// Tests to see of an atom is at the passed position with a given tolerance
    int atom_at_position(double *, double tol = 0.05) const;

    /// Computes center of mass of molecule (does not translate molecule)
    Vector3 center_of_mass() const;
    /// Computes nuclear repulsion energy
    double nuclear_repulsion_energy();
    /// Computes nuclear repulsion energy derivatives. Free with delete[]
    SimpleVector nuclear_repulsion_energy_deriv1();
    /// Computes nuclear repulsion energy second derivatives.
    SimpleMatrix* nuclear_repulsion_energy_deriv2();

    /// Returns the nuclear contribution to the dipole moment
    SimpleVector nuclear_dipole_contribution();
    /// Returns the nuclear contribution to the quadrupole moment
    SimpleVector nuclear_quadrupole_contribution();

    /// Translates molecule by r
    void translate(const Vector3& r);
    /// Moves molecule to center of mass
    void move_to_com();

    /// Compute inertia tensor.
    SimpleMatrix* inertia_tensor();
    
    /// Returns the number of irreps
    int nirrep() const { return nirreps_; }
    /// Sets the number of irreps
    void nirrep(int nirreps) { nirreps_ = nirreps; }

    /// Print the molecule
    void print(FILE *out = outfile);
};

}

#endif
