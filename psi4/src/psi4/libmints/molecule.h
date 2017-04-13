/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_molecule_h_
#define _psi_src_lib_libmints_molecule_h_

#include <vector>
#include <string>
#include <cstdio>
#include <map>
#include <memory>

#include "psi4/pybind11.h"

#define LINEAR_A_TOL 1.0E-2 //When sin(a) is below this, we consider the angle to be linear
#define DEFAULT_SYM_TOL 1.0E-8
#define FULL_PG_TOL 1.0e-8 // default

#include "typedefs.h"
#include "coordentry.h"

namespace psi {
class PointGroup;

enum RotorType {RT_ASYMMETRIC_TOP, RT_SYMMETRIC_TOP, RT_SPHERICAL_TOP, RT_LINEAR, RT_ATOM};
enum FullPointGroup {PG_ATOM, PG_Cinfv, PG_Dinfh, PG_C1, PG_Cs, PG_Ci, PG_Cn, PG_Cnv,
 PG_Cnh, PG_Sn, PG_Dn, PG_Dnd, PG_Dnh, PG_Td, PG_Oh, PG_Ih};

const std::string RotorTypeList[] = {"ASYMMETRIC_TOP", "SYMMETRIC_TOP",
 "SPHERICAL_TOP", "LINEAR", "ATOM"};

const std::string FullPointGroupList[] = {"ATOM", "C_inf_v", "D_inf_h", "C1", "Cs", "Ci", "Cn", "Cnv",
 "Cnh", "Sn", "Dn", "Dnd", "Dnh", "Td", "Oh", "Ih"};

/*! \ingroup MINTS
 *  \class Molecule
 *  \brief Molecule information class.
 */
class Molecule
{
public:
    /**
     * The type of geometry provided in the input
     */
    enum GeometryFormat {
        ZMatrix,     /*!< Z-matrix coordinates */
        Cartesian    /*!< Cartesian coordinates */
    };
    /**
     * The Units used to define the geometry
     */
    enum GeometryUnits {Angstrom, Bohr};
    /**
     * How to handle each fragment
     */
    enum FragmentType {
        Absent,  /*!< Neglect completely */
        Real,    /*!< Include, as normal */
        Ghost    /*!< Include, but with ghost atoms */
    };

    enum FragmentLevel {
        QMatom  = 1,    /*!< Quantum mechanical */
        EFPatom = 2,    /*!< Effective fragment potential */
        ALLatom = 3     /*!< All atom types */
    };

    typedef std::vector<std::shared_ptr<CoordEntry> > EntryVector;
    typedef EntryVector::iterator EntryVectorIter;

protected:
    /// Molecule (or fragment) name
    std::string name_;
    /// Atom info vector (no knowledge of dummy atoms)
    EntryVector atoms_;
    /// Atom info vector (includes dummy atoms)
    EntryVector full_atoms_;
    /// The charge of each fragment
    std::vector<int> fragment_charges_;
    /// The multiplicity of each fragment
    std::vector<int> fragment_multiplicities_;

    /// Reorient or not?
    bool fix_orientation_;
    /// Move to center of mass or not?
    bool move_to_com_;

    /// Whether the charge was given by the user
    bool charge_specified_;
    /// Whether the multiplicity was specified by the user
    bool multiplicity_specified_;
    /// The molecular charge
    int molecular_charge_;
    /// The multiplicity (defined as 2Ms + 1)
    int multiplicity_;
    /// The units used to define the geometry
    GeometryUnits units_;
    /// The conversion factor to take input units to Bohr
    double input_units_to_au_;
    /// A list of all variables known, whether they have been set or not.
    std::vector<std::string> all_variables_;
    /// Zero it out
    void clear();

    /**
     * Attempts to interpret a string as a double, if not it assumes it's a variable.
     *
     * @param str the string to interpret.
     * @return the CoordValue interpretation of the string.
     */
    CoordValue* get_coord_value(const std::string &str);

    /**
     * Attempts to interpret a string as an atom specifier in a zmatrix.
     *
     * @param str the string to interpret.
     * @param line the current line, for error message printing.
     * @return the atom number (adjusted to zero-based counting)
     */
    int get_anchor_atom(const std::string &str, const std::string &line);

    /// Point group to use with this molecule.
    std::shared_ptr<PointGroup> pg_;
    /// Full point group.
    FullPointGroup full_pg_;
    /// n of the highest rotational axis Cn
    int full_pg_n_;

    /// Number of unique atoms
    int nunique_;
    /// Number of equivalent atoms per unique atom (length nunique_)
    int *nequiv_;
    /// Equivalent atom mapping array
    int **equiv_;
    /// Atom to unique atom mapping array (length natom)
    int *atom_to_unique_;

    /// A listing of the variables used to define the geometries
    std::map<std::string, double> geometry_variables_;
    /// A list describing how to handle each fragment
    std::vector<FragmentType> fragment_types_;
//****AVC****//
// moved fragments_ to public
// moved fragment_levels_ to public
//****AVC****//
    /// Symmetry string from geometry specification
    std::string symmetry_from_input_;
    /// Reinterpret the coord entries or not
    /// Default is true, except for findif
    bool reinterpret_coordentries_;
    /// Nilpotence boolean (flagged upon first determination of symmetry frame, reset each time a substantiative change is made)
    bool lock_frame_;
    /// Whether this molecule has at least one zmatrix entry
    bool zmat_;

public:
//****AVC****//
    /// The list of atom ranges defining each fragment from parent molecule
    std::vector<std::pair<int, int> > fragments_;
//****AVC****//

    Molecule();
    /// Copy constructor.
    Molecule(const Molecule& other);
    virtual ~Molecule();

    Molecule clone(void) {
      Molecule new_obj(*this);
      return new_obj;
    }

    /// @{
    /// Operators
    /// Assignment operator.
    Molecule& operator=(const Molecule& other);
    /// Addition
    Molecule operator+(const Molecule& other);
    /// Subtraction
    Molecule operator-(const Molecule& other);
    /// Plus equals
    void operator+=(const Molecule& other);
    /// @}


    /**
     * Pull information from an XYZ file. Useful for debugging.
     * \param xyzfilename Filename of xyz file.
     */
    void init_with_xyz(const std::string& xyzfilename);

    /**
     * Add an atom to the molecule
     * \param Z atomic number
     * \param x cartesian coordinate
     * \param y cartesian coordinate
     * \param z cartesian coordinate
     * \param symb atomic symbol to use
     * \param mass mass to use if non standard
     * \param charge charge to use if non standard
     * \param lineno line number when taken from a string
     */
    void add_atom(int Z, double x, double y, double z,
                  const char *symb = "", double mass = 0.0,
                  double charge = 0.0, int lineno = -1);

    /// Whether the multiplicity was given by the user
    bool multiplicity_specified() const { return multiplicity_specified_; }
    /// Whether the charge was given by the user
    bool charge_specified() const { return charge_specified_; }
    /// The number of fragments in the molecule
    int nfragments() const { return fragments_.size();}
    /// The number of active fragments in the molecule
    int nactive_fragments();
    /// Returns the list of atoms belonging to a fragment.
    // Needed for EFP interface
    std::pair<int, int> fragment_atom_pair(int f) { return fragments_[f]; }

    /// Get molecule name
    const std::string name() const {return name_; }
    /// Returns the name of the basis set on the specified atom
    const std::string& basis_on_atom(int atom) const;
    /// Set molecule name
    void set_name(const std::string &_name) { name_ = _name; }
    /// Number of atoms
    int natom() const;
    /// Number of all atoms (includes dummies)
    int nallatom() const { return full_atoms_.size(); }
    /// Nuclear charge of atom
    const double& Z(int atom) const;
    /// Nuclear charge of atom
    double fZ(int atom) const;
    /// x position of atom
    double x(int atom) const;
    /// y position of atom
    double y(int atom) const;
    /// z position of atom
    double z(int atom) const;
    /// x position of atom
    double fx(int atom) const;
    /// y position of atom
    double fy(int atom) const;
    /// z position of atom
    double fz(int atom) const;
    /// Returns a Vector3 with x, y, z position of atom
    Vector3 xyz(int atom) const;
    Vector3 fxyz(int atom) const;
    /// Returns x, y, or z component of 'atom'
    double xyz(int atom, int _xyz) const;
    /// Returns mass atom atom
    double mass(int atom) const;

    /// Set the mass of a particular atom (good for isotopic substitutions)
    void set_mass(int atom, double mass) { full_atoms_[atom]->set_mass(mass); }
    /// Set the nuclear charge of an atom (primarily used in ECP calculations).
    void set_nuclear_charge(int atom, double newZ);

    /// Returns the cleaned up label of the atom (C2 => C, H4 = H)
    std::string symbol(int atom) const;
    /// Returns the cleaned up label of the atom (C2 => C, H4 = H)
    std::string fsymbol(int atom) const;
    /// Returns the original label of the atom as given in the input file (C2, H4).
    std::string label(int atom) const;
    /// Returns charge of atom
    double charge(int atom) const;
    /// Returns the true atomic number of an atom
    int true_atomic_number(int atom) const;
    int ftrue_atomic_number(int atom) const;
    /// Returns mass atom atom
    double fmass(int atom) const;
    /// Returns label of atom
    std::string flabel(int atom) const;
    /// Returns charge of atom
    double fcharge(int atom) const;
    /// Returns the CoordEntry for an atom
    const std::shared_ptr<CoordEntry>& atom_entry(int atom) const;

    void set_basis_all_atoms(const std::string& name, const std::string& type="BASIS");
    void set_basis_by_symbol(const std::string& symbol, const std::string& name, const std::string& type="BASIS");
    void set_basis_by_number(int number, const std::string& name, const std::string& type="BASIS");
    void set_basis_by_label(const std::string& label, const std::string& name, const std::string& type="BASIS");
    void set_shell_by_label(const std::string& label, const std::string& name, const std::string& type="BASIS");

    /// Number of frozen core for molecule given freezing state
    int nfrozen_core(const std::string& depth = "");

    /// @{
    /// Tests to see of an atom is at the passed position with a given tolerance
    int atom_at_position1(double *, double tol = 0.05) const;
    int atom_at_position2(Vector3&, double tol = 0.05) const;
    /// @}

    /// Do we reinterpret coordentries during a call to update_geometry?
    void set_reinterpret_coordentry(bool rc);

    /// Returns the geometry in a Matrix
    Matrix geometry() const;
    /// Returns the full (dummies included) in a Matrix
    Matrix full_geometry() const;

    /**
     * Sets the geometry, given a matrix of coordinates (in Bohr).
     * \param geom geometry matrix of dimension natom X 3
     */
    void set_geometry(double** geom);

    /**
     * Sets the geometry, given a Matrix of coordinates (in Bohr).
     */
    void set_geometry(const Matrix& geom);

    /**
     * Sets the full geometry, given a matrix of coordinates (in Bohr).
     */
    void set_full_geometry(double** geom);

    /**
     * Sets the full geometry, given a Matrix of coordinates (in Bohr).
     */
    void set_full_geometry(const Matrix& geom);

    /**
     * Rotates the molecule using rotation matrix R
     */
    void rotate(const Matrix& R);
    void rotate_full(const Matrix& R);

    /**
     * Reinterpret the fragments for reals/ghosts and build the atom list
     */
    void reinterpret_coordentries();

    /**
     * Reinterpret the fragments for QM/EFP and build the atom list
     */
    void reinterpret_fragments();

    /**
     * Find the nearest point group within the tolerance specified, and adjust
     * the coordinates to have that symmetry.
     */
    void symmetrize_to_abelian_group(double tol);

    /// Computes center of mass of molecule (does not translate molecule)
    Vector3 center_of_mass() const;
    /// Computes nuclear repulsion energy
    double nuclear_repulsion_energy() const;
    /// The dipole generated by the nuclear charges evaluated at the origin.
    Vector3 nuclear_dipole() const;
    /// The dipole generated by the nuclear charges evaluated using a given origin.
    Vector3 nuclear_dipole(const Vector3 &origin) const;
    /// Computes nuclear repulsion energy derivatives.
    Matrix nuclear_repulsion_energy_deriv1() const;
    /// Computes nuclear repulsion energy second derivatives.
    Matrix nuclear_repulsion_energy_deriv2() const;

    /// Computes the nuclear repuslion energy between this and another Molecule
    double pairwise_nuclear_repulsion_energy(std::shared_ptr<Molecule> other) const;

    /// Translates molecule by r
    void translate(const Vector3& r);
    /// Moves molecule to center of mass
    void move_to_com();
    /**
     *  Reorient molecule to standard frame. See input/reorient.cc
     *  If you want the molecule to be reoriented about the center of mass
     *  make sure you call move_to_com() prior to calling reorient()
     */
//    void reorient();

    /// Computes and returns a matrix depicting distances between atoms.
    Matrix distance_matrix() const;

    /// Compute inertia tensor.
    Matrix* inertia_tensor() const;

    /// Compute the rotational constants and return them in wavenumbers
    Vector rotational_constants(double tol = FULL_PG_TOL) const;

    /// Print the rotational constants
    void print_rotational_constants(void) const;
    /// Return the rotor type
    RotorType rotor_type(double tol = FULL_PG_TOL) const;

    /// Print the molecule
    void print() const;

    /// Print the molecule, adding spacers
    void print_cluster() const;

    /// Print full atom list
    void print_full() const;

    /// Print the molecule in the same format that the user provided
    void print_in_input_format() const;

    /// Print the molecule in Bohr
    void print_in_bohr() const;

    /// Print the molecule in Angstrom
    void print_in_angstrom() const;

    ///Print the geometrical parameters of the molecule
    void print_distances() const;
    void print_bond_angles() const;
    void print_dihedrals() const;
    void print_out_of_planes() const;

    /// Save an XYZ file
    void save_xyz_file(const std::string & filename, bool save_ghosts = true) const;
    /// Save an XYZ file to a string
    std::string save_string_xyz_file() const;

    /// Save an XYZ string
    std::string save_string_xyz() const;

    ///
    /// Unique details
    /// These routines assume the molecular point group has been determined.
    /// @{
    /// Return the number of unique atoms.
    int nunique() const { return nunique_; }
    /// Returns the overall number of the iuniq'th unique atom.
    int unique(int iuniq) const { return equiv_[iuniq][0]; }
    /// Returns the number of atoms equivalent to iuniq.
    int nequivalent(int iuniq) const { return nequiv_[iuniq]; }
    /// Returns the j'th atom equivalent to iuniq.
    int equivalent(int iuniq, int j) const { return equiv_[iuniq][j]; }
    /** Converts an atom number to the number of its generating unique atom.
        The return value is in [0, nunique). */
    int atom_to_unique(int iatom) const { return atom_to_unique_[iatom]; }
    /** Converts an atom number to the offset of this atom in the list of
        generated atoms. The unique atom itself is allowed offset 0. */
    int atom_to_unique_offset(int iatom) const;
    /** Returns the maximum number of equivalent atoms. */
    int max_nequivalent() const;
    /// @}

    ///
    /// Symmetry
    /// @{
    bool has_symmetry_element(Vector3& op, double tol=DEFAULT_SYM_TOL) const;
    std::shared_ptr<PointGroup> point_group() const;
    void set_point_group(std::shared_ptr<PointGroup> pg);
    /// Determine and set FULL point group
    void set_full_point_group(double tol=FULL_PG_TOL);
    /// Does the molecule have an inversion center at origin
    bool has_inversion(Vector3& origin, double tol=DEFAULT_SYM_TOL) const;
    /// Is a plane?
    bool is_plane(Vector3& origin, Vector3& uperp, double tol=DEFAULT_SYM_TOL) const;
    /// Is an axis?
    bool is_axis(Vector3& origin, Vector3& axis, int order, double tol=DEFAULT_SYM_TOL) const;
    /// Is the molecule linear, or planar?
    void is_linear_planar(bool& linear, bool& planar, double tol=DEFAULT_SYM_TOL) const;
    /// Find computational molecular point group, user can override this with the "symmetry" keyword
    std::shared_ptr<PointGroup> find_point_group(double tol=DEFAULT_SYM_TOL) const;
    /// Override symmetry from outside the molecule string
    void reset_point_group(const std::string& pgname);
    /// Find highest molecular point group
    std::shared_ptr<PointGroup> find_highest_point_group(double tol=DEFAULT_SYM_TOL) const;
    /// Determine symmetry reference frame. If noreorient is set, this is the rotation matrix
    /// applied to the geometry in update_geometry.
    std::shared_ptr<Matrix> symmetry_frame(double tol=DEFAULT_SYM_TOL);
    /// Release symmetry information
    void release_symmetry_information();
    /// Initialize molecular specific symemtry information
    /// Uses the point group object obtain by calling point_group()
    void form_symmetry_information(double tol=DEFAULT_SYM_TOL);
    /// Returns the symmetry label
    std::string sym_label();
    /// Returns the irrep labels
    char **irrep_labels();
    const std::string& symmetry_from_input() const { return symmetry_from_input_; }

    /**
     * Force the molecule to have the symmetry specified in pg_.
     * This is to handle noise coming in from optking.
     */
    void symmetrize(double tol=0.05, bool suppress_mol_print_in_exc=false);
    /// @}

    /**
     * Given a string (including newlines to separate lines), builds a new molecule
     * and wraps it in a smart pointer
     *
     * @param geom a string providing the user's input
     */
    static std::shared_ptr<Molecule> create_molecule_from_string(const std::string &geom);

    /**
     * Regenerates a input file molecule specification string
     * from the current state of the Molecule. Contains Cartesian
     * geometry info, fragmentation, charges and multiplicities,
     * and any frame restriction.
     */
    std::string create_psi4_string_from_molecule() const;

    /**
     * Sets all fragments in the molecule to be active.
     */
    void activate_all_fragments();

    /**
     * Sets all fragments in the molecule to be inactive.
     */
    void deactivate_all_fragments();

    /**
     * Sets the specified list of fragments to be real.
     * @param reals The list of real fragments.
     */
    void set_active_fragments(py::list reals);

    /**
     * Sets the specified fragment to be real.
     * @param fragment The fragment to set.
     */
    void set_active_fragment(int fragment);

    /**
     * Sets the specified list of fragments to be ghosts.
     * @param ghosts The list of ghosts fragments.
     */
    void set_ghost_fragments(py::list ghosts);

    /**
     * Sets the specified fragment to be a ghost.
     * @param fragment The fragment to set.
     */
    void set_ghost_fragment(int fragment);

    /**
     * Makes a copy of the molecule, returning a new ref counted molecule with
     * only certain fragment atoms present as either ghost or real atoms
     * @param real_list The list of fragments that should be present in the molecule as real atoms.
     * @param ghost_list The list of fragments that should be present in the molecule as ghosts.
     * @return The ref counted cloned molecule
     */
    std::shared_ptr<Molecule> extract_subsets(const std::vector<int> &real_list,
                                                const std::vector<int> &ghost_list) const;

    /**
     * A wrapper to extract_subsets, callable from Boost
     * @param reals A list containing the real atoms.
     * @param ghost A list containing the ghost atoms.
     * @return The ref counted cloned molecule.
     */
    std::shared_ptr<Molecule> py_extract_subsets_1(py::list reals,
                                                   py::list ghost);

    /**
     * A wrapper to extract_subsets, callable from Boost
     * @param reals A list containing the real atoms.
     * @param ghost An int containing the ghost atoms.
     * @return The ref counted cloned molecule.
     */
    std::shared_ptr<Molecule> py_extract_subsets_2(py::list reals,
                                                   int ghost = -1);

    /**
     * A wrapper to extract_subsets, callable from Boost
     * @param reals An int containing the real atoms.
     * @param ghost A list containing the ghost atoms.
     * @return The ref counted cloned molecule.
     */
    std::shared_ptr<Molecule> py_extract_subsets_3(int reals,
                                                   py::list ghost);

    /**
     * A wrapper to extract_subsets, callable from Boost
     * @param reals An int containing the real atoms.
     * @param ghost An int containing the ghost atoms.
     * @return The ref counted cloned molecule.
     */
    std::shared_ptr<Molecule> py_extract_subsets_4(int reals,
                                                   int ghost = -1);

    /**
     * A wrapper to extract_subsets, callable from Boost
     * @param reals A list containing the real atoms.
     * @return The ref counted cloned molecule.
     */
    std::shared_ptr<Molecule> py_extract_subsets_5(py::list reals);

    /**
     * A wrapper to extract_subsets, callable from Boost
     * @param reals An int containing the real atoms.
     * @return The ref counted cloned molecule.
     */
    std::shared_ptr<Molecule> py_extract_subsets_6(int reals);

    // => Fragment Composition <= //

    /// The list of atom ranges defining each fragment from parent molecule (fragments[frag_ind] = <Afirst,Alast+1>)
    const std::vector<std::pair<int, int> >& fragments() const { return fragments_; }
    /// A list describing how to handle each fragment
    const std::vector<FragmentType>& fragment_types() const { return fragment_types_; }
    /// The charge of each fragment
    const std::vector<int>& fragment_charges() const { return fragment_charges_; }
    /// The multiplicity of each fragment
    const std::vector<int>& fragment_multiplicities() const { return fragment_multiplicities_; }

    /// Sets whether this molecule contains at least one zmatrix entry
    void set_has_zmatrix(bool tf) {zmat_ = tf;}
    /// Whether this molecule has at least one zmatrix entry
    bool has_zmatrix() const {return zmat_;}
    /// Assigns the value val to the variable labelled string in the list of geometry variables.
    /// Also calls update_geometry()
    void set_variable(const std::string &str, double val);
    /// Checks to see if the variable str is in the list, sets it to val and returns
    /// true if it is, and returns false if not.
    double get_variable(const std::string &str);
    /// Checks to see if the variable str is in the list, returns
    /// true if it is, and returns false if not.
    bool is_variable(const std::string &str) const;

    /// Sets the molecular charge
    void set_molecular_charge(int charge) {charge_specified_ = true; molecular_charge_ = charge;}
    /// Gets the molecular charge
    int molecular_charge() const;
    /// Sets the multiplicity (defined as 2Ms + 1)
    void set_multiplicity(int mult) { multiplicity_specified_ = true; multiplicity_ = mult; }
    /// Get the multiplicity (defined as 2Ms + 1)
    int multiplicity() const;
    /// Sets the geometry units
    void set_units(GeometryUnits units) { units_ = units; }
    /// Gets the geometry units
    GeometryUnits units() const { return units_; }

    /// Get whether or not orientation is fixed
    bool orientation_fixed() const { return fix_orientation_; }
    /// Fix the orientation at its current frame
    void set_orientation_fixed(bool fix = true) { fix_orientation_ = fix;}
    /// Fix the center of mass at its current frame
    void set_com_fixed(bool fix = true) {move_to_com_ = !fix;}
    /// Returns the Schoenflies symbol
    std::string schoenflies_symbol() const;
    /// Check if current geometry fits current point group
    bool valid_atom_map(double tol = 0.05) const;
    /// Return point group name such as C3v or S8.
    std::string full_point_group() const;
    /// Return point group name such as Cnv or Sn.
    std::string full_point_group_with_n() const { return FullPointGroupList[full_pg_];}
    /// Return n in Cnv, etc.; If there is no n (e.g. Td) it's the highest-order rotation axis.
    int full_pg_n() { return full_pg_n_; }

    /**
     * Updates the geometry, by (re)interpreting the string used to create the molecule, and the current values
     * of the variables. The atoms list is cleared, and then rebuilt by this routine.
     */
    void update_geometry();
};

typedef std::shared_ptr<Molecule> SharedMolecule;

}

#endif
