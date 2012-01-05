#ifndef _psi_src_lib_libmints_molecule_h_
#define _psi_src_lib_libmints_molecule_h_

#include <vector>
#include <string>
#include <cstdio>
#include <map>

#include "typedefs.h"

#define LINEAR_A_TOL 1.0E-2 //When sin(a) is below this, we consider the angle to be linear

#include <boost/shared_ptr.hpp>  // something is going on requiring this header

// Forward declarations for boost.python used in the extract_subsets
namespace boost{
namespace python{
       class list;
}}

namespace psi {
extern FILE *outfile;

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

    typedef std::vector<boost::shared_ptr<CoordEntry> > EntryVector;
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
    /// Whether the user specified the charge, or default was used
    bool charge_specified_;
    /// Whether the user spefified the multiplicity, or default was used
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
    boost::shared_ptr<PointGroup> pg_;

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
    /// The list of atom ranges defining each fragment from parent molecule
    std::vector<std::pair<int, int> > fragments_;
    /// A list describing how to handle each fragment
    std::vector<FragmentType> fragment_types_;
    /// Symmetry string from geometry specification
    std::string symmetry_from_input_;
    /// Old previous symmetry frame (so one can fix to it, if desired)
    Matrix *old_symmetry_frame_;
    /// Old com displacement vector (so one can fix to it, if desired)
    Vector3 *old_com_vector_;

public:
    Molecule();
    /// Copy constructor.
    Molecule(const Molecule& other);
    virtual ~Molecule();

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
     * Pull information from a chkpt object created from psio
     * \param psio PSIO object to initialize with (will create Chkpt object).
     */
    void init_with_psio(boost::shared_ptr<PSIO> psio);

    /**
     * Pull information from the chkpt object passed
     * \param chkpt Chkpt object to initialize with
     */
    void init_with_chkpt(boost::shared_ptr<Chkpt> chkpt);

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

    /// The number of fragments in the molecule
    int nfragments() const { return fragments_.size();}
    /// Get molecule name
    const std::string name() const {return name_; }
    /// Set molecule name
    void set_name(const std::string &_name) { name_ = _name; }
    /// Number of atoms
    const int natom() const { return atoms_.size(); }
    /// Number of all atoms (includes dummies)
    const int nallatom() const { return full_atoms_.size(); }
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
    /// Returns the cleaned up label of the atom (C2 => C, H4 = H)
    std::string symbol(int atom) const;
    /// Returns the original label of the atom as given in the input file (C2, H4).
    std::string label(int atom) const;
    /// Returns charge of atom
    double charge(int atom) const;
    /// Returns mass atom atom
    double fmass(int atom) const;
    /// Returns label of atom
    std::string flabel(int atom) const;
    /// Returns charge of atom
    double fcharge(int atom) const;
    /// Returns the CoordEntry for an atom
    const boost::shared_ptr<CoordEntry>& atom_entry(int atom) const;

    void set_basis_all_atoms(const std::string& name, const std::string& type="BASIS");
    void set_basis_by_symbol(const std::string& symbol, const std::string& name, const std::string& type="BASIS");
    void set_basis_by_number(int number, const std::string& name, const std::string& type="BASIS");
    void set_basis_by_label(const std::string& label, const std::string& name, const std::string& type="BASIS");

    /// Number of frozen core for molecule given freezing state
    int nfrozen_core(const std::string& depth = "");

    /// @{
    /// Tests to see of an atom is at the passed position with a given tolerance
    int atom_at_position1(double *, double tol = 0.05) const;
    int atom_at_position2(Vector3&, double tol = 0.05) const;
    /// @}

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

    /// Computes center of mass of molecule (does not translate molecule)
    Vector3 center_of_mass() const;
    /// Computes nuclear repulsion energy
    double nuclear_repulsion_energy() const;
    /// Computes nuclear repulsion energy derivatives.
    Matrix nuclear_repulsion_energy_deriv1() const;
    /// Computes nuclear repulsion energy second derivatives.
    Matrix nuclear_repulsion_energy_deriv2() const;

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

    /// Returns true if the user specified the charge
    bool charge_specified() const { return charge_specified_; }
    /// Returns true if the user specified the multiplicity
    bool multiplicity_specified() const { return multiplicity_specified_; }

    /// Print the molecule
    void print() const;

    /// Print the molecule in Bohr
    void print_in_bohr() const;

    /// Save an XYZ file
    void save_xyz(const std::string & filename) const;

    /// Save an XYZ string
    std::string save_string_xyz() const;

    /// Save information to checkpoint file.
    void save_to_chkpt(boost::shared_ptr<Chkpt> chkpt, std::string prefix = "");

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
    bool has_symmetry_element(Vector3& op, double tol) const;
    boost::shared_ptr<PointGroup> point_group() const;
    void set_point_group(boost::shared_ptr<PointGroup> pg);
    /// Does the molecule have an inversion center at origin
    bool has_inversion(Vector3& origin, double tol = 0.05) const;
    /// Is a plane?
    bool is_plane(Vector3& origin, Vector3& uperp, double tol = 0.05) const;
    /// Is an axis?
    bool is_axis(Vector3& origin, Vector3& axis, int order, double tol=0.05) const;
    /// Is the molecule linear, or planar?
    void is_linear_planar(bool& linear, bool& planar, double tol) const;
    /// Find computational molecular point group, user can override this with the "symmetry" keyword
    boost::shared_ptr<PointGroup> find_point_group(double tol=1.0e-8) const;
    /// Override symmetry from outside the molecule string
    void reset_point_group(const std::string& pgname);
    /// Find highest molecular point group
    boost::shared_ptr<PointGroup> find_highest_point_group(double tol=1.0e-8) const;
    /// Release symmetry information
    void release_symmetry_information();
    /// Initialize molecular specific symemtry information
    /// Uses the point group object obtain by calling point_group()
    void form_symmetry_information(double tol=1.0e-8);
    /// Returns the symmetry label
    const char *sym_label();
    /// Returns the irrep labels
    char **irrep_labels();
    const std::string& symmetry_from_input() const { return symmetry_from_input_; }

    /**
     * Force the molecule to have the symmetry specified in pg_.
     * This is to handle noise coming in from optking.
     */
    void symmetrize();
    /// @}

    /**
     * Given a string (including newlines to separate lines), builds a new molecule
     * and wraps it in a smart pointer
     *
     * @param text: a string providing the user's input
     */
    static boost::shared_ptr<Molecule> create_molecule_from_string(const std::string &geom);

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
     * @param reals: The list of real fragments.
     */
    void set_active_fragments(boost::python::list reals);

    /**
     * Sets the specified fragment to be real.
     * @param fragment: The fragment to set.
     */
    void set_active_fragment(int fragment);

    /**
     * Sets the specified list of fragments to be ghosts.
     * @param ghosts: The list of ghosts fragments.
     */
    void set_ghost_fragments(boost::python::list ghosts);

    /**
     * Sets the specified fragment to be a ghost.
     * @param fragment: The fragment to set.
     */
    void set_ghost_fragment(int fragment);

    /**
     * Makes a copy of the molecule, returning a new ref counted molecule with
     * only certain fragment atoms present as either ghost or real atoms
     * @param real_list: The list of fragments that should be present in the molecule as real atoms.
     * @param ghost_list: The list of fragments that should be present in the molecule as ghosts.
     * @return The ref counted cloned molecule
     */
    boost::shared_ptr<Molecule> extract_subsets(const std::vector<int> &real_list,
                                                const std::vector<int> &ghost_list) const;

    /**
     * A wrapper to extract_subsets, callable from Boost
     * @param reals: A list containing the real atoms.
     * @param ghost: A list containing the ghost atoms.
     * @return The ref counted cloned molecule.
     */
    boost::shared_ptr<Molecule> py_extract_subsets_1(boost::python::list reals,
                                                   boost::python::list ghost);

    /**
     * A wrapper to extract_subsets, callable from Boost
     * @param reals: A list containing the real atoms.
     * @param ghost: An int containing the ghost atoms.
     * @return The ref counted cloned molecule.
     */
    boost::shared_ptr<Molecule> py_extract_subsets_2(boost::python::list reals,
                                                   int ghost = -1);

    /**
     * A wrapper to extract_subsets, callable from Boost
     * @param reals: An int containing the real atoms.
     * @param ghost: A list containing the ghost atoms.
     * @return The ref counted cloned molecule.
     */
    boost::shared_ptr<Molecule> py_extract_subsets_3(int reals,
                                                   boost::python::list ghost);

    /**
     * A wrapper to extract_subsets, callable from Boost
     * @param reals: An int containing the real atoms.
     * @param ghost: An int containing the ghost atoms.
     * @return The ref counted cloned molecule.
     */
    boost::shared_ptr<Molecule> py_extract_subsets_4(int reals,
                                                   int ghost = -1);

    /**
     * A wrapper to extract_subsets, callable from Boost
     * @param reals: A list containing the real atoms.
     * @return The ref counted cloned molecule.
     */
    boost::shared_ptr<Molecule> py_extract_subsets_5(boost::python::list reals);

    /**
     * A wrapper to extract_subsets, callable from Boost
     * @param reals: An int containing the real atoms.
     * @return The ref counted cloned molecule.
     */
    boost::shared_ptr<Molecule> py_extract_subsets_6(int reals);

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
    void set_molecular_charge(int charge) {molecular_charge_ = charge;}
    /// Gets the molecular charge
    int molecular_charge() const;
    /// Sets the multiplicity (defined as 2Ms + 1)
    void set_multiplicity(int mult) { multiplicity_ = mult; }
    /// Get the multiplicity (defined as 2Ms + 1)
    int multiplicity() const;
    /// Sets the geometry units
    void set_units(GeometryUnits units) { units_ = units; }
    /// Gets the geometry units
    GeometryUnits units() const { return units_; }

    /// Get whether or not orientation is fixed
    bool orientation_fixed() const { return fix_orientation_; }
    /// Fix the orientation at its current frame
    void set_orientation_fixed(bool _fix = true);
    /// Fix the center of mass at its current frame
    void set_com_fixed(bool _fix = true);

    /// Returns the Schoenflies symbol
    std::string schoenflies_symbol() const;

    /**
     * Updates the geometry, by (re)interpreting the string used to create the molecule, and the current values
     * of the variables. The atoms list is cleared, and then rebuilt by this routine.
     */
    void update_geometry();
};

typedef boost::shared_ptr<Molecule> SharedMolecule;

}

#endif
