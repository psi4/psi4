#ifndef _psi_src_lib_libmints_molecule_h_
#define _psi_src_lib_libmints_molecule_h_

#include <vector>
#include <string>
#include <cstdio>
#include <map>

#include <boost/regex.hpp>
#include <boost/shared_ptr.hpp>

#define LINEAR_A_TOL 1.0E-2 //When sin(a) is below this, we consider the angle to be linear
namespace psi {

class PSIO;
class Chkpt;
class PointGroup;
class Matrix;
class SimpleMatrix;
class Vector;
class SimpleVector;
class Vector3;

extern FILE *outfile;

/*! \ingroup MINTS
 *  \class Molecule
 *  \brief Molecule information class.
 */
class Molecule
{
public:
    struct atom_info {
        double x, y, z;
        int Z;				// if Z == dummy atom
        double charge;
        double mass;
        std::string label;
    };

    /**
     * The type of geometry provided in the input
     * None      - used to identify if the initial parsing has been done
     * Cartesian - Cartesian coordinates
     * ZMatrix   - Z matrix coordinates
     */
    enum GeometryFormat {None, ZMatrix, Cartesian};
    /**
     * The Units used to define the geometry
     */
    enum GeometryUnits {Angstrom, Bohr};

    static boost::shared_ptr<Molecule> create_molecule_from_string(const std::string &geom);
    /// Assigns the value val to the variable labelled string in the list of geometry variables.
    /// Also calls update_geometry()
    void set_variable(const std::string &str, double val);
    /// Checks to see if the variable str is in the list, sets it to val and returns
    /// true if it is, and returns false if not.
    double get_variable(const std::string &str);
    /// Checks to see if the variable str is in the list, returns
    /// true if it is, and returns false if not.
    bool is_variable(const std::string &str) const;

    /// Sets the geometry string
    void set_geometry_string(std::vector<std::string> strVec) { geometryString_ = strVec; }
    /// Sets the geometry format
    void set_geometry_format(GeometryFormat format) { geometryFormat_ = format; }
    /// Gets the geometry format
    GeometryFormat geometry_format() const { return geometryFormat_; }
    /// Sets the molcular charge
    void set_charge(int charge) {charge_ = charge;}
    /// Gets the molecular charge
    //int charge() const {return charge_;}
    /// Sets the multiplicity (defined as 2Ms + 1)
    void set_multiplicity(int mult) { multiplicity_ = mult; }
    /// Get the multiplicity (defined as 2Ms + 1)
    int multiplicity() const { return multiplicity_; }
    /// Sets the geometry units
    void set_units(GeometryUnits units) { units_ = units; }
    /// Gets the geometry units
    GeometryUnits units() const { return units_; }

    void update_geometry();
protected:
    /// Molecule (or fragment) name
    std::string name_;
    ///A regex match object to receive captured expressions in regex searches
    boost::smatch reMatches_;
    /// Atom info vector (no knowledge of dummy atoms)
    std::vector<atom_info> atoms_;
    /// Atom info vector (includes dummy atoms)
    std::vector<atom_info> full_atoms_;
    /// Each line of the string passed in to define the geometry
    std::vector<std::string> geometryString_;
    /// Symmetry information about the molecule
    int nirreps_;
    /// The molecular charge
    int charge_;
    /// The multiplicity (defined as 2Ms + 1)
    int multiplicity_;
    /// The units used to define the geometry
    GeometryUnits units_;
    /// Zero it out
    void clear();
    double get_value(const std::string &str, const std::string &line);
    int get_anchor_atom(const std::string &str, const std::vector<std::string> &atoms,
                              const std::string &line);

    // Point group to use with this molecule.
    boost::shared_ptr<PointGroup> pg_;

    /// Number of unique atoms
    int nunique_;
    int *nequiv_;
    int **equiv_;
    int *atom_to_unique_;

    /// A regular expression to test if a string looks like a floating point number
    static boost::regex realNumber_;
    /// A regular expression to test if a string looks like an integer
    static boost::regex integerNumber_;
    /// A regular expression to test if a string looks like an atom symbol
    static boost::regex atomSymbol_;
    /// A regular expression to test if a string looks like a variable assignment (e.g. x = 1.9)
    static boost::regex variableDefinition_;
    /// A regular expression to test if a string is just a blank line
    static boost::regex blankLine_;
    /// A regular expression to test if a string looks like a comment line
    static boost::regex commentLine_;
    /// A regular expression to test if a string looks like a command to specify the units
    static boost::regex unitLabel_;
    /// A regular expression to test if a string looks like a charge/multiplicty definition (e.g. -1 1)
    static boost::regex chargeAndMultiplicity_;

    /// The format of the geometry provided by the user
    GeometryFormat geometryFormat_;

    /// A listing of the variables used to define the geometries
    std::map<std::string, double> geometryVariables_;

public:
    Molecule();
    /// Copy constructor.
    Molecule(const Molecule& other);
    virtual ~Molecule();

    /// Assignment operator.
    Molecule& operator=(const Molecule& other);

    /// Pull information from a chkpt object created from psio
    void init_with_psio(boost::shared_ptr<PSIO> psio);
    /// Pull information from the chkpt object passed
    void init_with_chkpt(boost::shared_ptr<Chkpt> chkpt);
    /// Pull information from an XYZ file
    void init_with_xyz(const std::string& xyzfilename);

    /// Add an atom to the molecule
    void add_atom(int Z, double x, double y, double z,
                  const char * = 0, double mass = 0.0,
                  int have_charge = 0, double charge = 0.0);

    /// Get molecule name
    const std::string get_name() const {return name_; }
    /// Set molecule name
    void set_name(const std::string &_name) { name_ = _name; }

    /// Number of atoms
    int natom() const { return atoms_.size(); }
    /// Number of all atoms (includes dummies)
    int nallatom() const { return full_atoms_.size(); }
    /// Nuclear charge of atom
    int Z(int atom) const { return atoms_[atom].Z; }
    /// Nuclear charge of atom
    int fZ(int atom) const { return full_atoms_[atom].Z; }
    // x position of atom
    double x(int atom) const { return atoms_[atom].x; }
    // y position of atom
    double y(int atom) const { return atoms_[atom].y; }
    // z position of atom
    double z(int atom) const { return atoms_[atom].z; }
     // x position of atom
    double fx(int atom) const { return full_atoms_[atom].x; }
    // y position of atom
    double fy(int atom) const { return full_atoms_[atom].y; }
    // z position of atom
    double fz(int atom) const { return full_atoms_[atom].z; }
   /// Return reference to atom_info struct for atom
    const atom_info &r(int atom) const { return atoms_[atom]; }
    /// Return copy of atom_info for atom
    atom_info r(int atom) { return atoms_[atom]; }
   /// Return reference to atom_info struct for atom in full atoms
    const atom_info &fr(int atom) const { return full_atoms_[atom]; }
    /// Return copy of atom_info for atom in full atoms
    atom_info fr(int atom) { return full_atoms_[atom]; }
    /// Returns a Vector3 with x, y, z position of atom
    Vector3 xyz(int atom) const;
    Vector3 fxyz(int atom) const;
    /// Returns x, y, or z component of 'atom'
    double& xyz(int atom, int _xyz);
    const double& xyz(int atom, int _xyz) const;
    /// Returns mass atom atom
    double mass(int atom) const;
    /// Returns label of atom
    const std::string label(int atom) const;
    /// Returns charge of atom
    double charge(int atom) const { return atoms_[atom].charge; }
    /// Returns mass atom atom
    double fmass(int atom) const { return full_atoms_[atom].mass; }
    /// Returns label of atom
    const std::string flabel(int atom) const { return full_atoms_[atom].label; }
    /// Returns charge of atom
    double fcharge(int atom) const { return full_atoms_[atom].charge; }

    /// Number of frozen core for molecule given freezing state
    int nfrozen_core(std::string depth);

    /// Tests to see of an atom is at the passed position with a given tolerance
    int atom_at_position1(double *, double tol = 0.05) const;
    int atom_at_position2(Vector3&, double tol = 0.05) const;

    SimpleMatrix geometry();
    SimpleMatrix full_geometry();
    void set_geometry(SimpleMatrix& geom);
    void set_full_geometry(SimpleMatrix& geom);
    void rotate(SimpleMatrix& R);
    void rotate_full(SimpleMatrix& R);

    /// Computes center of mass of molecule (does not translate molecule)
    Vector3 center_of_mass() const;
    /// Computes nuclear repulsion energy
    double nuclear_repulsion_energy();
    /// Computes nuclear repulsion energy derivatives.
    SimpleMatrix nuclear_repulsion_energy_deriv1();
    /// Computes nuclear repulsion energy second derivatives.
    SimpleMatrix nuclear_repulsion_energy_deriv2();

    /// Returns the nuclear contribution to the dipole moment
    SimpleVector nuclear_dipole_contribution();
    /// Returns the nuclear contribution to the quadrupole moment
    SimpleVector nuclear_quadrupole_contribution();

    /// Translates molecule by r
    void translate(const Vector3& r);
    /// Moves molecule to center of mass
    void move_to_com();
    /** Reorient molecule to standard frame. See input/reorient.cc
     *  If you want the molecule to be reoriented about the center of mass
     *  make sure you call move_to_com() prior to calling reorient()
     */
    void reorient();

    /// Compute inertia tensor.
    SimpleMatrix* inertia_tensor();

    /// Returns the number of irreps
    int nirrep() const { return nirreps_; }
    /// Sets the number of irreps
    void nirrep(int nirreps) { nirreps_ = nirreps; }

    /// Print the molecule
    void print() const;

    /// Save information to checkpoint file.
    void save_to_chkpt(boost::shared_ptr<Chkpt> chkpt, std::string prefix = "");

    //
    // Unique details
    // These routines assume the molecular point group has been determined.
    //
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

    //
    // Symmetry
    //
    boost::shared_ptr<PointGroup> point_group() const { return pg_; }
    void set_point_group(boost::shared_ptr<PointGroup> pg) { pg_ = pg; }
    /// Does the molecule have an inversion center at origin
    bool has_inversion(Vector3& origin, double tol = 0.05) const;
    /// Is a plane?
    bool is_plane(Vector3& origin, Vector3& uperp, double tol = 0.05) const;
    /// Is an axis?
    bool is_axis(Vector3& origin, Vector3& axis, int order, double tol=0.05) const;
    /// Is the molecule linear, or planar?
    void is_linear_planar(bool& linear, bool& planar, double tol) const;
    /// Find highest molecular point group
    boost::shared_ptr<PointGroup> find_point_group(double tol=1.0e-8) const;
    /// Release symmetry information
    void release_symmetry_information();
    /// Initialize molecular specific symemtry information
    /// Uses the point group object obtain by calling point_group()
    void form_symmetry_information(double tol=1.0e-8);
    /// Returns the symmetry label
    const char *sym_label();
    // Returns the irrep labels
    char **irrep_labels();
};

typedef boost::shared_ptr<Molecule> SharedMolecule;

class AtomicRadii
{
protected:
    int count_;
    double *radii_;
public:
    AtomicRadii();
    virtual ~AtomicRadii();

    /** Get the atomic radii.
     *  @param Z atomic number of interest
     *  @returns the atomic radii, or 0 if not found.
     */
    double radii(int Z) const {
        if (Z < 0 || Z >= count_)
            return 0.0;
        return radii_[Z];
    }

    /** Modify the atomic radii. You can only change an existing radii
     *  you cannot add a new one.
     *  @param Z atomic number to modify
     *  @param val new value to use
     */
    void modify_radii(int Z, double val) {
        if (Z < 0 || Z >= count_)
            return;
        radii_[Z] = val;
    }
};

/// Treutler radial mapping radii (see Treutler 1995, Table 1, p. 348).
class TreutlerAtomicRadii : public AtomicRadii
{
public:
    TreutlerAtomicRadii();
};

/** Adjusted van der Waals radii (Angstrom) from atomic ROHF/TZV
 *  computations.
 *  S. Grimme, J. Comput. Chem. 27, 1787-799 (2006)
 */
class AdjustedVanDerWaalsAtomicRadii : public AtomicRadii
{
public:
    AdjustedVanDerWaalsAtomicRadii();
};

/// Bragg-Slater atomic radii
class BraggSlaterAtomicRadii : public AtomicRadii
{
public:
    BraggSlaterAtomicRadii();
};

/** SG1 atomic radii (for use with SG1 grid).
 *  P. M. W. Gill, B. G. Johnson, J. A. Pople, Chem. Phys. Lett.
 *  209, 506 (1993).
 */
class SG1AtomicRadii : public AtomicRadii
{
public:
    SG1AtomicRadii();
};

}

#endif
