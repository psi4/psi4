#include <cmath>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <limits>

#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>

#include <libmints/molecule.h>
#include <libmints/matrix.h>
#include <libmints/pointgrp.h>
#include <libciomr/libciomr.h>

#include <masses.h>
#include <physconst.h>
#include <element_to_Z.h>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/python.hpp>

using namespace std;
using namespace psi;
using namespace boost;

#include <string>
#include <sstream>
#include <iostream>

// the third parameter of from_string() should be
// one of std::hex, std::dec or std::oct
template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}

// TODO: These should probably be moved to psi4-def.h
#define ZERO_MOMENT_INERTIA 1.0E-10     /*Tolerance for degenerate rotational constants*/
#define ZERO 1.0E-14

namespace psi {

    boost::regex Molecule::realNumber_("(-?\\d+\\.\\d+)|(-?\\d+\\.)|(-?\\.\\d+)|(-?\\d+)", boost::regbase::normal | boost::regbase::icase);
    boost::regex Molecule::integerNumber_("(-?\\d+)", boost::regbase::normal | boost::regbase::icase);
    boost::regex Molecule::atomSymbol_("([A-Z]{1,2})\\d*", boost::regbase::normal | boost::regbase::icase);
    boost::regex Molecule::variableDefinition_("\\s*(\\w+)\\s*=\\s*((-?\\d+\\.\\d+)|(-?\\d+\\.)|(-?\\.\\d+)|(-?\\d+)|(tda))\\s*", boost::regbase::normal | boost::regbase::icase);
    boost::regex Molecule::blankLine_("[\\s%]*", boost::regbase::normal | boost::regbase::icase);
    boost::regex Molecule::commentLine_("\\s*[#%].*", boost::regbase::normal | boost::regbase::icase);
    boost::regex Molecule::unitLabel_("\\s*((ang)|(angstrom)|(bohr)|(au)|(a\\.u\\.))\\s*", boost::regbase::normal | boost::regbase::icase);
    boost::regex Molecule::chargeAndMultiplicity_("\\s*(-?\\d+)\\s+(\\d+)\\s*", boost::regbase::normal);
    boost::regex Molecule::fragmentMarker_("\\s*--\\s*", boost::regbase::normal);

    /**
     * Interprets a string as an integer, throwing if it's unsuccesful.
     */
    int
    str_to_int(const std::string& s)
    {
        int i;
        std::istringstream iss(s);
        if((iss >> std::dec >> i).fail())
            throw PSIEXCEPTION("Unable to convert " + s + " to an integer");
        return i;
    }

    /**
     * Interprets a string as an double, throwing if it's unsuccesful.
     */
    double
    str_to_double(const std::string& s)
    {
        double d;
        std::istringstream iss(s);
        if((iss >> std::dec >> d).fail())
            throw PSIEXCEPTION("Unable to convert " + s + " to a double");
        return d;
    }

    /**
     * Attempts to interpret a string as an atom specifier in a zmatrix.
     *
     * @param str: the string to interpret.
     * @param atoms: the list of atoms known so far
     * @param line: the current line, for error message printing.
     * @return the atom number (adjusted to zero-based counting)
     */
    int
    Molecule::get_anchor_atom(const std::string &str, const std::vector<std::string> &atoms,
                              const std::string &line)
    {
        if(regex_match(str, reMatches_, integerNumber_)){
            // This is just a number, return it
            return str_to_int(str) - 1;
        }else{
            // Look to see if this string is known
            std::vector<std::string>::const_iterator iter = find(atoms.begin(), atoms.end(), str);
            if(iter == atoms.end()){
                throw PSIEXCEPTION("Illegal value " + str + " in atom specification"
                                   + " on line " + line + "\n\n");
            }else{
                return distance(atoms.begin(), iter);
            }
        }
    }

    void Molecule::set_variable(const std::string &str, double val)
    {
        geometryVariables_[str] = val;

        try {
            update_geometry();
        }
        catch (...) {
            // Update geometry might have added some atoms, delete them to be safe.
            atoms_.clear();
        }
    }

    double Molecule::get_variable(const std::string &str)
    {
        if(geometryVariables_.count(str)){
            return geometryVariables_[str];
        }else{
            throw PSIEXCEPTION(str + " not known");
        }
    }

    bool Molecule::is_variable(const std::string &str) const
    {
        smatch reMatches;
        std::vector<std::string> splitLine;

        string neg = "-" + str;

        std::vector<std::string>::const_iterator line = geometryString_.begin();
        for(; line != geometryString_.end(); ++line){
            boost::split(splitLine, *line, boost::is_any_of("\t ,"),token_compress_on);

            vector<string>::iterator word = splitLine.begin();
            for (; word != splitLine.end(); ++word) {
                if (*word == str || *word == neg)
                    return true;
            }
        }

        return false;
    }

    /**
     * Attempts to interpret a string as a double.  If it's not succesful, it checks to
     * see if it's been defined as a variable and returns its value, or throws.
     *
     * @param str: the string to interpret.
     * @param line: the current line, for error message printing.
     * @return the numerical interpretation of the string.
     */
    double
    Molecule::get_value(const std::string &str, const std::string &line)
    {
        if(regex_match(str, reMatches_, realNumber_)){
            // This is already a number
            return str_to_double(str);
        }else{
            if(str[0] == '-'){
                // This is negative; ignore the leading '-' and return minus the value
                if(geometryVariables_.count(str.substr(1, str.size() - 1))){
                    return -1.0*geometryVariables_[str.substr(1, str.size() - 1)];
                }else{
                    throw PSIEXCEPTION("Illegal value in geometry specification: " + str
                                       + " on line " + line
                                       + ".\n Provide a number or define the variable\n");
                }
            }else{
                // This is positive; return the value using the string as-is
                if(geometryVariables_.count(str)){
                    return geometryVariables_[str];
                }else{
                    throw PSIEXCEPTION("Illegal value in geometry specification: " + str
                                       + " on line " + line
                                       + ".\n Provide a number or define the variable\n");
                }
            }
        }
    }

    void if_to_invert_axis(const Vector3& v1, int& must_invert, int& should_invert, double& maxproj)
    {
        int xyz, nzero;
        double vabs;

        maxproj = 0.0;
        must_invert = 0;
        should_invert = 0;

        nzero = 0;

        for(xyz=0; xyz<3; xyz++) {

            vabs = fabs(v1[xyz]);

            if (vabs < ZERO)
                nzero++;

            if (vabs > fabs(maxproj)) {
                maxproj = v1[xyz];
            }

        }

        if (nzero == 2) {
            if (maxproj < 0.0)
                must_invert = 1;
            else
                must_invert = 0;
        }
        else if (nzero < 2) {
            if (maxproj < 0.0)
                should_invert = 1;
            else
                should_invert = 0;
        }
    }
    extern FILE *outfile;
}

Molecule::Molecule():
    nirreps_(0),
    nunique_(0),
    nequiv_(0),
    equiv_(0),
    atom_to_unique_(0),
    multiplicity_(1),
    molecularCharge_(0),
    geometryFormat_(None)
{
    std::string def_name("default");
    name_ = def_name;
}

Molecule::Molecule(const Molecule& other)
{
    atoms_ = other.atoms_;
    nirreps_ = other.nirreps_;
    name_ = other.name_;
}

Molecule::~Molecule()
{
    clear();
}

Molecule& Molecule::operator=(const Molecule& other)
{
    // Self assignment is bad
    if (this == &other)
        return *this;

    atoms_ = other.atoms_;
    nirreps_ = other.nirreps_;

    return *this;
}

void Molecule::clear()
{
    nirreps_ = 0;
    atoms_.empty();
}

void Molecule::add_atom(int Z, double x, double y, double z,
                        const char *label, double mass, double charge)
{
    atom_info info;

    info.x = x;
    info.y = y;
    info.z = z;
    info.Z = Z;
    info.charge = charge;
    info.label = label;
    info.mass = mass;

    // Dummies go to full_atoms_, ghosts need to go to both.
    if(strcmp(label, "X") && strcmp(label, "x")) atoms_.push_back(info);
    full_atoms_.push_back(info);
}

double Molecule::mass(int atom) const
{
    if (atoms_[atom].mass != 0.0)
        return atoms_[atom].mass;

    return an2masses[atoms_[atom].Z];
}

const std::string Molecule::label(int atom) const
{
    return atoms_[atom].label;
}

int Molecule::atom_at_position1(double *xyz, double tol) const
{
    Vector3 b(xyz);
    for (int i=0; i < natom(); ++i) {
        Vector3 a(atoms_[i].x, atoms_[i].y, atoms_[i].z);
        if (b.distance(a) < tol)
            return i;
    }
    return -1;
}

int Molecule::atom_at_position2(Vector3& b, double tol) const
{
    for (int i=0; i < natom(); ++i) {
        Vector3 a(atoms_[i].x, atoms_[i].y, atoms_[i].z);
        if (b.distance(a) < tol)
            return i;
    }
    return -1;
}

Vector3 Molecule::center_of_mass() const
{
    Vector3 ret;
    double total_m;

    ret = 0.0;
    total_m = 0.0;

    for (int i=0; i<natom(); ++i) {
        double m = mass(i);
        ret += m * xyz(i);
        total_m += m;
    }

    ret *= 1.0/total_m;

    return ret;
}

double Molecule::nuclear_repulsion_energy()
{
    double e=0.0;

    for (int i=1; i<natom(); ++i) {
        for (int j=0; j<i; ++j) {
            e += Z(i) * Z(j) / (xyz(i).distance(xyz(j)));
        }
    }

    return e;
}

SimpleMatrix Molecule::nuclear_repulsion_energy_deriv1()
{
    SimpleMatrix de("Nuclear Repulsion Energy 1st Derivatives", natom(), 3);

    for (int i=0; i<natom(); ++i) {
        for (int j=0; j<natom(); ++j) {
            if (i != j) {
                double temp = pow((xyz(i).distance(xyz(j))), 3.0);
                de(i, 0) -= (x(i) - x(j)) * Z(i) * Z(j) / temp;
                de(i, 1) -= (y(i) - y(j)) * Z(i) * Z(j) / temp;
                de(i, 2) -= (z(i) - z(j)) * Z(i) * Z(j) / temp;
            }
        }
    }

    return de;
}

/*
    TODO Test nuclear_repulsion_energy_deriv2
*/
SimpleMatrix Molecule::nuclear_repulsion_energy_deriv2()
{
    SimpleMatrix hess("Nuclear Repulsion Energy 2nd Derivatives", 3*natom(), 3*natom());
    double sx, sy, sz, x2, y2, z2, r2, r, r5, pfac;

    for (int i=1; i<natom(); ++i) {
        int ix = 3*i;
        int iy = ix+1;
        int iz = iy+1;

        for (int j=0; j<i; ++j) {
            int jx = 3*j;
            int jy = jx+1;
            int jz = jy+1;

            sx = x(i) - x(j);
            sy = y(i) - y(j);
            sz = z(i) - z(j);

            x2 = sx*sx; y2 = sy*sy; z2 = sz*sz;
            r2 = x2 + y2 + z2;
            r = sqrt(r2);
            r5 = r2*r2*r;
            pfac = Z(i) * Z(j) / r5;

            hess.add(ix, ix, pfac * (3*x2 - r2));
            hess.add(iy, iy, pfac * (3*y2 - r2));
            hess.add(iz, iz, pfac * (3*z2 - r2));
            hess.add(ix, iy, pfac*3*sx*sy);
            hess.add(ix, iz, pfac*3*sx*sz);
            hess.add(iy, iz, pfac*3*sy*sz);

            hess.add(jx, jx, pfac * (3*x2 - r2));
            hess.add(jy, jy, pfac * (3*y2 - r2));
            hess.add(jz, jz, pfac * (3*z2 - r2));
            hess.add(jx, jy, pfac*3*sx*sy);
            hess.add(jx, jz, pfac*3*sx*sz);
            hess.add(jy, jz, pfac*3*sy*sz);

            hess.add(ix, jx, -pfac*(3*sx*sx-r2));
            hess.add(ix, jy, -pfac*(3*sx*sy));
            hess.add(ix, jz, -pfac*(3*sx*sz));
            hess.add(iy, jx, -pfac*(3*sy*sx));
            hess.add(iy, jy, -pfac*(3*sy*sy-r2));
            hess.add(iy, jz, -pfac*3*sy*sz);
            hess.add(iz, jx, -pfac*3*sz*sx);
            hess.add(iz, jy, -pfac*3*sz*sy);
            hess.add(iz, jz, -pfac*(3*sz*sz-r2));
        }
    }

    hess.element_add_mirror();

    return hess;
}

void Molecule::translate(const Vector3& r)
{
    for (int i=0; i<natom(); ++i) {
        atoms_[i].x += r[0];
        atoms_[i].y += r[1];
        atoms_[i].z += r[2];
    }
    for (int i=0; i<nallatom(); ++i) {
        full_atoms_[i].x += r[0];
        full_atoms_[i].y += r[1];
        full_atoms_[i].z += r[2];
    }
}

void Molecule::move_to_com()
{
    Vector3 com = -center_of_mass();
    translate(com);
}

SimpleMatrix Molecule::geometry()
{
    SimpleMatrix geom(natom(), 3);
    for (int i=0; i<natom(); ++i) {
        geom[i][0] = x(i);
        geom[i][1] = y(i);
        geom[i][2] = z(i);
    }

    return geom;
}
SimpleMatrix Molecule::full_geometry()
{
    SimpleMatrix geom(nallatom(), 3);
    for (int i=0; i<nallatom(); ++i) {
        geom[i][0] = fx(i);
        geom[i][1] = fy(i);
        geom[i][2] = fz(i);
    }

    return geom;
}

void Molecule::set_geometry(SimpleMatrix& geom)
{
    for (int i=0; i<natom(); ++i) {
        atoms_[i].x = geom[i][0];
        atoms_[i].y = geom[i][1];
        atoms_[i].z = geom[i][2];
    }
}
void Molecule::set_full_geometry(SimpleMatrix& geom)
{
    for (int i=0; i<nallatom(); ++i) {
        full_atoms_[i].x = geom[i][0];
        full_atoms_[i].y = geom[i][1];
        full_atoms_[i].z = geom[i][2];
    }
}

void Molecule::rotate(SimpleMatrix& R)
{
    SimpleMatrix new_geom(natom(), 3);
    SimpleMatrix geom = geometry();

    // Multiple the geometry by the rotation matrix.
    new_geom.gemm(false, false, 1.0, &geom, &R, 0.0);

    set_geometry(new_geom);
}
void Molecule::rotate_full(SimpleMatrix& R)
{
    SimpleMatrix new_geom(nallatom(), 3);
    SimpleMatrix geom = full_geometry();

    // Multiply the geometry by the rotation matrix.
    new_geom.gemm(false, false, 1.0, &geom, &R, 0.0);

    set_full_geometry(new_geom);
}

void Molecule::reorient()
{
    // Nothing for us to do.
    if (natom() <= 1)
        return;

    // Otherwise, do something.
    // Retrieve the inertia tensor.
    SimpleMatrix *itensor = inertia_tensor();
    SimpleMatrix itensor_axes(3, 3);
    SimpleVector itensor_moments(3);

    // Diagonalize the tensor matrix
    itensor->diagonalize(&itensor_axes, &itensor_moments);

    // Locate degeneracies
    int degen=0, deg_IM1=0, deg_IM2=0;
    int i, j;
    double abs, rel;
    for (i=0; i<2; ++i) {
        for (j=i+1; j<3; ++j) {
            abs = fabs(itensor_moments[i] - itensor_moments[j]);
            double tmp = (itensor_moments[i] > itensor_moments[j]) ? itensor_moments[i] : itensor_moments[j];
            if (abs > 1.0e-14)
                rel = abs / tmp;
            else
                rel = 0.0;
            if (rel < ZERO_MOMENT_INERTIA) {
                degen++;
                deg_IM1 = i;
                deg_IM2 = j;
            }
        }
    }

    Vector3 v1(itensor_axes(0, 1), itensor_axes(1, 1), itensor_axes(2, 1));
    Vector3 v2(itensor_axes(0, 2), itensor_axes(1, 2), itensor_axes(2, 2));
    Vector3 v3 = v1.cross(v2);
    itensor_axes(0, 0) = v3[0];
    itensor_axes(1, 0) = v3[1];
    itensor_axes(2, 0) = v3[2];

    int nmust = 0, nshould = 0, must_invert[3], should_invert[3];
    int axis;
    double maxproj[3];
    for (axis=0; axis<3; ++axis) {
        v1[0] = itensor_axes(0, axis);
        v1[1] = itensor_axes(1, axis);
        v1[2] = itensor_axes(2, axis);

        if_to_invert_axis(v1, must_invert[axis], should_invert[axis], maxproj[axis]);
        nmust += must_invert[axis];
        nshould += should_invert[axis];
    }

    SimpleMatrix R(3, 3);
    if (nmust == 2) {
        for (axis=0; axis<3; ++axis) {
            if (must_invert[axis])
                R[axis][axis] = -1.0;
            else
                R[axis][axis] = 1.0;
        }
    }
    else if (nmust == 1 && nshould > 0) {
        int axis1, axis2;
        if (nshould == 2) {
            for (axis=0; axis<3; ++axis) {
                if (should_invert[axis]) {
                    axis1 = axis;
                    axis++;
                    break;
                }
            }
            for (; axis<3; ++axis) {
                if (should_invert[axis]) {
                    axis2 = axis;
                    break;
                }
            }
            if (fabs(maxproj[axis1]) > fabs(maxproj[axis2])) {
                nshould = 1;
                should_invert[axis2] = 0;
            }
            else {
                nshould = 1;
                should_invert[axis1] = 0;
            }
        }
        for (axis=0; axis<3; ++axis) {
            if (must_invert[axis])
                R[axis][axis] = -1.0;
            else if (should_invert[axis])
                R[axis][axis] = -1.0;
            else
                R[axis][axis] = 1.0;
        }
    }
    else if (nmust == 3) {
        R[0][0] = -1.0;
        R[1][1] = -1.0;
        R[2][2] = 1.0;
    }
    else if (nmust == 0 && nshould > 1) {
        if (nshould == 3) {
            double tmp = fabs(maxproj[0]);
            i=0;
            for (axis=1; axis<3; ++axis) {
                if (fabs(maxproj[axis]) < fabs(tmp)) {
                    i = axis;
                    tmp = fabs(maxproj[axis]);
                }
            }
            should_invert[i] = 0;
            nshould = 2;
        }
        for (axis=0; axis<3; ++axis) {
            if (should_invert[axis])
                R[axis][axis] = -1.0;
            else
                R[axis][axis] = 1.0;
        }
    }
    else {
        R[0][0] = 1.0;
        R[1][1] = 1.0;
        R[2][2] = 1.0;
    }

    if (degen == 0) {
        rotate(itensor_axes);
        rotate_full(itensor_axes);
        rotate(R);
        rotate_full(R);
    }

    if (degen == 1) {
        int must_invert, should_invert, unique_axis;
        double maxproj, invert_pfac;

        if (deg_IM1 + deg_IM2 == 3)
            unique_axis = 0;
        else
            unique_axis = 2;

        v1[0] = itensor_axes[0][unique_axis];
        v1[1] = itensor_axes[1][unique_axis];
        v1[2] = itensor_axes[2][unique_axis];

        if_to_invert_axis(v1, must_invert, should_invert, maxproj);
        if (must_invert || should_invert)
            invert_pfac = 1.0;
        else
            invert_pfac = -1.0;

        v1 *= invert_pfac;

        double cos_theta = v1[2];
        double theta, sin_theta, v2norm, cos_phix, cos_phiy, phix;
        double sin_phix;
        if ( (1.0 - fabs(cos_theta)) > ZERO_MOMENT_INERTIA) {
            theta = acos(cos_theta);
            sin_theta = sin(theta);

            v3[0] = 0.0; v3[1] = 0.0; v3[2] = 1.0;
            v2 = v1.cross(v3);
            v2.normalize();

            cos_phix = v2[0];
            cos_phiy = v2[1];
            phix = acos(cos_phix);

            if (cos_phiy > 0.0) {
                phix *= -1.0;
            }
            sin_phix = sin(phix);

            R.zero();
            R[2][2] = 1.0;
            R[0][0] = cos_phix;
            R[1][1] = cos_phix;
            R[0][1] = sin_phix;
            R[1][0] = -sin_phix;
            rotate(R);
            rotate_full(R);

            R.zero();
            R[0][0] = 1.0;
            R[1][1] = cos_theta;
            R[2][2] = cos_theta;
            R[1][2] = sin_theta;
            R[2][1] = -sin_theta;
            rotate(R);
            rotate_full(R);

            R.zero();
            R[2][2] = 1.0;
            R[0][0] = cos_phix;
            R[1][1] = cos_phix;
            R[0][1] = -sin_phix;
            R[1][0] = sin_phix;
            rotate(R);
            rotate_full(R);
        }
    }

    // Delete the tensor matrix
    delete itensor;
}

int Molecule::nfrozen_core(std::string depth)
{
    if (depth == "FALSE") {
        return 0;
    }
    else if (depth == "TRUE" || depth == "SMALL") {
        int nfzc = 0;
        for (int A = 0; A < atoms_.size(); A++) {
            if (atoms_[A].Z > 2 && atoms_[A].Z <= 10)
                nfzc++;
            else if (atoms_[A].Z > 10)
                nfzc+=2;
        }
        return nfzc;
    }
    else if (depth == "LARGE") {
        int nfzc = 0;
        for (int A = 0; A < atoms_.size(); A++) {
            if (atoms_[A].Z > 2 && atoms_[A].Z <= 10)
                nfzc++;
            else if (atoms_[A].Z > 10)
                nfzc+=5;
        }
        return nfzc;
    }
    else {
        throw std::invalid_argument("Frozen core spec is not supported, options are {true, false, small, large}.");
    }
}

void Molecule::init_with_psio(shared_ptr<PSIO> psio)
{
    // User sent a psio object. Create a chkpt object based on it.
    shared_ptr<Chkpt> chkpt(new Chkpt(psio.get(), PSIO_OPEN_OLD));
    init_with_chkpt(chkpt);
}

void Molecule::init_with_chkpt(shared_ptr<Chkpt> chkpt)
{
    int atoms = chkpt->rd_natom();
    double *zvals = chkpt->rd_zvals();
    double **geom = chkpt->rd_geom();

    for (int i=0; i<atoms; ++i) {
        //fprintf(outfile,"  Atom %d, Z = %d, x = %14.10f,%14.10f,%14.10f, Label , Mass, \n",i,(int)zvals[i],geom[i][0],geom[i][1],geom[i][2]); fflush(outfile);
        add_atom((int)zvals[i], geom[i][0], geom[i][1], geom[i][2], atomic_labels[(int)zvals[i]], an2masses[(int)zvals[i]]);
    }

    nirreps_ = chkpt->rd_nirreps();

    Chkpt::free(zvals);
    Chkpt::free(geom);
}

void Molecule::init_with_xyz(const std::string& xyzfilename)
{
    Element_to_Z Z;
    Z.load_values();

    if (xyzfilename.empty())
        throw PSIEXCEPTION("Molecule::init_with_xyz: given filename is blank.");

    //printf("%s",xyzfilename.c_str());

    ifstream infile(xyzfilename.c_str());
    string line, natom_str;
    const string bohr("bohr"), au("au");
    bool angstrom_in_file = true;

    if (!infile)
        throw PSIEXCEPTION("Molecule::init_with_xyz: Unable to open xyz file.");

    // Read in first line
    getline(infile, line);

    // This is what we should match on the first line
    boost::regex rx("(\\d+)\\s*(bohr|au)?", boost::regbase::normal | boost::regbase::icase);
    boost::smatch what;

    int natom;
    // Try to match the first line
    if (regex_match(line, what, rx)) {
        // matched
        // Convert the matches to what we need.
        if (!from_string<int>(natom, what[1], std::dec))
            throw PSIEXCEPTION("Molecule::init_with_xyz: Unable to convert number of atoms from xyz file.");

//        cout << "init_with_xyz: " << what.size() << endl;
        if (what.size() == 3) {
            string s(what[2].first, what[2].second);
            if (boost::iequals(bohr, s) || boost::iequals(au, s)) {
                angstrom_in_file = false;
            }
        }
    }
    else
        throw PSIEXCEPTION("Molecule::init_with_xyz: Malformed first line\n"+line);

    // Next line is a comment line, ignore it
    getline(infile, line);

    // Next line begins the useful information.
    // This is the regex for the remaining lines
    rx.assign("(?:\\s*)([A-Z](?:[a-z])?)(?:\\s+)(-?\\d+\\.\\d+)(?:\\s+)(-?\\d+\\.\\d+)(?:\\s+)(-?\\d+\\.\\d+)(?:\\s*)",
        boost::regbase::normal | boost::regbase::icase);
    for (int i=0; i<natom; ++i) {
        // Get an atom info line.
        getline(infile, line);

        // Try to match it
        if (regex_match(line, what, rx))
        {
            // First is a string
            string atomSym(what[1].first, what[1].second);
            transform(atomSym.begin(), atomSym.end(), atomSym.begin(), ::toupper);

            // Then the coordinates:
            double x, y, z;
            if (!from_string<double>(x, what[2], std::dec))
                throw PSIEXCEPTION("Molecule::init_with_xyz: Unable to convert x coordinate.\n" + line);
            if (!from_string<double>(y, what[3], std::dec))
                throw PSIEXCEPTION("Molecule::init_with_xyz: Unable to convert y coordinate.\n" + line);
            if (!from_string<double>(z, what[4], std::dec))
                throw PSIEXCEPTION("Molecule::init_with_xyz: Unable to convert z coordinate.\n" + line);

            if (angstrom_in_file) {
                // Coordinates in Molecule must be bohr.
                x /= _bohr2angstroms;
                y /= _bohr2angstroms;
                z /= _bohr2angstroms;
            }

            // Add it to the molecule.
            add_atom((int)Z[atomSym], x, y, z, atomSym.c_str(), atomic_masses[(int)Z[atomSym]]);
        }
        else {
            throw PSIEXCEPTION("Molecule::init_with_xyz: Malformed atom information line.\n"+line);
        }
    }
}

/**
 * Given a string (including newlines to separate lines), builds a new molecule
 * and wraps it in a smart pointer
 *
 * @param text: a string providing the user's input
 */
shared_ptr<Molecule>
Molecule::create_molecule_from_string(const std::string &text)
{
    smatch reMatches;
    // Split the input at newlines, storing the result in "lines"
    std::vector<std::string> lines;
    boost::split(lines, text, boost::is_any_of("\n"));

    shared_ptr<Molecule> mol(new Molecule);
    mol->set_units(Angstrom);
    /*
     * Time to look for lines that look like they describe charge and multiplicity,
     * a variable, units, comment lines, and blank lines.  When found, process them
     * and remove them so that only the raw geometry remains.  Iterated backwards, as
     * elements are deleted as they are processed.
     */
    for(int lineNumber = lines.size() - 1 ; lineNumber >= 0; --lineNumber){
        if(regex_match(lines[lineNumber], reMatches, variableDefinition_)){
            // A variable definition
            double value = (boost::iequals("tda", reMatches[2].str()) ?
                               360.0*atan(sqrt(2))/M_PI : str_to_double(reMatches[2]));
            mol->set_variable(reMatches[1].str(), value);
            lines.erase(lines.begin() + lineNumber);
        }else if(regex_match(lines[lineNumber], reMatches, blankLine_)){
            // A blank line, nuke it
            lines.erase(lines.begin() + lineNumber);
        }else if(regex_match(lines[lineNumber], reMatches, commentLine_)){
            // A comment line, nuke it
            lines.erase(lines.begin() + lineNumber);
        }else if(regex_match(lines[lineNumber], reMatches, unitLabel_)){
            // A comment line, nuke it
            if(   boost::iequals("ang", reMatches[1].str())
               || boost::iequals("angstrom",   reMatches[1].str())){
                mol->set_units(Angstrom);
            }else{
                mol->set_units(Bohr);
            }
            lines.erase(lines.begin() + lineNumber);
        }else if(regex_match(lines[lineNumber], reMatches, chargeAndMultiplicity_)){
            if(lineNumber && !regex_match(lines[lineNumber-1], reMatches, fragmentMarker_)){
                // As long as this does not follow a "--", it's a global charge/multiplicity
                // specifier, so we process it, then nuke it
                mol->set_molecular_charge(str_to_int(reMatches[1]));
                mol->set_multiplicity(str_to_int(reMatches[2]));
                lines.erase(lines.begin() + lineNumber);
            }
        }
    }

    // Now go through the rest of the lines looking for fragment markers
    unsigned int firstAtom  = 0;
    mol->fragmentMultiplicities_.push_back(mol->multiplicity_);
    mol->fragmentCharges_.push_back(mol->molecularCharge_);
    unsigned int atomCount = 0;
    for(unsigned int lineNumber = 0; lineNumber < lines.size(); ++lineNumber){
        if(regex_match(lines[lineNumber], reMatches, fragmentMarker_)){
            // Check that there are more lines remaining
            if(lineNumber == lines.size() - 1)
                throw PSIEXCEPTION("Nothing specified after the final \"--\" in geometry");
            // Now we process the atom markers
            mol->fragments_.push_back(std::make_pair(firstAtom, atomCount));
            firstAtom = atomCount;
            // Figure out how to handle the multiplicity
            if(regex_match(lines[lineNumber+1], reMatches, chargeAndMultiplicity_)){
                // The user specified a charge/multiplicity for this fragment
                mol->fragmentCharges_.push_back(str_to_int(reMatches[1]));
                mol->fragmentMultiplicities_.push_back(str_to_int(reMatches[2]));
                // Don't forget to skip over the charge multiplicity line..
                ++lineNumber;
            }else{
                // The user didn't give us charge/multiplicity - use the molecule default
                mol->fragmentCharges_.push_back(mol->molecularCharge_);
                mol->fragmentMultiplicities_.push_back(mol->multiplicity_);
            }
        }else{
            ++atomCount;
        }
    }
    mol->fragments_.push_back(std::make_pair(firstAtom, atomCount));

    // Clean up the "--" and charge/multiplicity specifiers - they're no longer needed
    for(int lineNumber = lines.size() - 1 ; lineNumber >= 0; --lineNumber){
        if(   regex_match(lines[lineNumber], reMatches, fragmentMarker_)
           || regex_match(lines[lineNumber], reMatches, chargeAndMultiplicity_))
           lines.erase(lines.begin() + lineNumber);
    }

//    for(int i = 0; i < mol->fragments_.size(); ++i)
//        fprintf(outfile, "Fragment %d, atom %d to %d, charge %d mult %d\n",
//                i,mol->fragments_[i].first,mol->fragments_[i].second, mol->fragmentCharges_[i], mol->fragmentMultiplicities_[i]);

    if(!lines.size()) throw PSIEXCEPTION("No geometry specified");

    // For now, we can count the entries on the first line to determine the type
    std::vector<std::string> splitLine;
    // Trim leading and trailing whitespace
    boost::algorithm::trim(lines[0]);
    boost::split(splitLine, lines[0], boost::is_any_of("\t ,"),token_compress_on);
    if(splitLine.size() == 4){
        mol->set_geometry_format(Cartesian);
    }else if(splitLine.size() == 1){
        mol->set_geometry_format(ZMatrix);
    }else{
        throw PSIEXCEPTION("Illegal geometry specification line : " + lines[0] +
                       ".  You should provide either Z-Matrix or Cartesian input");
    };

    mol->set_geometry_string(lines);
//    mol->update_geometry();
    return mol;
}

/**
 * Updates the geometry, by (re)interpreting the string used to create the molecule, and the current values
 * of the variables. The atoms list is cleared, and then rebuilt by this routine.
 */
void
Molecule::update_geometry()
{
    // Nothing to do.
    if (geometryString_.empty())
        return;

    smatch reMatches;
    Element_to_Z zVals;
    zVals.load_values();
    std::vector<std::string> splitLine;
    std::vector<int> zIn;
    for(int i = 0; i < atoms_.size(); ++i) zIn.push_back(Z(i));
    
    // Start over, even if we have atoms already
    atoms_.clear();
    // Have the atoms been added already?  If so, we only need to update positions.
//    bool firstRun = natom();
    // Begin by assuming that the user provided angstrom input
    double conversionFactor = (units_ == Angstrom ? 1.0 / _bohr2angstroms : 1.0);
    // Now it's time to interpret the geometry string
    std::string atomSym;
    double x, y, z;
    unsigned int currentAtom = 0;
    if(geometryFormat_ == Cartesian){
        std::vector<std::string>::iterator line = geometryString_.begin();
        for(; line != geometryString_.end(); ++line){
            //fprintf(outfile,"GeometryString size %d\n",geometryString_.size());
            // Trim leading and trailing whitespace
            boost::algorithm::trim(*line);
            boost::split(splitLine, *line, boost::is_any_of("\t ,"),token_compress_on);
            if(splitLine.size() != 4){
                throw PSIEXCEPTION("Incorrect number of entries in geometry specification :" + *line);
            }
            // Check that the atom symbol is valid
            if(!regex_match(splitLine[0], reMatches, atomSymbol_))
                throw PSIEXCEPTION("Illegal atom symbol in geometry specification: " + splitLine[0]
                                   + " on line\n" + *(line));
            atomSym = boost::to_upper_copy(reMatches[1].str());
            x = get_value(splitLine[1], *line) * conversionFactor;
            y = get_value(splitLine[2], *line) * conversionFactor;
            z = get_value(splitLine[3], *line) * conversionFactor;
            if(zIn.empty()){
                add_atom((int)zVals[atomSym], x, y, z, atomSym.c_str(), atomic_masses[(int)zVals[atomSym]]);
            }else{
                add_atom(zIn[currentAtom], x, y, z, atomSym.c_str(), atomic_masses[(int)zVals[atomSym]]);
            }
            ++currentAtom;
        }
    }else if(geometryFormat_ == ZMatrix){
        std::vector<std::string>::iterator line = geometryString_.begin();
        double r, a, d; // The bond length, valence bond angle, and torsion angle, respectively
        int rTo, aTo, dTo; // The "anchor" atoms, used to define r, a, and d
        std::vector<std::string> atoms;
        for(; line != geometryString_.end(); ++line){
            // Trim leading and trailing whitespace
            boost::algorithm::trim(*line);
            boost::split(splitLine, *line, boost::is_any_of("\t ,"),token_compress_on);
            if(currentAtom == 0){
                // This is the first line
                if(splitLine.size() != 1){
                    throw PSIEXCEPTION("Incorrect number of entries in geometry specification :" + *line);
                }
                // Check that the atom symbol is valid
                if(!regex_match(splitLine[0], reMatches, atomSymbol_))
                    throw PSIEXCEPTION("Illegal atom symbol in geometry specification: " + splitLine[0]
                                       + " on line\n" + *(line));
                atoms.push_back(splitLine[0]);
                atomSym = boost::to_upper_copy(reMatches[1].str());
                x = y = z = 0.0;
            }else if(currentAtom == 1){
                // This is the second line
                if(splitLine.size() != 3){
                    throw PSIEXCEPTION("Incorrect number of entries in geometry specification :" + *line);
                }
                // Check that the atom symbol is valid
                if(!regex_match(splitLine[0], reMatches, atomSymbol_))
                    throw PSIEXCEPTION("Illegal atom symbol in geometry specification: " + splitLine[0]
                                       + " on line\n" + *(line));
                atomSym = boost::to_upper_copy(reMatches[1].str());
                atoms.push_back(splitLine[0]);
                rTo = get_anchor_atom(splitLine[1], atoms, *line);
                if(rTo >= currentAtom)
                     throw PSIEXCEPTION("Error on geometry input line " + *line + "\nAtom "
                                         + splitLine[1] + " has not been defined yet.");
                r = get_value(splitLine[2], *line) * conversionFactor;
                x = y = 0.0;
                z = r;
            }else if(currentAtom == 2){
                // This is the third line
                if(splitLine.size() != 5){
                    throw PSIEXCEPTION("Incorrect number of entries in geometry specification :" + *line);
                }
                // Check that the atom symbol is valid
                if(!regex_match(splitLine[0], reMatches, atomSymbol_))
                    throw PSIEXCEPTION("Illegal atom symbol in geometry specification: " + splitLine[0]
                                       + " on line\n" + *line);
                atomSym = boost::to_upper_copy(reMatches[1].str());
                atoms.push_back(splitLine[0]);
                rTo = get_anchor_atom(splitLine[1], atoms, *line);
                if(rTo >= currentAtom)
                     throw PSIEXCEPTION("Error on geometry input line " + *line + "\nAtom "
                                         + splitLine[1] + " has not been defined yet.");
                r = get_value(splitLine[2], *line) * conversionFactor;
                aTo = get_anchor_atom(splitLine[3], atoms, *line);
                if(aTo >= currentAtom)
                     throw PSIEXCEPTION("Error on geometry input line " + *line + "\nAtom "
                                         + splitLine[3] + " has not been defined yet.");
                a = get_value(splitLine[4], *line);
                if(aTo == rTo)
                     throw PSIEXCEPTION("Atom used multiple times on line " + *line);

                a *= M_PI / 180.0;
                if(rTo == 0){
                    x = r*sin(a);
                    y = 0.0;
                    z = r*cos(a);
                }else{
                    x = r*sin(a);
                    y = 0.0;
                    z = atoms_[1].z - r*cos(a);
                }
            }else{
                // This is line 4 onwards
                if(splitLine.size() != 7){
                    throw PSIEXCEPTION("Incorrect number of entries in geometry specification :" + *line);
                }
                // Check that the atom symbol is valid
                if(!regex_match(splitLine[0], reMatches, atomSymbol_))
                    throw PSIEXCEPTION("Illegal atom symbol in geometry specification: " + splitLine[0]
                                       + " on line\n" + *(line));
                atomSym = boost::to_upper_copy(reMatches[1].str());
                atoms.push_back(splitLine[0]);
                rTo = get_anchor_atom(splitLine[1], atoms, *line);
                if(rTo >= currentAtom)
                     throw PSIEXCEPTION("Error on geometry input line " + *line + "\nAtom "
                                         + splitLine[1] + " has not been defined yet.");
                r = get_value(splitLine[2], *line) * conversionFactor;
                aTo = get_anchor_atom(splitLine[3], atoms, *line);
                if(aTo >= currentAtom)
                     throw PSIEXCEPTION("Error on geometry input line " + *line + "\nAtom "
                                         + splitLine[3] + " has not been defined yet.");
                a = get_value(splitLine[4], *line);
                dTo = get_anchor_atom(splitLine[5], atoms, *line);
                if(dTo >= currentAtom)
                     throw PSIEXCEPTION("Error on geometry input line " + *line + "\nAtom "
                                         + splitLine[5] + " has not been defined yet.");
                d = get_value(splitLine[6], *line);
                if(aTo == rTo || rTo == dTo /* for you star wars fans */ || aTo == dTo)
                     throw PSIEXCEPTION("Atom used multiple times on line " + *line);

                a *= M_PI / 180.0;
                d *= M_PI / 180.0;

                /*
                 * The atom specification is
                 *      this       rTo   rVal  aTo  aVal   dTo   dVal
                 *        D         C           B           A
                 * which allows us to define the vector from B->C (eBC) as the +z axis, and eAB
                 * lies in the xz plane.  Then eX, eY and eZ (=eBC) are the x, y, and z axes, respecively.
                 */
                Vector3 eAB(atoms_[aTo].x - atoms_[dTo].x, atoms_[aTo].y - atoms_[dTo].y, atoms_[aTo].z - atoms_[dTo].z);
                Vector3 eBC(atoms_[rTo].x - atoms_[aTo].x, atoms_[rTo].y - atoms_[aTo].y, atoms_[rTo].z - atoms_[aTo].z);
                eAB.normalize();
                eBC.normalize();
                double cosBCD = cos(a);
                double sinBCD = sin(a);
                double cosABCD = cos(d);
                double sinABCD = sin(d);
                double cosABC = -eAB.dot(eBC);
                double sinABC = sqrt(1.0 - cosABC * cosABC);
                if(sinABC < LINEAR_A_TOL || sinBCD < LINEAR_A_TOL)
                    throw PSIEXCEPTION("Line " + *line + " defines a dihedral using a linear bond angle");
                Vector3 eY = eAB.perp_unit(eBC);
                Vector3 eX = eY.perp_unit(eBC);

                x = atoms_[rTo].x + r * ( -eBC[0] * cosBCD + eX[0] * sinBCD * cosABCD + eY[0] * sinBCD * sinABCD);
                y = atoms_[rTo].y + r * ( -eBC[1] * cosBCD + eX[1] * sinBCD * cosABCD + eY[1] * sinBCD * sinABCD);
                z = atoms_[rTo].z + r * ( -eBC[2] * cosBCD + eX[2] * sinBCD * cosABCD + eY[2] * sinBCD * sinABCD);
            }
            if(zIn.empty()){
                add_atom((int)zVals[atomSym], x, y, z, atomSym.c_str(), atomic_masses[(int)zVals[atomSym]]);
            }else{
                add_atom(zIn[currentAtom], x, y, z, atomSym.c_str(), atomic_masses[(int)zVals[atomSym]]);
            }
            ++currentAtom;
        }
    }

    fprintf(outfile, "moving to center of mass.\n");
    move_to_com();
    fprintf(outfile, "reorienting molecule.\n");
    reorient();
}


/**
 * A wrapper to extract_subsets, callable from Boost
 * @param reals: A list containing the real atoms.
 * @param ghost: A list containing the ghost atoms.
 * @return The ref counted cloned molecule.
 */
boost::shared_ptr<Molecule>
Molecule::py_extract_subsets(boost::python::list & reals,
                             boost::python::list & ghosts)
{
    std::vector<int> realVec;
    for(int i = 0; i < boost::python::len(reals); ++i)
        realVec.push_back(boost::python::extract<int>(reals[i] - 1));
    std::vector<int> ghostVec;
    for(int i = 0; i < boost::python::len(ghosts); ++i)
        ghostVec.push_back(boost::python::extract<int>(ghosts[i] - 1));

    return extract_subsets(realVec, ghostVec);
}


/**
 * Makes a copy of the molecule, returning a new ref counted molecule with
 * only certain fragment atoms present as either ghost or real atoms
 * @param real_list: The list of fragments that should be present in the molecule as real atoms.
 * @param ghost_list: The list of fragments that should be present in the molecule as ghosts.
 * @return The ref counted cloned molecule
 */
boost::shared_ptr<Molecule>
Molecule::extract_subsets(const std::vector<int> &real_list, const std::vector<int> &ghost_list) const
{
    if(ghost_list.size() + real_list.size() > fragments_.size())
        throw PSIEXCEPTION("The sum of real- and ghost-atom subsets is greater than the number of subsets");

    boost::shared_ptr<Molecule> clone(new Molecule);
    std::vector<std::pair<int, int> >list;
    // 0 is a ghost atom, 1 is a real atom
    for(int i = 0; i < ghost_list.size(); ++i){
        if(ghost_list[i] > fragments_.size())
            throw PSIEXCEPTION("Invalid ghost atom subset requested");
        list.push_back(std::make_pair(ghost_list[i], 0));
    }
    for(int i = 0; i < real_list.size(); ++i){
        if(real_list[i] > fragments_.size())
            throw PSIEXCEPTION("Invalid real atom subset requested");
        list.push_back(std::make_pair(real_list[i], 1));
    }
    sort(list.begin(), list.end());

    int fragCharge = 0;
    int multiplicity = 1;
    for(int fragment = 0; fragment < list.size(); ++fragment){
        int fragmentNum = list[fragment].first;
        bool isReal     = list[fragment].second;
        for(int atom = fragments_[fragmentNum].first; atom < fragments_[fragmentNum].second; ++atom){
            int Z             = full_atoms_[atom].Z;
            double x          = full_atoms_[atom].x;
            double y          = full_atoms_[atom].y;
            double z          = full_atoms_[atom].z;
            const char *label = full_atoms_[atom].label.c_str();
            double mass       = full_atoms_[atom].mass;
            double charge     = full_atoms_[atom].charge;
            if(isReal){
                clone->add_atom(Z, x, y, z, label, mass, charge);
            }else{
                clone->add_atom(0, x, y, z, label, mass, charge);
            }
            clone->geometryString_.push_back(geometryString_[atom]);
        }
        if (isReal) {
            fragCharge += fragmentCharges_[fragmentNum];
            multiplicity += fragmentMultiplicities_[fragmentNum] - 1;
        }
    }
    clone->set_molecular_charge(fragCharge);
    clone->set_multiplicity(multiplicity);
    clone->geometryFormat_ = geometryFormat_;
    clone->units_ = units_;

    return clone;
}


void Molecule::save_to_chkpt(shared_ptr<Chkpt> chkpt, std::string prefix)
{
    // Save the current prefix
    string pre = chkpt->get_prefix();
    // If needed switch the prefix in the chkpt file.
    if (!prefix.empty()) {
        chkpt->set_prefix(prefix.c_str());
    }

    // Need to save natom, zvals, geom
    chkpt->wt_natom(natom());
    chkpt->wt_nallatom(nallatom());

    double *zvals = new double[natom()];
    double **geom = block_matrix(natom(), 3);
    double **fgeom = block_matrix(nallatom(), 3);
    int *dummyflags = new int[nallatom()];

    for (int i=0; i<natom(); ++i) {
        zvals[i] = static_cast<double>(Z(i));
        geom[i][0] = x(i); geom[i][1] = y(i); geom[i][2] = z(i);
    }

    for (int i=0; i<nallatom(); ++i) {
    fgeom[i][0] = fx(i); geom[i][1] = fy(i); geom[i][2] = fz(i);
    dummyflags[i] = fZ(i) > 0 ? 0 : 1;
    }

    chkpt->wt_zvals(zvals);
    chkpt->wt_atom_dummy(dummyflags);
    chkpt->wt_fgeom(fgeom);

    chkpt->wt_enuc(nuclear_repulsion_energy());

    // Reset the prefix
    if (!prefix.empty()) {
        chkpt->set_prefix(pre.c_str());
    }

    delete[]dummyflags;
    delete[]zvals;
    free_block(geom);
    free_block(fgeom);
}

void Molecule::print() const
{
    if (natom()) {
        fprintf(outfile,"       Center              X                  Y                   Z           Charge\n");
        fprintf(outfile,"    ------------   -----------------  -----------------  -----------------  ---------\n");

        for(int i = 0; i < natom(); ++i){
            Vector3 geom = xyz(i);
            fprintf(outfile, "    %8s%4s ",label(i).c_str(),Z(i) ? "" : "(Gh)"); fflush(outfile);
            for(int j = 0; j < 3; j++)
                fprintf(outfile, "  %17.12f", geom[j]);
            fprintf(outfile,"  %17d",Z(i));
            fprintf(outfile,"\n");
        }
        fprintf(outfile,"\n");
        fflush(outfile);

        // Print symmetry information, if available
        if (nunique_) {
            fprintf(outfile, "    Number of unique atoms: %d\n\n", nunique_);
            fprintf(outfile, "    Atoms equivalency:\n");
            for (int i=0; i<nunique_; ++i) {
                fprintf(outfile, "       unique atom %d: ", i);
                for (int j=0; j<nequiv_[i]; ++j) {
                    fprintf(outfile, "%d ", equiv_[i][j]);
                }
                fprintf(outfile, "\n");
            }
            fprintf(outfile, "\n");
            fflush(outfile);
        }
    }
    else {
        fprintf(outfile, "  No atoms in this molecule.\n");
    }
}

SimpleVector Molecule::nuclear_dipole_contribution()
{
    SimpleVector ret(3);

    for(int i=0; i<natom(); ++i) {
        Vector3 geom = xyz(i);
        ret[0] += Z(i) * geom[0];
        ret[1] += Z(i) * geom[1];
        ret[2] += Z(i) * geom[2];
    }

    return ret;
}

SimpleVector Molecule::nuclear_quadrupole_contribution()
{
    SimpleVector ret(6);
    double xx, xy, xz, yy, yz, zz;

    xx = xy = xz = yy = yz = zz = 0.0;

    for (int i=0; i<natom(); ++i) {
        Vector3 geom = xyz(i);
        ret[0] += Z(i) * geom[0] * geom[0]; // xx
        ret[1] += Z(i) * geom[0] * geom[1]; // xy
        ret[2] += Z(i) * geom[0] * geom[2]; // xz
        ret[3] += Z(i) * geom[1] * geom[1]; // yy
        ret[4] += Z(i) * geom[1] * geom[2]; // yz
        ret[5] += Z(i) * geom[2] * geom[2]; // zz
    }

    return ret;
}

SimpleMatrix* Molecule::inertia_tensor()
{
    int i;
    SimpleMatrix* tensor = new SimpleMatrix("Inertia Tensor", 3, 3);

    for (i = 0; i < natom(); i++) {
        // I(alpha, alpha)
        tensor->add(0, 0, mass(i) * (pow(y(i), 2) + pow(z(i), 2)));
        tensor->add(1, 1, mass(i) * (pow(x(i), 2) + pow(z(i), 2)));
        tensor->add(2, 2, mass(i) * (pow(x(i), 2) + pow(y(i), 2)));

        // I(alpha, beta)
        tensor->add(0, 1, -mass(i) * x(i) * y(i));
        tensor->add(0, 2, -mass(i) * x(i) * z(i));
        tensor->add(1, 2, -mass(i) * y(i) * z(i));
        //    mirror
        tensor->add(1, 0, -mass(i) * x(i) * y(i));
        tensor->add(2, 0, -mass(i) * x(i) * z(i));
        tensor->add(2, 1, -mass(i) * y(i) * z(i));
    }

    return tensor;
}

double& Molecule::xyz(int atom, int _xyz)
{
    if (_xyz == 0)
        return atoms_[atom].x;
    else if (_xyz == 1)
        return atoms_[atom].y;
    else
        return atoms_[atom].z;
}

const double& Molecule::xyz(int atom, int _xyz) const
{
    if (_xyz == 0)
        return atoms_[atom].x;
    else if (_xyz == 1)
        return atoms_[atom].y;
    else
        return atoms_[atom].z;
}

//
// Symmetry
//
bool Molecule::has_inversion(Vector3& origin, double tol) const
{
    for (int i=0; i<natom(); ++i) {
        Vector3 inverted = origin-(xyz(i) - origin);
        int atom = atom_at_position2(inverted, tol);
        if (atom < 0 || Z(atom) != Z(i)) {
            return false;
        }
    }
    return true;
}

bool Molecule::is_plane(Vector3& origin, Vector3& uperp, double tol) const
{
    for (int i=0; i<natom(); ++i) {
        Vector3 A = xyz(i)-origin;
        Vector3 Apar = uperp.dot(A)*origin;
        Vector3 Aperp = A - Apar;
        A = (Aperp- Apar) + origin;
        int atom = atom_at_position2(A, tol);
        if (atom < 0 || Z(atom) != Z(i)) {
            return false;
        }
    }
    return true;
}

bool Molecule::is_axis(Vector3& origin, Vector3& axis, int order, double tol) const
{
    for (int i=0; i<natom(); ++i) {
        Vector3 A = xyz(i) - origin;
        for (int j=1; j<order; ++j) {
            Vector3 R = A;
            R.rotate(j*2.0*M_PI/order, axis);
            R += origin;
            int atom = atom_at_position2(R, tol);
            if (atom < 0 || Z(atom) != Z(i)) {
                return false;
            }
        }
    }
    return true;
}

enum AxisName { XAxis, YAxis, ZAxis };

static AxisName like_world_axis(Vector3& axis, const Vector3& worldxaxis, const Vector3& worldyaxis, const Vector3& worldzaxis)
{
    AxisName like;
    double xlikeness = fabs(axis.dot(worldxaxis));
    double ylikeness = fabs(axis.dot(worldyaxis));
    double zlikeness = fabs(axis.dot(worldzaxis));
    if (xlikeness > ylikeness && xlikeness > zlikeness) {
        like = XAxis;
        if (axis.dot(worldxaxis) < 0) axis = - axis;
    }
    else if (ylikeness > zlikeness) {
        like = YAxis;
        if (axis.dot(worldyaxis) < 0) axis = - axis;
    }
    else {
        like = ZAxis;
        if (axis.dot(worldzaxis) < 0) axis = - axis;
    }
    return like;
}

void Molecule::is_linear_planar(bool& linear, bool& planar, double tol) const
{
    if (natom() < 3) {
        linear = true;
        planar = true;
        return;
    }

    // find three atoms not on the same line
    Vector3 A = xyz(0);
    Vector3 B = xyz(1);
    Vector3 BA = B-A;
    BA.normalize();
    Vector3 CA;

    int i;
    double min_BAdotCA = 1.0;
    for (i=2; i<natom(); ++i) {
        Vector3 tmp = xyz(i) - A;
        tmp.normalize();
        if (fabs(BA.dot(tmp)) < min_BAdotCA) {
            CA = tmp;
            min_BAdotCA = fabs(BA.dot(tmp));
        }
    }
    if (min_BAdotCA >= 1.0 - tol) {
        linear = true;
        planar = true;
        return;
    }

    linear = false;
    if (natom() < 4) {
        planar = true;
        return;
    }

    // check for nontrivial planar molecules
    Vector3 BAxCA = BA.cross(CA);
    BAxCA.normalize();
    for (i=2; i<natom(); ++i) {
        Vector3 tmp = xyz(i)-A;
        if (fabs(tmp.dot(BAxCA)) > tol) {
            planar = false;
            return;
        }
    }
    planar = true;
}

boost::shared_ptr<PointGroup> Molecule::find_point_group(double tol) const
{
    int i, j;

    Vector3 com = center_of_mass();

    Vector3 worldxaxis(1.0, 0.0, 0.0);
    Vector3 worldyaxis(0.0, 1.0, 0.0);
    Vector3 worldzaxis(0.0, 0.0, 1.0);

    bool linear, planar;
    is_linear_planar(linear, planar, tol);

    bool have_inversion = has_inversion(com, tol);

    // check for C2 axis
    Vector3 c2axis;
    bool have_c2axis = false;
    if (natom() < 2) {
        have_c2axis = true;
        c2axis = Vector3(0.0, 0.0, 1.0);
    }
    else if (linear) {
        have_c2axis = true;
        c2axis = xyz(1) - xyz(0);
        c2axis.normalize();
    }
    else if (planar && have_inversion) {
        // there is a c2 axis that won't be found using the usual
        // algorithm. fine two noncolinear atom-atom vectors (we know
        // that linear == 0)
        Vector3 BA = xyz(1) - xyz(0);
        BA.normalize();
        for (i=2; i<natom(); ++i) {
            Vector3 CA = xyz(i) - xyz(0);
            CA.normalize();
            Vector3 BAxCA = BA.cross(CA);
            if (BAxCA.norm() > tol) {
                have_c2axis = true;
                BAxCA.normalize();
                c2axis = BAxCA;
                break;
            }
        }
    }
    else {
        // loop through pairs of atoms o find c2 axis candidates
        for (i=0; i<natom(); ++i) {
            Vector3 A = xyz(i) - com;
            double AdotA = A.dot(A);
            for (j=0; j<=i; ++j) {
                // the atoms must be identical
                if (Z(i) != Z(j)) continue;
                Vector3 B = xyz(j)-com;
                // the atoms must be the same distance from the com
                if (fabs(AdotA - B.dot(B)) > tol) continue;
                Vector3 axis = A+B;
                // atoms colinear with the com don't work
                if (axis.norm() < tol) continue;
                axis.normalize();
                if (is_axis(com, axis, 2, tol)) {
                    have_c2axis = true;
                    c2axis = axis;
                    goto found_c2axis;
                }
            }
        }
    }
found_c2axis:

    AxisName c2like = ZAxis;
    if (have_c2axis) {
        // try to make the sign of the axis correspond to one of the
        // world axes
        c2like = like_world_axis(c2axis, worldxaxis, worldyaxis, worldzaxis);
    }

    // check for c2 axis perp to first c2 axis
    Vector3 c2axisperp;
    bool have_c2axisperp = false;
    if (have_c2axis) {
        if (natom() < 2) {
            have_c2axisperp = true;
            c2axisperp = Vector3(1.0, 0.0, 0.0);
        }
        else if (linear) {
            if (have_inversion) {
                have_c2axisperp = true;
                c2axisperp = c2axis.perp_unit(Vector3(0.0,0.0,1.0));
            }
        }
        else {
            // loop through paris of atoms to find c2 axis candidates
            for (i=0; i<natom(); ++i) {
                Vector3 A = xyz(i) - com;
                double AdotA = A.dot(A);
                for (j=0; j<i; ++j) {
                    // the atoms must be identical
                    if (Z(i) != Z(j) || fabs(mass(i) - mass(j)) > tol) continue;
                    Vector3 B = xyz(i) - com;
                    // the atoms must be the same distance from the com
                    if (fabs(AdotA - B.dot(B)) > tol) continue;
                    Vector3 axis= A+B;
                    // atoms colinear with the com don't work
                    if (axis.norm() < tol) continue;
                    axis.normalize();
                    // if axis is not perp continue
                    if (fabs(axis.dot(c2axis)) > tol) continue;
                    if (is_axis(com, axis, 2, tol)) {
                        have_c2axisperp = true;
                        c2axisperp = axis;
                        goto found_c2axisperp;
                    }
                }
            }
        }
    }
found_c2axisperp:

    AxisName c2perplike;
    if (have_c2axisperp) {
        // try to make the sign of the axis correspond to one of
        // the world axes
        c2perplike = like_world_axis(c2axisperp, worldxaxis, worldyaxis, worldzaxis);

        // try to make c2axis the z axis
        if (c2perplike == ZAxis) {
            Vector3 tmpv = c2axisperp;
            tmpv = c2axisperp; c2axisperp = c2axis; c2axis = tmpv;
            c2perplike = c2like;
            c2like = ZAxis;
        }
        if (c2like != ZAxis) {
            if (c2like == XAxis) c2axis = c2axis.cross(c2axisperp);
            else c2axis = c2axisperp.cross(c2axis);
            c2like = like_world_axis(c2axis, worldxaxis, worldyaxis, worldzaxis);
        }
        // try to make c2axisperplike the x axis
        if (c2perplike == YAxis) {
            c2axisperp = c2axisperp.cross(c2axis);
            c2perplike = like_world_axis(c2axisperp, worldxaxis, worldyaxis, worldzaxis);
        }
    }

    // Check for vertical plane
    bool have_sigmav = false;
    Vector3 sigmav;
    if (have_c2axis) {
        if (natom() < 2) {
            have_sigmav = true;
            sigmav = c2axisperp;
        }
        else if (linear) {
            have_sigmav = true;
            if (have_c2axisperp) {
                sigmav = c2axisperp;
            }
            else {
                sigmav = c2axis.perp_unit(Vector3(0.0, 0.0, 1.0));
            }
        }
        else {
            // loop through pairs of atoms to find sigma v plane
            // candidates
            for (i=0; i<natom(); ++i) {
                Vector3 A = xyz(i) - com;
                double AdotA = A.dot(A);
                // the second atom can equal i because i might be
                // in the plane
                for (j=0; j<=i; ++j) {
                    // the atoms must be identical
                    if (Z(i) != Z(j) || fabs(mass(i) - mass(j)) > tol) continue;
                    Vector3 B = xyz(j) - com;
                    // the atoms must be the same distance from the com
                    if (fabs(AdotA - B.dot(B)) > tol) continue;
                    Vector3 inplane = B+A;
                    double norm_inplane = inplane.norm();
                    if (norm_inplane < tol) continue;
                    inplane *= 1.0/norm_inplane;
                    Vector3 perp = c2axis.cross(inplane);
                    double norm_perp = perp.norm();
                    if (norm_perp < tol) continue;
                    perp *= 1.0/norm_perp;
                    if (is_plane(com, perp, tol)) {
                        have_sigmav = true;
                        sigmav = perp;
                        goto found_sigmav;
                    }
                }
            }
        }
    }

found_sigmav:
    if (have_sigmav) {
        // try to make the sign of the oop vec correspond to one of
        // the world axes
        int sigmavlike = like_world_axis(sigmav, worldxaxis, worldyaxis, worldzaxis);

        // Choose sigmav to be the world x axis, if possible
        if (c2like == ZAxis && sigmavlike == YAxis) {
            sigmav = sigmav.cross(c2axis);
        }
        else if (c2like == YAxis && sigmavlike == ZAxis) {
            sigmav = c2axis.cross(sigmav);
        }
    }

    // under certain conditions i need to know if there is any sigma
    // plane
    bool have_sigma = false;
    Vector3 sigma;
    if (!have_inversion && !have_c2axis) {
        if (planar) {
            // find two noncolinear atom-atom vectors
            // we know that linear==0 since !have_c2axis
            Vector3 BA = xyz(1) - xyz(0);
            BA.normalize();
            for (i=2; i<natom(); ++i) {
                Vector3 CA = xyz(i) - xyz(0);
                CA.normalize();
                Vector3 BAxCA = BA.cross(CA);
                if (BAxCA.norm() > tol) {
                    have_sigma = true;
                    BAxCA.normalize();
                    sigma = BAxCA;
                    break;
                }
            }
        }
        else {
            // loop through pairs of atoms to contruct trial planes
            for (i=0; i<natom(); ++i) {
                Vector3 A = xyz(i) - com;
                double AdotA = A.dot(A);
                for (j=0; j<i; ++j) {
                    // the atomsmust be identical
                    if (Z(i) != Z(j) || fabs(mass(i)-mass(j)) > tol) continue;
                    Vector3 B = xyz(j)-com;
                    double BdotB = B.dot(B);
                    // the atoms must be the same distance from the com
                    if (fabs(AdotA - BdotB) > tol) continue;
                    Vector3 perp = B-A;
                    double norm_perp = perp.norm();
                    if (norm_perp < tol) continue;
                    perp *= 1.0 / norm_perp;
                    if (is_plane(com, perp, tol)) {
                        have_sigma = true;
                        sigma = perp;
                        goto found_sigma;
                    }
                }
            }
        }
    }
found_sigma:

    if (have_sigma) {
        // try to make the sign of the oop vec correspond to one of
        // the world axes
        double xlikeness = fabs(sigma.dot(worldxaxis));
        double ylikeness = fabs(sigma.dot(worldyaxis));
        double zlikeness = fabs(sigma.dot(worldzaxis));

        if (xlikeness > ylikeness && xlikeness > zlikeness) {
            if (sigma.dot(worldxaxis) < 0) sigma = -sigma;
        }
        else if (ylikeness > zlikeness) {
            if (sigma.dot(worldyaxis) < 0) sigma = -sigma;
        }
        else {
            if (sigma.dot(worldzaxis) < 0) sigma = -sigma;
        }
    }

//    fprintf(outfile, "find point group:\n");
//    fprintf(outfile, "  linear          = %s\n", linear          ? "true" : "false");
//    fprintf(outfile, "  planar          = %s\n", planar          ? "true" : "false");
//    fprintf(outfile, "  have_inversion  = %s\n", have_inversion  ? "true" : "false");
//    fprintf(outfile, "  have_c2axis     = %s\n", have_c2axis     ? "true" : "false");
//    fprintf(outfile, "  have_c2axisperp = %s\n", have_c2axisperp ? "true" : "false");
//    fprintf(outfile, "  have_sigmav     = %s\n", have_sigmav     ? "true" : "false");
//    fprintf(outfile, "  have_sigma      = %s\n", have_sigma      ? "true" : "false");
//
//    if (have_c2axis)
//        fprintf(outfile, "  c2axis          = %s\n", c2axis.to_string().c_str());
//    if (have_c2axisperp)
//        fprintf(outfile, "  c2axisperp      = %s\n", c2axisperp.to_string().c_str());
//    if (have_sigmav)
//        fprintf(outfile, "  sigmav          = %s\n", sigmav.to_string().c_str());
//    if (have_sigma)
//        fprintf(outfile, "  sigma           = %s\n", sigma.to_string().c_str());

    // Find the three axes for the symmetry frame
    Vector3 xaxis = worldxaxis;
    Vector3 yaxis;
    Vector3 zaxis = worldzaxis;
    if (have_c2axis) {
        zaxis = c2axis;
        if (have_sigmav) {
            xaxis = sigmav;
        }
        else if (have_c2axisperp) {
            xaxis = c2axisperp;
        }
        else {
            // any axis orthogonal to the zaxis will do
            xaxis = zaxis.perp_unit(zaxis);
        }
    }
    else if (have_sigma) {
        zaxis = sigma;
        xaxis = zaxis.perp_unit(zaxis);
    }
    // the y is then -x cross z
    yaxis = -xaxis.cross(zaxis);

//    fprintf(outfile, "  X: %s\n", xaxis.to_string().c_str());
//    fprintf(outfile, "  Y: %s\n", yaxis.to_string().c_str());
//    fprintf(outfile, "  Z: %s\n", zaxis.to_string().c_str());

    SymmetryOperation frame;
    Vector3 origin;
    for (i=0; i<3; ++i) {
        frame(i,0) = xaxis[i];
        frame(i,1) = yaxis[i];
        frame(i,2) = zaxis[i];
        origin[i] = com[i];
    }

    boost::shared_ptr<PointGroup> pg;
    if (have_inversion) {
        if (have_c2axis) {
            if (have_sigmav) {
                pg = shared_ptr<PointGroup>(new PointGroup("d2h", frame, origin));
            }
            else {
                pg = shared_ptr<PointGroup>(new PointGroup("c2h", frame, origin));
            }
        }
        else {
            pg = shared_ptr<PointGroup>(new PointGroup("ci", frame, origin));
        }
    }
    else {
        if (have_c2axis) {
            if (have_sigmav) {
                pg = shared_ptr<PointGroup>(new PointGroup("c2v", frame, origin));
            }
            else {
                if (have_c2axisperp) {
                    pg = shared_ptr<PointGroup>(new PointGroup("d2", frame, origin));
                }
                else {
                    pg = shared_ptr<PointGroup>(new PointGroup("c2", frame, origin));
                }
            }
        }
        else {
            if (have_sigma) {
                pg = shared_ptr<PointGroup>(new PointGroup("cs", frame, origin));
            }
            else {
                pg = shared_ptr<PointGroup>(new PointGroup("c1", frame, origin));
            }
        }
    }

    //fprintf(outfile, "\n  Molecular point group: %s\n\n", pg->symbol());

    return pg;
}

void Molecule::release_symmetry_information()
{
    for (int i=0; i<nunique_; ++i) {
    delete[] equiv_[i];
    }
    delete[] equiv_;
    delete[] nequiv_;
    delete[] atom_to_unique_;
    nunique_ = 0;
    equiv_   = 0;
    nequiv_  = 0;
    atom_to_unique_ = 0;
}

void Molecule::form_symmetry_information(double tol)
{
    if (equiv_)
    release_symmetry_information();

    if (natom() == 0) {
        nunique_ = 0;
        equiv_   = 0;
        nequiv_  = 0;
        atom_to_unique_ = 0;
        return;
    }

    nequiv_         = new int[natom()];
    atom_to_unique_ = new int[natom()];
    equiv_          = new int*[natom()];

    if (!strcmp(point_group()->symbol(), "c1")) {
        nunique_ = natom();
        for (int i=0; i<natom(); ++i) {
            nequiv_[i] = 1;
            equiv_[i] = new int[1];
            equiv_[i][0] = i;
            atom_to_unique_[i] = i;
        }
        return;
    }

    // The first atom is always unique
    nunique_           = 1;
    nequiv_[0]         = 1;
    equiv_[0]          = new int[1];
    equiv_[0][0]       = 0;
    atom_to_unique_[0] = 0;

    CharacterTable ct  = point_group()->char_table();

    Vector3 ac;
    SymmetryOperation so;
    Vector3 np;

    // Find the equivalent atoms
    int i;
    for (i=1; i<natom(); ++i) {
        ac = xyz(i);
        int i_is_unique = 1;
        int i_equiv = 0;

        // Apply all symmetry ops in the group to the atom
        for (int g=0; g<ct.order(); ++g) {
            so = ct.symm_operation(g);
            for (int ii=0; ii<3; ++ii) {
                np[ii] = 0;
                for (int jj=0; jj<3; ++jj)
                    np[ii] += so(ii, jj) * ac[jj];
            }

            // See if the transformed atom is equivalent to a
            // unique atom
            for (int j=0; j<nunique_; ++j) {
                int unique = equiv_[j][0];
                Vector3 aj(xyz(unique));
                if (np.distance(aj) < tol
                        && Z(unique) == Z(i)
                        && fabs(mass(unique)-mass(i)) < tol) {
                    i_is_unique = 0;
                    i_equiv = j;
                    break;
                }
            }
        }
        if (i_is_unique) {
            nequiv_[nunique_] = 1;
            equiv_[nunique_] = new int[1];
            equiv_[nunique_][0] = i;
            atom_to_unique_[i] = nunique_;
            nunique_++;
        }
        else {
            int *tmp = new int[nequiv_[i_equiv]+1];
            memcpy(tmp, equiv_[i_equiv],nequiv_[i_equiv]*sizeof(int));
            delete[] equiv_[i_equiv];
            equiv_[i_equiv] = tmp;
            equiv_[i_equiv][nequiv_[i_equiv]] = i;
            nequiv_[i_equiv]++;
            atom_to_unique_[i] = i_equiv;
        }
    }

    // The first atom in the equiv list is considered the primary
    // unique atom. Just to make things look pretty, amek the
    // atom with the most zeros in its x, y, z coordinate the
    // unique atom. Nothing else should rely on this being done.
    double ztol=1.0e-5;
    for (i=0; i<nunique_; ++i) {
        int maxzero = 0;
        int jmaxzero = 0;
        for (int j=0; j<nequiv_[i]; ++j) {
            int nzero = 0;
            for (int k=0; k<3; ++k)
                if (fabs(xyz(equiv_[i][j], k)) < ztol)
                    nzero++;
            if (nzero > maxzero) {
                maxzero = nzero;
                jmaxzero = j;
            }
        }
        int tmp = equiv_[i][jmaxzero];
        equiv_[i][jmaxzero] = equiv_[i][0];
        equiv_[i][0] = tmp;
    }
}

const char* Molecule::sym_label()
{
  if (pg_==NULL) set_point_group(find_point_group());
  const char *symlabel;
  symlabel = pg_->symbol();
  return symlabel;
}

char** Molecule::irrep_labels()
{
  if (pg_==NULL) set_point_group(find_point_group());
  int nirreps = pg_->char_table().nirrep();
  char **irreplabel = (char **) malloc(sizeof(char *)*nirreps);
  for (int i=0; i<nirreps; i++) {
    irreplabel[i] = (char *) malloc(sizeof(char)*5);
    strcpy(irreplabel[i],pg_->char_table().gamma(i).symbol());
  }
  return irreplabel;
}

Vector3 Molecule::xyz(int atom) const
{
    return Vector3(atoms_[atom].x, atoms_[atom].y, atoms_[atom].z);
}

Vector3 Molecule::fxyz(int atom) const
{
    return Vector3(full_atoms_[atom].x, full_atoms_[atom].y, full_atoms_[atom].z);
}

//////////////////////////////////////////////////////////////////////////////
//
// Atomic Radii
//
//////////////////////////////////////////////////////////////////////////////

AtomicRadii::AtomicRadii()
    : count_(0), radii_(0)
{

}

AtomicRadii::~AtomicRadii()
{
    if (radii_)
        delete[] radii_;
    radii_ = 0;
    count_ = 0;
}

TreutlerAtomicRadii::TreutlerAtomicRadii()
    : AtomicRadii()
{
    count_ = 37;
    radii_ = new double[count_];

    radii_[0]  = 1.0; /* Ghost Atom */
    radii_[1]  = 0.8; /* H */
    radii_[2]  = 0.9; /* He */
    radii_[3]  = 1.8; /* Li */
    radii_[4]  = 1.4; /* Be */
    radii_[5]  = 1.3; /* B */
    radii_[6]  = 1.1; /* C */
    radii_[7]  = 0.9; /* N */
    radii_[8]  = 0.9; /* O */
    radii_[9]  = 0.9; /* F */
    radii_[10] = 0.9; /* Ne */
    radii_[11] = 1.4; /* Na */
    radii_[12] = 1.3; /* Mg */
    radii_[13] = 1.3; /* Al */
    radii_[14] = 1.2; /* Si */
    radii_[15] = 1.1; /* P */
    radii_[16] = 1.0; /* S */
    radii_[17] = 1.0; /* Cl */
    radii_[18] = 1.0; /* Ar */
    radii_[19] = 1.5; /* K */
    radii_[20] = 1.4; /* Ca */
    radii_[21] = 1.3; /* Sc */
    radii_[22] = 1.2; /* Ti */
    radii_[23] = 1.2; /* V */
    radii_[24] = 1.2; /* Cr */
    radii_[25] = 1.2; /* Mn */
    radii_[26] = 1.2; /* Fe */
    radii_[27] = 1.2; /* Co */
    radii_[28] = 1.1; /* Ni */
    radii_[29] = 1.1; /* Cu */
    radii_[30] = 1.1; /* Zn */
    radii_[31] = 1.1; /* Ga */
    radii_[32] = 1.0; /* Ge */
    radii_[33] = 0.9; /* As */
    radii_[34] = 0.9; /* Se */
    radii_[35] = 0.9; /* Br */
    radii_[36] = 0.9; /* Kr */
}

AdjustedVanDerWaalsAtomicRadii::AdjustedVanDerWaalsAtomicRadii()
{
    count_ = 55;
    radii_ = new double[count_];

    radii_[0]  = 1.000; /* Default element or ghost */
    radii_[1]  = 1.001; /* H  (untouched) */
    radii_[2]  = 0.5*1.012; /* He */
    radii_[3]  = 0.5*0.825; /* Li */
    radii_[4]  = 0.5*1.408; /* Be */
    radii_[5]  = 0.5*1.485; /* B  */
    radii_[6]  = 0.5*1.452; /* C  */
    radii_[7]  = 0.5*1.397; /* N  */
    radii_[8]  = 0.5*1.342; /* O  */
    radii_[9]  = 0.5*1.287; /* F  */
    radii_[10] = 0.5*1.243; /* Ne */
    radii_[11] = 0.5*1.144; /* Na */
    radii_[12] = 0.5*1.364; /* Mg */
    radii_[13] = 0.5*1.639; /* Al */
    radii_[14] = 0.5*1.716; /* Si */
    radii_[15] = 0.5*1.705; /* P  */
    radii_[16] = 0.5*1.683; /* S  */
    radii_[17] = 0.5*1.639; /* Cl */
    radii_[18] = 0.5*1.595; /* Ar */
    radii_[19] = 0.5*1.485; /* K  */
    radii_[20] = 0.5*1.474; /* Ca */
    radii_[21] = 0.5*1.562; /* Sc */
    radii_[22] = 0.5*1.562; /* Ti */
    radii_[23] = 0.5*1.562; /* V  */
    radii_[24] = 0.5*1.562; /* Cr */
    radii_[25] = 0.5*1.562; /* Mn */
    radii_[26] = 0.5*1.562; /* Fe */
    radii_[27] = 0.5*1.562; /* Co */
    radii_[28] = 0.5*1.562; /* Ni */
    radii_[29] = 0.5*1.562; /* Cu */
    radii_[30] = 0.5*1.562; /* Zn */
    radii_[31] = 0.5*1.650; /* Ga */
    radii_[32] = 0.5*1.727; /* Ge */
    radii_[33] = 0.5*1.760; /* As */
    radii_[34] = 0.5*1.771; /* Se */
    radii_[35] = 0.5*1.749; /* Br */
    radii_[36] = 0.5*1.727; /* Kr */
    radii_[37] = 0.5*1.628; /* Rb */
    radii_[38] = 0.5*1.606; /* Sr */
    radii_[39] = 0.5*1.639; /* Y  */
    radii_[40] = 0.5*1.639; /* Zr */
    radii_[41] = 0.5*1.639; /* Nb */
    radii_[42] = 0.5*1.639; /* Mo */
    radii_[43] = 0.5*1.639; /* Tc */
    radii_[44] = 0.5*1.639; /* Ru */
    radii_[45] = 0.5*1.639; /* Rh */
    radii_[46] = 0.5*1.639; /* Pd */
    radii_[47] = 0.5*1.639; /* Ag */
    radii_[48] = 0.5*1.639; /* Cd */
    radii_[49] = 0.5*1.672; /* In */
    radii_[50] = 0.5*1.804; /* Sn */
    radii_[51] = 0.5*1.881; /* Sb */
    radii_[52] = 0.5*1.892; /* Te */
    radii_[53] = 0.5*1.892; /* I  */
    radii_[54] = 0.5*1.881; /* Xe */
}

BraggSlaterAtomicRadii::BraggSlaterAtomicRadii()
    : AtomicRadii()
{
    count_ = 0;
    radii_ = new double[count_];

}

SG1AtomicRadii::SG1AtomicRadii()
    : AtomicRadii()
{
    count_ = 19;
    radii_ = new double[count_];

    radii_[0]  = 2.0000; /* Ghost */
    radii_[1]  = 1.0000; /* H */
    radii_[2]  = 0.5882; /* He */
    radii_[3]  = 3.0769; /* Li */
    radii_[4]  = 2.0513; /* Be */
    radii_[5]  = 1.5385; /* B */
    radii_[6]  = 1.2308; /* C */
    radii_[7]  = 1.0256; /* N */
    radii_[8]  = 0.8791; /* O */
    radii_[9]  = 0.7692; /* F */
    radii_[10] = 0.6838; /* Ne */
    radii_[11] = 4.0909; /* Na */
    radii_[12] = 3.1579; /* Mg */
    radii_[13] = 2.5714; /* Al */
    radii_[14] = 2.1687; /* Si */
    radii_[15] = 1.8750; /* P */
    radii_[16] = 1.6514; /* S */
    radii_[17] = 1.4754; /* Cl */
    radii_[18] = 1.3333; /* Ar */
}

