#include <cmath>
#include <cstdio>

#include <libmints/molecule.h>
#include <libmints/matrix.h>
#include <libmints/pointgrp.h>
#include <libciomr/libciomr.h>

#include <masses.h>
#include <physconst.h>

using namespace std;
using namespace psi;
using namespace boost;

// TODO: These should probably be moved to psi4-def.h
#define ZERO_MOMENT_INERTIA 1.0E-10     /*Tolerance for degenerate rotational constants*/
#define ZERO 1.0E-14

namespace psi {
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
    nirreps_(0)
{

}

Molecule::Molecule(const Molecule& other)
{
    atoms_ = other.atoms_;
    nirreps_ = other.nirreps_;
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
                        const char *label, double mass,
                        int have_charge, double charge)
{
    atom_info info;

    info.x = x;
    info.y = y;
    info.z = z;
    info.Z = Z;
    info.charge = charge;
    if (label != NULL)
        info.label = label;
    else
        info.label = "Gh";
    info.mass = mass;

    if (Z > 0)
        atoms_.push_back(info);
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

SimpleVector Molecule::nuclear_repulsion_energy_deriv1()
{
    SimpleVector de(3*natom());

    for (int i=1; i<natom(); ++i) {
        for (int j=0; j<natom(); ++j) {
            if (i != j) {
                double temp = pow((xyz(i).distance(xyz(j))), 3.0);
                de[3*i+0] -= (x(i) - x(j)) * Z(i) * Z(j) / temp;
                de[3*i+1] -= (y(i) - y(j)) * Z(i) * Z(j) / temp;
                de[3*i+2] -= (z(i) - z(j)) * Z(i) * Z(j) / temp;
            }
        }
    }

    return de;
}

/*
    TODO Test nuclear_repulsion_energy_deriv2
*/
SimpleMatrix* Molecule::nuclear_repulsion_energy_deriv2()
{
	SimpleMatrix *hess = new SimpleMatrix("Nuclear Repulsion Energy 2nd Derivatives", 3*natom(), 3*natom());
	double sx, sy, sz, x2, y2, z2, r2, r, r5, pfac;

	for (int i=1; i<natom(); ++i) {
		int ix = 3*i;
		int iy = ix+1;
		int iz = iy+1;

		for (int j=0; j<i; ++j) {
			int jx = 3*j;
			int jy = jx+j;
			int jz = jy+j;

			sx = x(i) - x(j);
			sy = y(i) - y(j);
			sz = z(i) - z(j);

			x2 = sx*sx; y2 = sy*sy; z2 = sz*sz;
			r2 = x2 + y2 + z2;
			r = sqrt(r2);
			r5 = r2*r2*r;
			pfac = Z(i) * Z(j) / r5;

			hess->add(ix, ix, pfac * (3*x2 - r2));
			hess->add(iy, iy, pfac * (3*y2 - r2));
			hess->add(iz, iz, pfac * (3*z2 - r2));
			hess->add(ix, iy, pfac*3*sx*sy);
			hess->add(ix, iz, pfac*3*sx*sz);
			hess->add(iy, iz, pfac*3*sy*sz);

			hess->add(jx, jx, pfac * (3*x2 - r2));
			hess->add(jy, jy, pfac * (3*y2 - r2));
			hess->add(jz, jz, pfac * (3*z2 - r2));
			hess->add(jx, jy, pfac*3*sx*sy);
			hess->add(jx, jz, pfac*3*sx*sz);
			hess->add(jy, jz, pfac*3*sy*sz);

			hess->add(ix, jx, -pfac*(3*sx*sx-r2));
			hess->add(ix, jy, -pfac*(3*sx*sy));
			hess->add(ix, jz, -pfac*(3*sx*sz));
			hess->add(iy, jx, -pfac*(3*sy*sx));
			hess->add(iy, jy, -pfac*(3*sy*sy-r2));
			hess->add(iy, jz, -pfac*3*sy*sz);
			hess->add(iz, jx, -pfac*3*sz*sx);
			hess->add(iz, jy, -pfac*3*sz*sy);
			hess->add(iz, jz, -pfac*(3*sz*sz-r2));
		}
	}
	return hess;
}

void Molecule::translate(const Vector3& r)
{
    for (int i=0; i<natom(); ++i) {
        atoms_[i].x += r[0];
        atoms_[i].y += r[1];
        atoms_[i].z += r[2];
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

void Molecule::set_geometry(SimpleMatrix& geom)
{
    for (int i=0; i<natom(); ++i) {
        atoms_[i].x = geom[i][0];
        atoms_[i].y = geom[i][1];
        atoms_[i].z = geom[i][2];
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

void Molecule::reorient()
{
    // If no atoms are present throw
    if (natom() <= 0) {
        throw PSIEXCEPTION("Molecule::reorient called with no atoms.");
    }

    // Nothing for us to do.
    if (natom() == 1)
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
        rotate(R);
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

            R.zero();
            R[0][0] = 1.0;
            R[1][1] = cos_theta;
            R[2][2] = cos_theta;
            R[1][2] = sin_theta;
            R[2][1] = -sin_theta;
            rotate(R);

            R.zero();
            R[2][2] = 1.0;
            R[0][0] = cos_phix;
            R[1][1] = cos_phix;
            R[0][1] = -sin_phix;
            R[1][0] = sin_phix;
            rotate(R);
        }
    }

    // Delete the tensor matrix
    delete itensor;
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
        add_atom((int)zvals[i], geom[i][0], geom[i][1], geom[i][2], atomic_labels[(int)zvals[i]], an2masses[(int)zvals[i]]);
    }

    nirreps_ = chkpt->rd_nirreps();

    Chkpt::free(zvals);
    Chkpt::free(geom);
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
        fprintf(outfile,"       Center              X                  Y                   Z\n");
        fprintf(outfile,"    ------------   -----------------  -----------------  -----------------\n");

        for(int i = 0; i < natom(); ++i){
            Vector3 geom = xyz(i);
            fprintf(outfile, "    %12s ",label(i).c_str()); fflush(outfile);
            for(int j = 0; j < 3; j++)
                fprintf(outfile, "  %17.12f", geom[j]);
            fprintf(outfile,"\n");
        }
        fprintf(outfile,"\n");
        fflush(outfile);
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

    // Print the molecule we're working with
    fprintf(outfile, "natom() = %d\n", natom());
    print();

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

    fprintf(outfile, "find point group:\n");
    fprintf(outfile, "  linear          = %s\n", linear          ? "true" : "false");
    fprintf(outfile, "  planar          = %s\n", planar          ? "true" : "false");
    fprintf(outfile, "  have_inversion  = %s\n", have_inversion  ? "true" : "false");
    fprintf(outfile, "  have_c2axis     = %s\n", have_c2axis     ? "true" : "false");
    fprintf(outfile, "  have_c2axisperp = %s\n", have_c2axisperp ? "true" : "false");
    fprintf(outfile, "  have_sigmav     = %s\n", have_sigmav     ? "true" : "false");
    fprintf(outfile, "  have_sigma      = %s\n", have_sigma      ? "true" : "false");

    if (have_c2axis)
        fprintf(outfile, "  c2axis          = %s\n", c2axis.to_string().c_str());
    if (have_c2axisperp)
        fprintf(outfile, "  c2axisperp      = %s\n", c2axisperp.to_string().c_str());
    if (have_sigmav)
        fprintf(outfile, "  sigmav          = %s\n", sigmav.to_string().c_str());
    if (have_sigma)
        fprintf(outfile, "  sigma           = %s\n", sigma.to_string().c_str());

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

    fprintf(outfile, "X: %s\n", xaxis.to_string().c_str());
    fprintf(outfile, "Y: %s\n", yaxis.to_string().c_str());
    fprintf(outfile, "Z: %s\n", zaxis.to_string().c_str());

    SymmetryOperation frame;
    Vector3 origin;
    for (i=0; i<3; ++i) {
        frame(i,0) = xaxis[i];
        frame(i,1) = yaxis[i];
        frame(i,2) = zaxis[i];
        origin[i] = com[i];
    }

    fprintf(outfile, "frame:\n");
    frame.print(outfile);
    fprintf(outfile, "origin: %s\n", origin.to_string().c_str());

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

    if (!strcmp(point_group))
}