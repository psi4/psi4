#include <cmath>
#include <cstdio>

#include <libmints/molecule.h>
#include <libmints/matrix.h>

#include <masses.h>
#include <physconst.h>

using namespace std;
using namespace psi;

Molecule::Molecule():
    natoms_(0), nirreps_(0)
{

}

Molecule::~Molecule()
{
    clear();
}

void Molecule::clear()
{
    natoms_ = 0;
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

    natoms_++;
    atoms_.push_back(info);
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

int Molecule::atom_at_position(double *xyz, double tol) const
{
    Vector3 b(xyz);
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

void Molecule::init_with_chkpt(shared_ptr<PSIO> psio)
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

void Molecule::print(FILE *out)
{
    if (natom()) {
        fprintf(out,"       Center              X                  Y                   Z\n");
        fprintf(out,"    ------------   -----------------  -----------------  -----------------\n");

        for(int i = 0; i < natom(); ++i){
            Vector3 geom = xyz(i);
            fprintf(out, "    %12s ",label(i).c_str()); fflush(out);
            for(int j = 0; j < 3; j++)
                fprintf(out, "  %17.12f", geom[j]);
            fprintf(out,"\n");
        }
        fprintf(out,"\n");
        fflush(out);
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
