#include <libciomr/libciomr.h>

#include "mints.h"
#include "extern.h"

using namespace boost;
using namespace psi;

ExternalPotential::ExternalPotential()
{
}
ExternalPotential::~ExternalPotential()
{
}
void ExternalPotential::clear()
{
    charges_.clear();
    dipoles_.clear();
    quadrupoles_.clear();
}
void ExternalPotential::addCharge(double Z, double x, double y, double z)
{
    charges_.push_back(make_tuple(Z,x,y,z));
}
void ExternalPotential::addDipole(double mx, double my, double mz, double x, double y, double z)
{
    dipoles_.push_back(make_tuple(mx,my,mz,x,y,z));
}
void ExternalPotential::addQuadrupole(double qxx, double qxy, double qxz, double qyy, double qyz, double qzz, double x, double y, double z)
{
    quadrupoles_.push_back(make_tuple(qxx,qxy,qxz,qyy,qyz,qzz,x,y,z));
}
void ExternalPotential::print(FILE* out) const
{
    fprintf(out, "  ==> External Potential Field: %s <==\n\n", name_.c_str());

    // Charges
    fprintf(out, "   => Charges [a.u.] <=\n\n");
    fprintf(out, "    %8s %8s %8s %8s\n","Z","x","y","z");
    fprintf(out, "    ---------------------------------------------\n");
    for (int i = 0 ; i < charges_.size(); i++) {
        fprintf(out, "   %8.5f %8.5f %8.5f %8.5f\n",
            get<0>(charges_[i]), get<1>(charges_[i]), get<2>(charges_[i]), get<3>(charges_[i]));
    }
    fprintf(out,"\n");

    // Dipoles
    //TODO

    // Quadrupoles
    // TODO

    fflush(out);
}
void ExternalPotential::translate(double dx, double dy, double dz)
{
    for (int i = 0; i < charges_.size(); i++) {
        get<1> (charges_[i]) += dx;
        get<2> (charges_[i]) += dy;
        get<3> (charges_[i]) += dz;
    }
    for (int i = 0; i < dipoles_.size(); i++) {
        get<3> (dipoles_[i]) += dx;
        get<4> (dipoles_[i]) += dy;
        get<5> (dipoles_[i]) += dz;
    }
    for (int i = 0; i < quadrupoles_.size(); i++) {
        get<6> (quadrupoles_[i]) += dx;
        get<7> (quadrupoles_[i]) += dy;
        get<8> (quadrupoles_[i]) += dz;
    }
}
void ExternalPotential::rotate(shared_ptr<Matrix> R)
{
    //TODO
    throw FeatureNotImplemented("psi::ExternalPotential::rotate", "Rotation of multipoles not yet implemented", __FILE__,__LINE__);
}
double ExternalPotential::computePotentialPoint(double x, double y, double z)
{
    double V = 0.0;

    // Charges
    for (int i = 0; i < charges_.size(); i++) {
        double Z = get<0>(charges_[i]);

        double dx = x - get<1>(charges_[i]);
        double dy = y - get<2>(charges_[i]);
        double dz = z - get<3>(charges_[i]);

        double R  = sqrt(dx*dx + dy*dy + dz*dz);

        V -= Z / R;
    }

    // Dipoles
    for (int i = 0; i < dipoles_.size(); i++) {

        double mx = get<0>(dipoles_[i]);
        double my = get<1>(dipoles_[i]);
        double mz = get<2>(dipoles_[i]);

        double dx = x - get<3>(dipoles_[i]);
        double dy = y - get<4>(dipoles_[i]);
        double dz = z - get<5>(dipoles_[i]);

        double R  = sqrt(dx*dx + dy*dy + dz*dz);

        V -= (mx*dx + my*dy + mz*dz) / (R*R*R);
    }

    // Quadrupoles
    for (int i = 0; i < dipoles_.size(); i++) {

        double qxx = get<0>(quadrupoles_[i]);
        double qxy = get<1>(quadrupoles_[i]);
        double qxz = get<2>(quadrupoles_[i]);
        double qyy = get<3>(quadrupoles_[i]);
        double qyz = get<4>(quadrupoles_[i]);
        double qzz = get<5>(quadrupoles_[i]);

        double dx = x - get<6>(quadrupoles_[i]);
        double dy = y - get<7>(quadrupoles_[i]);
        double dz = z - get<8>(quadrupoles_[i]);

        double R  = sqrt(dx*dx + dy*dy + dz*dz);

        V -= (0.5 * qxx * dx * dx
            + 0.5 * qyy * dy * dy
            + 0.5 * qzz * dz * dz
                  + qxy * dx * dy
                  + qxz * dx * dz
                  + qyz * dy * dz)
            / (R*R*R*R*R);
    }

    return V;
}
shared_ptr<Matrix> ExternalPotential::computePotentialMatrix(shared_ptr<BasisSet> basis)
{
    int n = basis->nbf();
    shared_ptr<Matrix> V(new Matrix("External Potential",n,n));
    shared_ptr<IntegralFactory> fact(new IntegralFactory(basis,basis,basis,basis));

    // Monopoles
    shared_ptr<Matrix> V_charge(new Matrix("External Potential (Charges)", n, n));

    shared_ptr<Matrix> Zxyz(new Matrix("Charges (Z,x,y,z)", charges_.size(), 4));
    double** Zxyzp = Zxyz->pointer();
    for (int i = 0; i < charges_.size(); i++) {
        Zxyzp[i][0] = get<0>(charges_[i]);
        Zxyzp[i][1] = get<1>(charges_[i]);
        Zxyzp[i][2] = get<2>(charges_[i]);
        Zxyzp[i][3] = get<3>(charges_[i]);
    }

    shared_ptr<PotentialInt> pot(static_cast<PotentialInt*>(fact->ao_potential()));
    pot->set_charge_field(Zxyz);
    pot->compute(V_charge);

    V->add(V_charge);
    V_charge.reset();
    pot.reset();

    // Dipoles
    // TODO
    if (dipoles_.size())
        throw FeatureNotImplemented("psi::ExternalPotential::computePotentialMatrix","Dipole integrals not implemented",__FILE__,__LINE__);

    // Quadrupoles
    // TODO
    if (quadrupoles_.size())
        throw FeatureNotImplemented("psi::ExternalPotential::computePotentialMatrix","Quadrupole integrals not implemented",__FILE__,__LINE__);

    return V;
}


