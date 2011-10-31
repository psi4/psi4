#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <utility>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.hpp>
#include <libciomr/libciomr.h>
#include "mints.h"
#include <libqt/qt.h>
#include <psi4-dec.h>
#include <physconst.h>

using namespace boost;
using namespace psi;
using namespace std;

namespace psi {

Prop::Prop(boost::shared_ptr<Wavefunction> wfn) : wfn_(wfn)
{
    if (wfn_.get() == NULL)
        throw PSIEXCEPTION("Prop: Wavefunction is null");
    common_init();
}
Prop::~Prop()
{
}
void Prop::common_init()
{
    tasks_.clear();

    basisset_ = wfn_->basisset();
    restricted_ = wfn_->restricted();

    integral_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));

    boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral_));
    AO2USO_ = pet->aotoso();
    factory_ = wfn_->matrix_factory();

    // For now wavefunction has SO quantities only, so we'll use those
    Ca_so_ = wfn_->Ca();
    Da_so_ = wfn_->Da();

    if (restricted_) {
        Cb_so_ = Ca_so_;
        Db_so_ = Da_so_;
    } else {
        Cb_so_ = wfn_->Cb();
        Db_so_ = wfn_->Db();
    }
}
void Prop::add(const std::string& prop)
{
    tasks_.insert(prop);
}
void Prop::add(std::vector<std::string> props)
{
    for (int i = 0; i < props.size(); i++) {
        tasks_.insert(props[i]);
    }
}
SharedMatrix Prop::Da_ao()
{
    int nao = basisset_->nbf();
    SharedMatrix D = SharedMatrix(new Matrix("Da (AO basis)", basisset_->nbf(), basisset_->nbf()));
    double** Dp = D->pointer();
    for (int h = 0; h < AO2USO_->nirrep(); h++) {
        int nso = AO2USO_->colspi()[h];
        if (nso == 0) continue;
        double** Da = Da_so_->pointer(h);
        double** X = AO2USO_->pointer(h);
        double** Temp = block_matrix(nao,nso);
        C_DGEMM('N','N',nao,nso,nso,1.0,X[0],nso,Da[0],nso,0.0,Temp[0],nso);
        C_DGEMM('N','T',nao,nao,nso,1.0,Temp[0],nso,X[0],nso,1.0,Dp[0],nao);
        free_block(Temp);
    }
    return D;
}
SharedMatrix Prop::Db_ao()
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, asking for Db makes no sense");

    int nao = basisset_->nbf();
    SharedMatrix D = SharedMatrix(new Matrix("Db (AO basis)", basisset_->nbf(), basisset_->nbf()));
    double** Dp = D->pointer();
    for (int h = 0; h < AO2USO_->nirrep(); h++) {
        int nso = AO2USO_->colspi()[h];
        if (nso == 0) continue;
        double** Db = Db_so_->pointer(h);
        double** X = AO2USO_->pointer(h);
        double** Temp = block_matrix(nao,nso);
        C_DGEMM('N','N',nao,nso,nso,1.0,X[0],nso,Db[0],nso,0.0,Temp[0],nso);
        C_DGEMM('N','T',nao,nao,nso,1.0,Temp[0],nso,X[0],nso,1.0,Dp[0],nao);
        free_block(Temp);
    }
    return D;
}
SharedMatrix Prop::Ca_ao()
{
    // TODO reorder by eigenvalue
    int nao = basisset_->nbf();
    SharedMatrix C = SharedMatrix(new Matrix("Ca (AO basis)", basisset_->nbf(), basisset_->nbf()));
    double** Cp = C->pointer();
    int counter = 0;
    for (int h = 0; h < AO2USO_->nirrep(); h++) {
        int nso = AO2USO_->colspi()[h];
        int nmo = C->colspi()[h];
        if (nso == 0 || nmo == 0) continue;
        double** Ca = Ca_so_->pointer(h);
        double** X = AO2USO_->pointer(h);

        C_DGEMM('N','N',nao,nmo,nso,1.0,X[0],nso,Ca[0],nmo,0.0,&Cp[0][counter],nao);

        counter += nmo;
    }
    return C;
}
SharedMatrix Prop::Cb_ao()
{
    // TODO reorder by eigenvalue
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, asking for Cb makes no sense");

    int nao = basisset_->nbf();
    SharedMatrix C = SharedMatrix(new Matrix("Cb (AO basis)", basisset_->nbf(), basisset_->nbf()));
    double** Cp = C->pointer();
    int counter = 0;
    for (int h = 0; h < AO2USO_->nirrep(); h++) {
        int nso = AO2USO_->colspi()[h];
        int nmo = C->colspi()[h];
        if (nso == 0 || nmo == 0) continue;
        double** Ca = Ca_so_->pointer(h);
        double** X = AO2USO_->pointer(h);

        C_DGEMM('N','N',nao,nmo,nso,1.0,X[0],nso,Ca[0],nmo,0.0,&Cp[0][counter],nao);

        counter += nmo;
    }
    return C;
}
SharedMatrix Prop::Da_mo()
{
    // MO D are nso x nso, zeros padding if nmo < nso
    SharedMatrix D = factory_->create_shared_matrix("Da (MO Basis)");
    for (int h = 0; h < D->nirrep(); h++) {
        int nso = D->colspi()[h];
        int nmo = Ca_so_->colspi()[h];
        if (nso == 0 || nmo == 0) continue;
        double** Dso = Da_so_->pointer(h);
        double** Cso = Ca_so_->pointer(h);
        double** Dmo = D->pointer(h);
        double** Temp = block_matrix(nso,nmo);

        C_DGEMM('N','N',nso,nmo,nso,1.0,Dso[0],nso,Cso[0],nmo,0.0,Temp[0],nmo);
        C_DGEMM('T','N',nmo,nmo,nso,1.0,Cso[0],nso,Temp[0],nmo,0.0,Dmo[0],nso);

        free_block(Temp);
    }
    return D;
}
SharedMatrix Prop::Db_mo()
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, asking for Db makes no sense");

    // MO D are nso x nso, zeros padding if nmo < nso
    SharedMatrix D = factory_->create_shared_matrix("Db (MO Basis)");
    for (int h = 0; h < D->nirrep(); h++) {
        int nso = D->colspi()[h];
        int nmo = Ca_so_->colspi()[h];
        if (nso == 0 || nmo == 0) continue;
        double** Dso = Db_so_->pointer(h);
        double** Cso = Cb_so_->pointer(h);
        double** Dmo = D->pointer(h);
        double** Temp = block_matrix(nso,nmo);

        C_DGEMM('N','N',nso,nmo,nso,1.0,Dso[0],nso,Cso[0],nmo,0.0,Temp[0],nmo);
        C_DGEMM('T','N',nmo,nmo,nso,1.0,Cso[0],nso,Temp[0],nmo,0.0,Dmo[0],nso);

        free_block(Temp);
    }
    return D;
}
void Prop::set_Da_so(SharedMatrix D)
{
    Da_so_ = D;
}
void Prop::set_Db_so(SharedMatrix D)
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Db makes no sense");

    Db_so_ = D;
}
void Prop::set_Ca_so(SharedMatrix C)
{
    Ca_so_ = C;
}
void Prop::set_Cb_so(SharedMatrix C)
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Cb makes no sense");

    Cb_so_ = C;
}
void Prop::set_Da_ao(SharedMatrix D)
{
    throw FeatureNotImplemented("Prop", "Advanced set methods not implemented", __FILE__, __LINE__);
}
void Prop::set_Db_ao(SharedMatrix D)
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Db makes no sense");

    throw FeatureNotImplemented("Prop", "Advanced set methods not implemented", __FILE__, __LINE__);
}
void Prop::set_Ca_ao(SharedMatrix C)
{
    throw FeatureNotImplemented("Prop", "Advanced set methods not implemented", __FILE__, __LINE__);
}
void Prop::set_Cb_ao(SharedMatrix C)
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Cb makes no sense");

    throw FeatureNotImplemented("Prop", "Advanced set methods not implemented", __FILE__, __LINE__);
}
void Prop::set_Da_mo(SharedMatrix D)
{
    throw FeatureNotImplemented("Prop", "Advanced set methods not implemented", __FILE__, __LINE__);
}
void Prop::set_Db_mo(SharedMatrix D)
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Db makes no sense");

    throw FeatureNotImplemented("Prop", "Advanced set methods not implemented", __FILE__, __LINE__);
}

OEProp::OEProp(boost::shared_ptr<Wavefunction> wfn) : Prop(wfn_)
{
    common_init();
}
OEProp::OEProp() : Prop(Process::environment.reference_wavefunction())
{
    common_init();
}
OEProp::~OEProp()
{
}
void OEProp::common_init()
{
}
void OEProp::print_header()
{
    fprintf(outfile, "\n OEPROP: One-electron properties/analyses.\n");
    fprintf(outfile, "  by Rob Parrish and Justin Turney.\n");
    fprintf(outfile, "  built on LIBMINTS.\n\n");
}
void OEProp::compute()
{
    print_header();
    if (tasks_.count("DIPOLE"))
        compute_dipole();
    if (tasks_.count("QUADRUPOLE"))
        compute_quadrupole();
    if (tasks_.count("OCTUPOLE"))
        compute_octupole();
    if (tasks_.count("HEXADECAPOLE"))
        compute_hexadecapole();
    if (tasks_.count("MO_EXTENTS"))
        compute_mo_extents();
    if (tasks_.count("MULLIKEN_CHARGES"))
        compute_mulliken_charges();
    if (tasks_.count("LOWDIN_CHARGES"))
        compute_lowdin_charges();
    if (tasks_.count("MAYER_INDICES"))
        compute_mayer_indices();
    if (tasks_.count("WIBERG_LOWDIN_INDICES"))
        compute_wiberg_lowdin_indices();
}
void OEProp::compute_dipole()
{
    boost::shared_ptr<Molecule> mol = basisset_->molecule();
    OperatorSymmetry dipsymm (1, mol, integral_, factory_);
    std::vector<SharedMatrix> so_dipole = dipsymm.create_matrices("SO Dipole");
    boost::shared_ptr<OneBodySOInt> sodOBI (integral_->so_dipole());
    sodOBI->compute(so_dipole);
    Vector3 de;
    SharedMatrix Da;
    SharedMatrix Db;

    if (restricted_) {
        Da = Da_so_;
        Db = Da;
    } else {
        Da = Da_so_;
        Db = Db_so_;
    }

    de[0] = Da->vector_dot(so_dipole[0]) + Db->vector_dot(so_dipole[0]);
    de[1] = Da->vector_dot(so_dipole[1]) + Db->vector_dot(so_dipole[1]);
    de[2] = Da->vector_dot(so_dipole[2]) + Db->vector_dot(so_dipole[2]);

    SharedVector ndip = mol->nuclear_dipole_contribution();
    de[0] += ndip->get(0, 0);
    de[1] += ndip->get(0, 1);
    de[2] += ndip->get(0, 2);

    fprintf(outfile," Dipole Moment: (a.u.)\n");
    fprintf(outfile,"     X: %10.4lf      Y: %10.4lf      Z: %10.4lf     Total: %10.4lf\n", \
       de[0], de[1], de[2], de.norm());
    fprintf(outfile, "\n");

    double dfac = _dipmom_au2debye;
    fprintf(outfile," Dipole Moment: (Debye)\n");
    fprintf(outfile,"     X: %10.4lf      Y: %10.4lf      Z: %10.4lf     Total: %10.4lf\n", \
       de[0]*dfac, de[1]*dfac, de[2]*dfac, de.norm()*dfac);
    fprintf(outfile, "\n");

    // Dipole components in Debye
    Process::environment.globals["DIPOLE X"] = de[0]*dfac;
    Process::environment.globals["DIPOLE Y"] = de[1]*dfac;
    Process::environment.globals["DIPOLE Z"] = de[2]*dfac;

    fflush(outfile);
}
void OEProp::compute_quadrupole()
{
    boost::shared_ptr<Molecule> mol = basisset_->molecule();
    SharedMatrix Da;
    SharedMatrix Db;

    if (restricted_) {
        Da = Da_so_;
        Db = Da;
    } else {
        Da = Da_so_;
        Db = Db_so_;
    }

    // Form the one-electron integral matrices from the matrix factory
    //    parameters: 1 (multipole order: 1=dipole, 2=quadrupole, etc.)
    OperatorSymmetry quadsymm(2, mol, integral_, factory_);

    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> so_Qpole = quadsymm.create_matrices("SO Quadrupole");

    // Form the one-electron integral objects from the integral factory
    boost::shared_ptr<OneBodySOInt> soqOBI(integral_->so_quadrupole());

    // Compute multipole moment integrals
    soqOBI->compute(so_Qpole);

    // Each multipole integral needs to be dotted with the SO Density
    SimpleVector qe(6);
    qe[0] = Da->vector_dot(so_Qpole[0]) + Db->vector_dot(so_Qpole[0]);
    qe[1] = Da->vector_dot(so_Qpole[1]) + Db->vector_dot(so_Qpole[1]);
    qe[2] = Da->vector_dot(so_Qpole[2]) + Db->vector_dot(so_Qpole[2]);
    qe[3] = Da->vector_dot(so_Qpole[3]) + Db->vector_dot(so_Qpole[3]);
    qe[4] = Da->vector_dot(so_Qpole[4]) + Db->vector_dot(so_Qpole[4]);
    qe[5] = Da->vector_dot(so_Qpole[5]) + Db->vector_dot(so_Qpole[5]);

    // Add in nuclear contribution
    SharedVector nquad = mol->nuclear_quadrupole_contribution();
    qe[0] += nquad->get(0, 0);
    qe[1] += nquad->get(0, 1);
    qe[2] += nquad->get(0, 2);
    qe[3] += nquad->get(0, 3);
    qe[4] += nquad->get(0, 4);
    qe[5] += nquad->get(0, 5);

    // Print multipole components
    double dfac = _dipmom_au2debye * _bohr2angstroms;
    fprintf(outfile, " Quadrupole Moment: (Debye Ang)\n");
    fprintf(outfile, "    XX: %10.4lf     YY: %10.4lf     ZZ: %10.4lf\n", \
       qe[0]*dfac, qe[3]*dfac, qe[5]*dfac);
    fprintf(outfile, "    XY: %10.4lf     XZ: %10.4lf     YZ: %10.4lf\n", \
       qe[1]*dfac, qe[2]*dfac, qe[4]*dfac);
    fprintf(outfile, "\n");

    double dtrace = (1.0 / 3.0) * (qe[0] + qe[3] + qe[5]);
    fprintf(outfile, " Traceless Quadrupole Moment: (Debye Ang)\n");
    fprintf(outfile, "    XX: %10.4lf     YY: %10.4lf     ZZ: %10.4lf\n", \
       (qe[0]-dtrace)*dfac, (qe[3]-dtrace)*dfac, (qe[5]-dtrace)*dfac);
    fprintf(outfile, "    XY: %10.4lf     XZ: %10.4lf     YZ: %10.4lf\n", \
       qe[1]*dfac, qe[2]*dfac, qe[4]*dfac);
    fprintf(outfile, "\n");

    // Quadrupole components in Debye Ang
    Process::environment.globals["QUADRUPOLE XX"] = qe[0]*dfac;
    Process::environment.globals["QUADRUPOLE YY"] = qe[3]*dfac;
    Process::environment.globals["QUADRUPOLE ZZ"] = qe[5]*dfac;
    Process::environment.globals["QUADRUPOLE XY"] = qe[1]*dfac;
    Process::environment.globals["QUADRUPOLE XZ"] = qe[2]*dfac;
    Process::environment.globals["QUADRUPOLE YZ"] = qe[4]*dfac;

    fflush(outfile);
}
void OEProp::compute_octupole()
{
    throw FeatureNotImplemented("OEProp::compute_octupole", "Octupole expectation value not implemented", __FILE__, __LINE__);

    fprintf(outfile, " OCTUPOLE ANALYSIS [a.u.]:\n\n");

    // Awesome code goes here.

    fflush(outfile);
}
void OEProp::compute_hexadecapole()
{
    throw FeatureNotImplemented("OEProp::compute_hexadecapole", "Hexadecapole expectation value not implemented", __FILE__, __LINE__);

    fprintf(outfile, " HEXADECAPOLE ANALYSIS [a.u.]:\n\n");

    // Awesome code goes here.

    fflush(outfile);
}
void OEProp::compute_mo_extents()
{
    throw FeatureNotImplemented("OEProp::compute_mo_extents", "MO Extents not implemented", __FILE__, __LINE__);

    fprintf(outfile, " MO Extents (<r^2>) [a.u.]:\n\n");

    boost::shared_ptr<Molecule> mol = basisset_->molecule();
    SharedMatrix Ca;
    SharedMatrix Cb;

    if (restricted_) {
        Ca = Ca_so_;
        Cb = Ca;
    } else {
        Ca = Ca_so_;
        Cb = Cb_so_;
    }

    // Form the one-electron integral matrices from the matrix factory parameters
    //    (multipole order: 1=dipole, 2=quadrupole, etc.)
    OperatorSymmetry diplsymm(1, mol, integral_, factory_);
    OperatorSymmetry quadsymm(2, mol, integral_, factory_);

    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> so_Dpole = diplsymm.create_matrices("SO Dipole");
    std::vector<SharedMatrix> so_Qpole = quadsymm.create_matrices("SO Quadrupole");

    // Form the one-electron integral objects from the integral factory
    boost::shared_ptr<OneBodySOInt> sodOBI(integral_->so_dipole());
    boost::shared_ptr<OneBodySOInt> soqOBI(integral_->so_quadrupole());

    // Compute multipole moment integrals
    sodOBI->compute(so_Dpole);
    soqOBI->compute(so_Qpole);


    //so_Dpole[0]->print();
    //so_Dpole[1]->print();
    //so_Dpole[2]->print();

    //so_Qpole[0]->print();
    //so_Qpole[3]->print();
    //so_Qpole[5]->print();

    SharedMatrix temp = factory_->create_shared_matrix("Temporary Matrix");

    // Compute directional expectation values for dipole contributions by molecular orbital
    //SharedMatrix CDCax = factory_->create_shared_matrix("CDxC alpha");
    //temp->gemm(false, false, 1.0, so_Dpole[0], Ca, 0.0);
    //CDCax->gemm(true, false, 1.0, Ca, temp, 0.0);
    //SharedMatrix CDCbx = factory_->create_shared_matrix("CDxC beta");
    //temp->gemm(false, false, 1.0, so_Dpole[0], Cb, 0.0);
    //CDCbx->gemm(true, false, 1.0, Cb, temp, 0.0);

    //SharedMatrix CDCay = factory_->create_shared_matrix("CDyC alpha");
    //temp->gemm(false, false, 1.0, so_Dpole[1], Ca, 0.0);
    //CDCay->gemm(true, false, 1.0, Ca, temp, 0.0);
    //SharedMatrix CDCby = factory_->create_shared_matrix("CDyC beta");
    //temp->gemm(false, false, 1.0, so_Dpole[1], Cb, 0.0);
    //CDCby->gemm(true, false, 1.0, Cb, temp, 0.0);

    SharedMatrix CDCaz = factory_->create_shared_matrix("CDzC alpha");
    temp->gemm(false, false, 1.0, so_Dpole[2], Ca, 0.0);
    CDCaz->gemm(true, false, 1.0, Ca, temp, 0.0);
    SharedMatrix CDCbz = factory_->create_shared_matrix("CDzC beta");
    temp->gemm(false, false, 1.0, so_Dpole[2], Cb, 0.0);
    CDCbz->gemm(true, false, 1.0, Cb, temp, 0.0);

    // Compute directional expectation values for quadrupole contributions by molecular orbital
    SharedMatrix CQCaxx = factory_->create_shared_matrix("CQxxC alpha");
    temp->gemm(false, false, 1.0, so_Qpole[0], Ca, 0.0);
    CQCaxx->gemm(true, false, 1.0, Ca, temp, 0.0);
    SharedMatrix CQCbxx = factory_->create_shared_matrix("CQxxC beta");
    temp->gemm(false, false, 1.0, so_Qpole[0], Cb, 0.0);
    CQCbxx->gemm(true, false, 1.0, Cb, temp, 0.0);

    SharedMatrix CQCayy = factory_->create_shared_matrix("CQyyC alpha");
    temp->gemm(false, false, 1.0, so_Qpole[3], Ca, 0.0);
    CQCayy->gemm(true, false, 1.0, Ca, temp, 0.0);
    SharedMatrix CQCbyy = factory_->create_shared_matrix("CQyyC beta");
    temp->gemm(false, false, 1.0, so_Qpole[3], Cb, 0.0);
    CQCbyy->gemm(true, false, 1.0, Cb, temp, 0.0);

    SharedMatrix CQCazz = factory_->create_shared_matrix("CQzzC alpha");
    temp->gemm(false, false,  1.0, so_Qpole[5], Ca, 0.0);
    CQCazz->gemm(true, false, 1.0, Ca, temp, 0.0);
    SharedMatrix CQCbzz = factory_->create_shared_matrix("CQzzC beta");
    temp->gemm(false, false,  1.0, so_Qpole[5], Cb, 0.0);
    CQCbzz->gemm(true, false, 1.0, Cb, temp, 0.0);

    // Accumulate registers for dipole
    double xi, yi, zi;
    char **labels = mol->irrep_labels();
    fprintf(outfile, " Molecular Orbital Polarities: (Bohr)\n");
    fprintf(outfile, "  Symmetry    MO        < x >      < y >      < z >          < r >\n");
    for (int h = 0; h < Ca->nirrep(); h++) {
       for (int i = 0; i < Ca->rowspi()[h]; i++) {

    //      xi = CDCax->get(h,i,i);
    //      yi = CDCay->get(h,i,i);
          zi = CDCaz->get(h,i,i);
          fprintf(outfile,"      %4s  %4d   %10.4lf %10.4lf %10.4lf     %10.4lf\n", \
             labels[h], i+1, xi, yi, zi, xi+yi+zi);
       }
    }
    fprintf(outfile, "\n");
    ////for(int h = 0; h < Ca->nirrep(); h++) free(labels[h]); free(labels);

    // Accumulate registers for quadrupole
    double x2i, y2i, z2i;
    ////char **labels = mol->irrep_labels();
    fprintf(outfile, " Molecular Orbital Extents: (Bohr^2)\n");
    fprintf(outfile, "  Symmetry    MO      < x^2 >    < y^2 >    < z^2 >        < r^2 >\n");
    for (int h = 0; h < Ca->nirrep(); h++) {
       for (int i = 0; i < Ca->rowspi()[h]; i++) {

          x2i = -CQCaxx->get(h,i,i);
          y2i = -CQCayy->get(h,i,i);
          z2i = -CQCazz->get(h,i,i);
          fprintf(outfile,"      %4s  %4d   %10.4lf %10.4lf %10.4lf     %10.4lf\n", \
             labels[h], i+1, x2i, y2i, z2i, x2i+y2i+z2i);
       }
    }
    fprintf(outfile, "\n");
    for(int h = 0; h < Ca->nirrep(); h++) free(labels[h]); free(labels);

    fflush(outfile);
}
void OEProp::compute_mulliken_charges()
{
    fprintf(outfile, " Mulliken Charges: (a.u.)\n");

    boost::shared_ptr<Molecule> mol = basisset_->molecule();

    double* Qa = new double[mol->natom()];
    double* PSa = new double[basisset_->nbf()];
    double suma = 0.0;

    double* Qb = new double[mol->natom()];
    double* PSb = new double[basisset_->nbf()];
    double sumb = 0.0;

    ::memset(Qa, '\0', mol->natom()*sizeof(double));
    ::memset(Qb, '\0', mol->natom()*sizeof(double));

    SharedMatrix Da;
    SharedMatrix Db;

//    Get the Density Matrices for alpha and beta spins

    if (restricted_) {
        Da = Da_ao();
        Db = Da;
    } else {
        Da = Da_ao();
        Db = Db_ao();
    }

//    Compute the overlap matrix

    boost::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    SharedMatrix S(new Matrix("S",basisset_->nbf(),basisset_->nbf()));
    overlap->compute(S);

//    Form the idempotent D*S matrix

    SharedMatrix PSam(new Matrix("PSa",basisset_->nbf(),basisset_->nbf()));
    PSam->gemm(false,false,1.0,Da,S,0.0);
    SharedMatrix PSbm(new Matrix("PSb",basisset_->nbf(),basisset_->nbf()));
    PSbm->gemm(false,false,1.0,Db,S,0.0);

//     Accumulate registers

    for (int mu = 0; mu < basisset_->nbf(); mu++) {
        PSa[mu] = PSam->get(0,mu,mu);
        PSb[mu] = PSbm->get(0,mu,mu);

        int shell = basisset_->function_to_shell(mu);
        int A = basisset_->shell_to_center(shell);

        Qa[A] += PSa[mu];
        Qb[A] += PSb[mu];

        suma += PSa[mu];
        sumb += PSb[mu];
    }

//    Print out the Mulliken populations and charges

    fprintf(outfile, "   Center  Symbol    Alpha    Beta     Spin     Total\n");
    double nuc = 0.0;
    for (int A = 0; A < mol->natom(); A++) {
        double Qs = Qa[A] - Qb[A];
        double Qt = mol->Z(A) - (Qa[A] + Qb[A]);
        fprintf(outfile,"   %5d    %2s    %8.5f %8.5f %8.5f %8.5f\n", A+1,mol->label(A).c_str(), \
            Qa[A], Qb[A], Qs, Qt);
        nuc += (double) mol->Z(A);
   }

    fprintf(outfile, "\n   Total alpha = %8.5f, Total beta = %8.5f, Total charge = %8.5f\n", \
        suma, sumb, nuc - suma - sumb);

//    Free memory
    delete[] Qa;
    delete[] Qb;
    delete[] PSa;
    delete[] PSb;

    fflush(outfile);
}
void OEProp::compute_lowdin_charges()
{
    fprintf(outfile, "\n\n Lowdin Charges [a.u.]:\n\n");

    boost::shared_ptr<Molecule> mol = basisset_->molecule();

    double* Qa = new double[mol->natom()];
    double suma = 0.0;

    double* Qb = new double[mol->natom()];
    double sumb = 0.0;

    ::memset(Qa, '\0', mol->natom()*sizeof(double));
    ::memset(Qb, '\0', mol->natom()*sizeof(double));

    SharedMatrix Da;
    SharedMatrix Db;
    SharedMatrix evecs(new Matrix("Eigenvectors of S matrix",basisset_->nbf(),basisset_->nbf()));
    SharedMatrix temp(new Matrix("Temporary matrix",basisset_->nbf(),basisset_->nbf()));
    SharedMatrix SDSa(new Matrix("S_12 * D * S_12 alpha matrix",basisset_->nbf(),basisset_->nbf()));
    SharedMatrix SDSb(new Matrix("S_12 * D * S_12 beta matrix",basisset_->nbf(),basisset_->nbf()));
    boost::shared_ptr<Vector> evals(new Vector(basisset_->nbf()));

//    Get the Density Matrices for alpha and beta spins

    if (restricted_) {
        Da = Da_ao();
        Db = Da;
    } else {
        Da = Da_ao();
        Db = Db_ao();
    }

//    Compute the overlap matrix

    boost::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    SharedMatrix S(new Matrix("S",basisset_->nbf(),basisset_->nbf()));
    overlap->compute(S);

//    Form the S^(1/2) matrix

    S->power(1.0/2.0);

//    Compute the S^(1/2)*D*S^(1/2) matrix

    temp->gemm(false,false,1.0,Da,S,0.0);
    SDSa->gemm(false,false,1.0,S,temp,0.0);
    temp->gemm(false,false,1.0,Db,S,0.0);
    SDSb->gemm(false,false,1.0,S,temp,0.0);

//    Accumulate AO populations for each atom

    for (int mu = 0; mu < basisset_->nbf(); mu++) {
        int shell = basisset_->function_to_shell(mu);
        int A = basisset_->shell_to_center(shell);

        Qa[A] += SDSa->get(0,mu,mu);
        Qb[A] += SDSb->get(0,mu,mu);

        suma += SDSa->get(0,mu,mu);
        sumb += SDSb->get(0,mu,mu);
    }

//    Print out the populations and charges

    fprintf(outfile, "   Center  Symbol    Alpha    Beta     Spin     Total\n");
    double nuc = 0.0;
    for (int A = 0; A < mol->natom(); A++) {
        double Qs = Qa[A] - Qb[A];
        double Qt = mol->Z(A) - (Qa[A] + Qb[A]);
        fprintf(outfile,"   %5d    %2s    %8.5f %8.5f %8.5f %8.5f\n", A+1,mol->label(A).c_str(), \
            Qa[A], Qb[A], Qs, Qt);
        nuc += (double) mol->Z(A);
    }

    fprintf(outfile, "\n  Total alpha = %8.5f, Total beta = %8.5f, Total charge = %8.5f\n", \
        suma, sumb, nuc - suma - sumb);

    delete[] Qa;
    delete[] Qb;

    fflush(outfile);
}
void OEProp::compute_mayer_indices()
{
    fprintf(outfile, "\n\n Mayer Bond Indices:\n\n");

    boost::shared_ptr<Molecule> mol = basisset_->molecule();

    int nbf = basisset_->nbf();

    SharedMatrix Da;      // Density matrix for alpha spin
    SharedMatrix Db;      // Density matrix for beta spin

    SharedMatrix DSa(new Matrix("D * S alpha matrix",nbf,nbf));
    SharedMatrix DSb(new Matrix("D * S beta matrix",nbf,nbf));

//    Get the Density Matrices for alpha and beta spins

    if (restricted_) {
        Da = Da_ao();
        Db = Da;
    } else {
        Da = Da_ao();
        Db = Db_ao();
    }

//    Compute the overlap matrix

    boost::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    SharedMatrix S(new Matrix("S matrix",nbf,nbf));
    overlap->compute(S);

//    Form the idempotent D*S matrix

    DSa->gemm(false,false,1.0,Da,S,0.0);
    DSb->gemm(false,false,1.0,Db,S,0.0);

//     Compute Mayer bond indices

    int natom = mol->natom();

    SharedMatrix MBI_total(new Matrix(natom,natom));
    SharedMatrix MBI_alpha;
    SharedMatrix MBI_beta;

    if (!restricted_) {
        MBI_alpha = SharedMatrix (new Matrix(natom,natom));
        MBI_beta = SharedMatrix (new Matrix(natom,natom));
    }

    for (int mu = 0; mu < nbf; mu++) {
        for (int nu = 0; nu < mu; nu++) {
            int shell_mu = basisset_->function_to_shell(mu);
            int shell_nu = basisset_->function_to_shell(nu);
            int atom_mu = basisset_->shell_to_center(shell_mu);
            int atom_nu = basisset_->shell_to_center(shell_nu);
            if (atom_mu == atom_nu) continue;

            double alpha = DSa->get(0, mu, nu) * DSa->get(0, nu, mu) ;
            double beta = DSb->get(0, mu, nu) * DSb->get(0, nu, mu);
            MBI_total->add(0, atom_mu, atom_nu, 2 * (alpha + beta));
            MBI_total->add(0, atom_nu, atom_mu, 2 * (alpha + beta));

            if (!restricted_) {
                MBI_alpha->add(0, atom_mu, atom_nu, 2 * (alpha));
                MBI_alpha->add(0, atom_nu, atom_mu, 2 * (alpha));
                MBI_beta->add(0, atom_mu, atom_nu, 2 * (beta));
                MBI_beta->add(0, atom_nu, atom_mu, 2 * (beta));
            }
        }
    }

//    Compute valences

    boost::shared_ptr<Vector> MBI_valence(new Vector(natom));

    for (int iat = 0; iat < natom; iat++) {
        for (int jat = 0; jat < natom; jat++) {
            double valence = MBI_valence->get(0, iat);
            MBI_valence->set(0, iat, valence + MBI_total->get(0, iat, jat));
        }
    }

//    Print out the bond index matrix and valences

//    Note: The computed Mayer bond indices (MBI) will be different from the MBI values computed using
//    some other program packages for the unrestricted case. The reason is that these programs left out
//    the spin density contribution in the MBI equation for the unrestricted wavefunctions. As the result,
//    the MBI value will be underestimated. For example, the MBI value for the H–H bond of H2+
//    is calculated to be 0.25 using the NBO program. The equation coded above gives the correct value of 0.5.
//    For reference, see IJQC 29 (1986) P. 73 and IJQC 29 (1986) P. 477.

//    A nicer output is needed ...

    if (restricted_) {
        MBI_total->print();
        fprintf(outfile, "  Atomic Valences: \n");
        MBI_valence->print();
    }
    else {
        fprintf(outfile, "  Total Bond Index: \n");
        MBI_total->print();
        fprintf(outfile, "  Alpha Contribution: \n");
        MBI_alpha->print();
        fprintf(outfile, "  Beta Contribution: \n");
        MBI_beta->print();
        fprintf(outfile, "  Atomic Valences: \n");
        MBI_valence->print();
    }

    fflush(outfile);
}
void OEProp::compute_wiberg_lowdin_indices()
{
    fprintf(outfile, "\n\n Wiberg Bond Indices using Orthogonal Lowdin Orbitals:\n\n");

//    We may wanna get rid of these if we have NAOs...

    boost::shared_ptr<Molecule> mol = basisset_->molecule();

    int nbf = basisset_->nbf();

    SharedMatrix Da;
    SharedMatrix Db;
    SharedMatrix evecs(new Matrix("Eigenvectors of S matrix",nbf,nbf));
    SharedMatrix temp(new Matrix("Temporary matrix",nbf,nbf));
    SharedMatrix SDSa(new Matrix("S_12 * D * S_12 alpha matrix",nbf,nbf));
    SharedMatrix SDSb(new Matrix("S_12 * D * S_12 beta matrix",nbf,nbf));
    boost::shared_ptr<Vector> evals(new Vector(nbf));

//    Get the Density Matrices for alpha and beta spins

    if (restricted_) {
        Da = Da_ao();
        Db = Da;
    } else {
        Da = Da_ao();
        Db = Db_ao();
    }

//    Compute the overlap matrix

    boost::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    SharedMatrix S(new Matrix("S",basisset_->nbf(),basisset_->nbf()));
    overlap->compute(S);

//    Form the S^(1/2) matrix

    S->diagonalize(evecs,evals);
    S->zero();
    for (int p = 0; p < basisset_->nbf(); ++p) S->set(0, p, p, sqrt(evals->get(0,p)));
    S->back_transform(evecs);

//    Compute the S^(1/2)*D*S^(1/2) matrix

    temp->gemm(false,false,1.0,Da,S,0.0);
    SDSa->gemm(false,false,1.0,S,temp,0.0);
    temp->gemm(false,false,1.0,Db,S,0.0);
    SDSb->gemm(false,false,1.0,S,temp,0.0);

//    Compute Wiberg bond indices

    int natom = mol->natom();

    SharedMatrix WBI_total(new Matrix(natom,natom));
    SharedMatrix WBI_alpha;
    SharedMatrix WBI_beta;

    if (!restricted_) {
        WBI_alpha = SharedMatrix (new Matrix(natom,natom));
        WBI_beta = SharedMatrix (new Matrix(natom,natom));
    }

    for (int mu = 0; mu < nbf; mu++) {
        for (int nu = 0; nu < mu; nu++) {
            int shell_mu = basisset_->function_to_shell(mu);
            int shell_nu = basisset_->function_to_shell(nu);
            int atom_mu = basisset_->shell_to_center(shell_mu);
            int atom_nu = basisset_->shell_to_center(shell_nu);
            if (atom_mu == atom_nu) continue;

            double alpha = SDSa->get(0, mu, nu) * SDSa->get(0, nu, mu) ;
            double beta = SDSb->get(0, mu, nu) * SDSb->get(0, nu, mu);
            WBI_total->add(0, atom_mu, atom_nu, 2 * (alpha + beta));
            WBI_total->add(0, atom_nu, atom_mu, 2 * (alpha + beta));

            if (!restricted_) {
                WBI_alpha->add(0, atom_mu, atom_nu, 2 * (alpha));
                WBI_alpha->add(0, atom_nu, atom_mu, 2 * (alpha));
                WBI_beta->add(0, atom_mu, atom_nu, 2 * (beta));
                WBI_beta->add(0, atom_nu, atom_mu, 2 * (beta));
            }
        }
    }

//    Compute valences

        boost::shared_ptr<Vector> WBI_valence(new Vector(natom));

        for (int iat = 0; iat < natom; iat++) {
            for (int jat = 0; jat < natom; jat++) {
                double valence = WBI_valence->get(0, iat);
                WBI_valence->set(0, iat, valence + WBI_total->get(0, iat, jat));
            }
        }

//    Print out the bond index matrix
//    A nicer output is needed ...

    if (restricted_) {
        WBI_total->print();
        fprintf(outfile, "  Atomic Valences: \n");
        WBI_valence->print();
    }
    else {
        fprintf(outfile, "  Total Bond Index: \n");
        WBI_total->print();
        fprintf(outfile, "  Alpha Contribution: \n");
        WBI_alpha->print();
        fprintf(outfile, "  Beta Contribution: \n");
        WBI_beta->print();
        fprintf(outfile, "  Atomic Valences: \n");
        WBI_valence->print();
    }

    fflush(outfile);
}
GridProp::GridProp(boost::shared_ptr<Wavefunction> wfn) : filename_("out.grid"), Prop(wfn)
{
    common_init();
}
GridProp::GridProp() : filename_("out.grid"), Prop(Process::environment.reference_wavefunction())
{
    common_init();
}
GridProp::~GridProp()
{
    reset();
    free_block(temp_tens_);
}
void GridProp::common_init()
{
    initialized_ = false;
    format_ = "DF3";

    n_[0] = 40;
    n_[1] = 40;
    n_[2] = 40;

    l_[0] = 5.0;
    l_[1] = 5.0;
    l_[2] = 5.0;

    o_[0] = 0.0;
    o_[1] = 0.0;
    o_[2] = 0.0;

    block_size_= 5000;

    irrep_offsets_[0] = 0;
    for (int h = 0; h < Ca_so_->nirrep() - 1; h++)
        irrep_offsets_[h + 1] = irrep_offsets_[h] + Ca_so_->colspi()[h];

    temp_tens_ = block_matrix(block_size_, basisset_->nbf());
}
void GridProp::add_alpha_mo(int irrep, int index)
{
    alpha_mos_.push_back(make_pair(irrep,index));
}
void GridProp::add_beta_mo(int irrep, int index)
{
    beta_mos_.push_back(make_pair(irrep,index));
}
void GridProp::add_basis_fun(int irrep, int index)
{
    basis_funs_.push_back(make_pair(irrep,index));
}
void GridProp::print_header()
{
    fprintf(outfile, "\n GRIDPROP: One-electron grid properties.\n");
    fprintf(outfile, "  by Rob Parrish and Justin Turney.\n");
    fprintf(outfile, "  built on LIBMINTS.\n\n");
}
double*** GridProp::block_grid(int nx, int ny, int nz)
{
    double*** grid = new double**[nx];

    double** pointers = new double*[nx*(unsigned long int)ny];

    double* memory = new double[nx*(unsigned long int)ny*nz];
    memset(static_cast<void*>(memory), '\0', sizeof(double)*nx*ny*nz);

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            pointers[i*(unsigned long int)ny + j] = &memory[i*(unsigned long int)ny*nz + j*(unsigned long int)nz];

    for (int i = 0; i < nx; i++)
        grid[i] = &pointers[i*(unsigned long int)ny];

    return grid;
}
void GridProp::free_grid(double*** grid)
{
    delete[] grid[0][0];
    delete[] grid[0];
    delete[] grid;
}
void GridProp::build_grid_overages(double over)
{
    boost::shared_ptr<Molecule> mol = basisset_->molecule();

    double min_x = mol->x(0);
    double min_y = mol->y(0);
    double min_z = mol->z(0);
    double max_x = mol->x(0);
    double max_y = mol->y(0);
    double max_z = mol->z(0);

    for (int A = 0; A < mol->natom(); A++) {
        if (mol->x(A) <= min_x)
            min_x = mol->x(A);
        if (mol->x(A) >= max_x)
            max_x = mol->x(A);
        if (mol->y(A) <= min_y)
            min_y = mol->y(A);
        if (mol->y(A) >= max_y)
            max_y = mol->y(A);
        if (mol->z(A) <= min_z)
            min_z = mol->z(A);
        if (mol->z(A) >= max_z)
            max_z = mol->z(A);
    }

    min_x -= over;
    min_y -= over;
    min_z -= over;
    max_x += over;
    max_y += over;
    max_z += over;

    o_[0] = 0.5*(min_x + max_x);
    o_[1] = 0.5*(min_y + max_y);
    o_[2] = 0.5*(min_z + max_z);
    l_[0] = (-min_x + max_x);
    l_[1] = (-min_y + max_y);
    l_[2] = (-min_z + max_z);

    caxis_[0] = 0.0;
    caxis_[1] = 1.0;

    build_grid();
}
void GridProp::build_grid()
{
    int nx = n_[0] + 1;
    int ny = n_[1] + 1;
    int nz = n_[2] + 1;

    grid_["x"] = block_grid(nx,ny,nz);
    grid_["y"] = block_grid(nx,ny,nz);
    grid_["z"] = block_grid(nx,ny,nz);

    double* x = new double[nx];
    double* y = new double[nx];
    double* z = new double[nx];

    double*** xg = grid_["x"];
    double*** yg = grid_["y"];
    double*** zg = grid_["z"];

    if (nx == 0)
        x[0] = 0.0;
    else
        for (int i = 0; i < nx; i++)
            x[i] = ((double) i) / (double (nx - 1));
    if (ny == 0)
        y[0] = 0.0;
    else
        for (int i = 0; i < ny; i++)
           y[i] = ((double) i) / (double (ny - 1));
    if (nz == 0)
        z[0] = 0.0;
    else
        for (int i = 0; i < nz; i++)
           z[i] = ((double) i) / (double (nz - 1));

    for (int i = 0; i < nx; i++) {
        x[i] = l_[0] * (x[i] - 0.5) + o_[0];
    }
    for (int i = 0; i < ny; i++) {
        y[i] = l_[1] * (y[i] - 0.5) + o_[1];
    }
    for (int i = 0; i < nz; i++) {
        z[i] = l_[2] * (z[i] - 0.5) + o_[2];
    }

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                xg[i][j][k] = x[i];
                yg[i][j][k] = y[j];
                zg[i][j][k] = z[k];
            }

    delete[] x;
    delete[] y;
    delete[] z;
}
void GridProp::allocate_arrays()
{
    int nx = n_[0] + 1;
    int ny = n_[1] + 1;
    int nz = n_[2] + 1;

    for (std::set<std::string>::iterator it = tasks_.begin(); it != tasks_.end(); it++) {
        if ((*it) == "MOS") {
            // TODO
        } else if ((*it) == "BASIS_FUNS") {
            // Also TODO
        } else {
            grid_[(*it)] = block_grid(nx,ny,nz);
        }
    }
}
void GridProp::compute()
{
#if 0
    reset();

    initialized_ = true;

    Da_ao_ = Da_ao();
    Ca_ao_ = Ca_ao();
    if (restricted_) {
        Db_ao_ = Da_ao_;
        Cb_ao_ = Ca_ao_;
    } else {
        Db_ao_ = Db_ao();
        Cb_ao_ = Cb_ao();
    }

    print_header();
    build_grid();
    allocate_arrays();

    int nx = n_[0] + 1;
    int ny = n_[1] + 1;
    int nz = n_[2] + 1;
    ULI ngrid = nx*(ULI)ny*nz;
    int nblock = ngrid / block_size_;
    if (ngrid % block_size_ != 0)
        nblock++;

    double*** xp = grid_["x"];
    double*** yp = grid_["y"];
    double*** zp = grid_["z"];

    // Basis points object (heavy lifting)
    points_ = boost::shared_ptr<BasisPoints>(new BasisPoints(basisset_, block_size_));
    if (tasks_.count("GAMMA_AA") || tasks_.count("GAMMA_BB") || tasks_.count("GAMMA_AB") \
        || tasks_.count("TAU_A") || tasks_.count("TAU_B"))
        points_->setToComputeGradients(true);

    // Grid block traversal object
    boost::shared_ptr<GridBlock> gridblock(new GridBlock());
    gridblock->setMaxPoints(block_size_);

    for (int block = 0; block < nblock; block++) {
        // Indexing
        int size = block_size_;
        if (block*(ULI)block_size_ >= ngrid)
            size = ngrid - block*(ULI)block_size_;

        ULI offset = block*(ULI)block_size_;

        // Line up gridblock pointers
        // Last xp is a dirty hack b/c w is not needed for points
        gridblock->setGrid(&xp[0][0][offset],&yp[0][0][offset],&zp[0][0][offset],&xp[0][0][offset]);
        gridblock->setTruePoints(size);

        // Compute basis functions/gradients
        points_->computePoints(gridblock);

        // Call compute routines
        if (tasks_.count("MOS"))
            compute_mos(gridblock, offset);
        if (tasks_.count("BASIS_FUNS"))
            compute_basis_funs(gridblock, offset);
        if (tasks_.count("RHO"))
            compute_rho(gridblock, &grid_["RHO"][0][0][offset]);
        if (tasks_.count("RHO_S"))
            compute_rho_s(gridblock, &grid_["RHO_S"][0][0][offset]);
        if (tasks_.count("RHO_A"))
            compute_rho_a(gridblock, &grid_["RHO_A"][0][0][offset]);
        if (tasks_.count("RHO_B"))
            compute_rho_b(gridblock, &grid_["RHO_B"][0][0][offset]);
        if (tasks_.count("GAMMA_AA"))
            compute_gamma_aa(gridblock, &grid_["GAMMA_AA"][0][0][offset]);
        if (tasks_.count("GAMMA_AB"))
            compute_gamma_ab(gridblock, &grid_["GAMMA_AB"][0][0][offset]);
        if (tasks_.count("GAMMA_BB"))
            compute_gamma_bb(gridblock, &grid_["GAMMA_BB"][0][0][offset]);
        if (tasks_.count("TAU_A"))
            compute_rho_b(gridblock, &grid_["TAU_A"][0][0][offset]);
        if (tasks_.count("TAU_B"))
            compute_rho_b(gridblock, &grid_["TAU_B"][0][0][offset]);
    }

    // ESP is special we think
    if (tasks_.count("ESP"))
        compute_ESP();

    if (format_ == "DF3")
        write_df3_grid();
    else
        write_data_grid();

#endif
}
void GridProp::write_data_grid()
{
    int nx = n_[0] + 1;
    int ny = n_[1] + 1;
    int nz = n_[2] + 1;

    for (std::map<std::string, double***>::iterator it = grid_.begin(); it != grid_.end(); it++) {
        std::string key = (*it).first;
        double*** data = (*it).second;

        /* Write it to a file */
        int i,j,k;
        std::string file = filename_ + "." + key + ".dat";
        FILE* fptr = fopen(file.c_str(),"w");
        fprintf(fptr,"%d %d %d\n\n", nx,ny,nz);
        for (k=0;k<nz;k++) {
           for (j=0;j<ny;j++) {
              for (i=0;i<nx;i++) {
                    fprintf(fptr,"%24.16f ", data[i][j][k]);
                }
            fprintf(fptr,"\n");
           }
           fprintf(fptr,"\n");
        }
        fclose(fptr);
    }
}
void GridProp::write_df3_grid()
{
    int nx = n_[0] + 1;
    int ny = n_[1] + 1;
    int nz = n_[2] + 1;

    for (std::map<std::string, double***>::iterator it = grid_.begin(); it != grid_.end(); it++) {
        std::string key = (*it).first;
        double*** data = (*it).second;

        double v;
        double themin = data[0][0][0];
        double themax = data[0][0][0];

        /* Write it to a file */
        std::string file = filename_ + "." + key + ".df3";
        FILE* fptr = fopen(file.c_str(),"w");
        fputc(nx >> 8,fptr);
        fputc(nx & 0xff,fptr);
        fputc(ny >> 8,fptr);
        fputc(ny & 0xff,fptr);
        fputc(nz >> 8,fptr);
        fputc(nz & 0xff,fptr);
        int i,j,k;
        for (k=0;k<nz;k++) {
           for (j=0;j<ny;j++) {
              for (i=0;i<nx;i++) {
                 if (data[i][j][k] > caxis_[1] )
                    v = 255;
                 else if (data[i][j][k] < caxis_[0])
                    v = 0;
                 else
                    v = 255 * (data[i][j][k]-caxis_[0])/(caxis_[1] - caxis_[0]);
                 fputc((int)v,fptr);
              }
           }
        }
        fclose(fptr);
    }
}
void GridProp::reset()
{
    if (!initialized_)
        return;

    // Free the points object
    //points_.reset();

    // Free the grids
    for (std::map<std::string, double***>::iterator it = grid_.begin(); it != grid_.end(); ++it) {
        free_grid((*it).second);
    }
}
#if 0
void GridProp::compute_mos(boost::shared_ptr<GridBlock> g, ULI offset)
{
    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
}
void GridProp::compute_basis_funs(boost::shared_ptr<GridBlock> g, ULI offset)
{
    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
}
void GridProp::compute_rho(boost::shared_ptr<GridBlock> g, double* results)
{
    int npoints = g->getTruePoints();
    int nbf = basisset_->nbf();
    double** points = points_->getPoints();
    double** Da = Da_ao_->pointer();
    double** Db = Db_ao_->pointer();

    // rho_a_
    // rho_a^Q = phi_m^Q * Da_mn * phi_n^Q
    C_DGEMM('N', 'N', npoints, nbf, nbf, 1.0, &points[0][0], nbf, &Da[0][0], nbf, \
        0.0, &temp_tens_[0][0], nbf);

    for (int Q = 0; Q < npoints; Q++) {
        results[Q] = C_DDOT(nbf, &temp_tens_[Q][0], 1, &points[Q][0], 1);
        //printf(" Q = %d, rho = %14.10E\n", Q, rho_a_[Q]);
    }

    if (!restricted_) {

        // rho_b^Q = phi_m^Q * Db_mn * phi_n^Q
        C_DGEMM('N', 'N', npoints, nbf, nbf, 1.0, &points[0][0], nbf, &Db[0][0], nbf, \
            0.0, &temp_tens_[0][0], nbf);

        for (int Q = 0; Q < npoints; Q++) {
            results[Q] += C_DDOT(nbf, &temp_tens_[Q][0], 1, &points[Q][0], 1);
            //printf(" Q = %d, rho = %14.10E\n", Q, rho_b_[Q]);
        }

    } else {
        C_DSCAL(npoints,2.0,results,1);
    }
}
void GridProp::compute_rho_s(boost::shared_ptr<GridBlock> g, double* results)
{
    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
}
void GridProp::compute_rho_a(boost::shared_ptr<GridBlock> g, double* results)
{
    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
}
void GridProp::compute_rho_b(boost::shared_ptr<GridBlock> g, double* results)
{
    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
}
void GridProp::compute_gamma_aa(boost::shared_ptr<GridBlock> g, double* results)
{
    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
}
void GridProp::compute_gamma_ab(boost::shared_ptr<GridBlock> g, double* results)
{
    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
}
void GridProp::compute_gamma_bb(boost::shared_ptr<GridBlock> g, double* results)
{
    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
}
void GridProp::compute_tau_a(boost::shared_ptr<GridBlock> g, double* results)
{
    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
}
void GridProp::compute_tau_b(boost::shared_ptr<GridBlock> g, double* results)
{
    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
}
#endif
void GridProp::compute_ESP()
{
    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);
}

}
