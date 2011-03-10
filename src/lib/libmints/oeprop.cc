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

using namespace boost;
using namespace psi;
using namespace std;

namespace psi {

Prop::Prop(shared_ptr<Wavefunction> wfn) : wfn_(wfn) 
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
    bool restricted_ = wfn_->restricted();
    
    integral_ = shared_ptr<IntegralFactory>(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));

    shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral_));
    AO2USO_ = shared_ptr<Matrix>(pet->aotoso());
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
shared_ptr<Matrix> Prop::Da_ao()
{
    int nao = basisset_->nbf();
    shared_ptr<Matrix> D = shared_ptr<Matrix>(new Matrix("Da (AO basis)", basisset_->nbf(), basisset_->nbf()));
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
shared_ptr<Matrix> Prop::Db_ao()
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, asking for Db makes no sense");

    int nao = basisset_->nbf();
    shared_ptr<Matrix> D = shared_ptr<Matrix>(new Matrix("Db (AO basis)", basisset_->nbf(), basisset_->nbf()));
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
shared_ptr<Matrix> Prop::Ca_ao()
{
    // TODO reorder by eigenvalue
    int nao = basisset_->nbf();
    shared_ptr<Matrix> C = shared_ptr<Matrix>(new Matrix("Ca (AO basis)", basisset_->nbf(), basisset_->nbf()));
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
shared_ptr<Matrix> Prop::Cb_ao()
{
    // TODO reorder by eigenvalue
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, asking for Cb makes no sense");

    int nao = basisset_->nbf();
    shared_ptr<Matrix> C = shared_ptr<Matrix>(new Matrix("Cb (AO basis)", basisset_->nbf(), basisset_->nbf()));
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
shared_ptr<Matrix> Prop::Da_mo()
{
    // MO D are nso x nso, zeros padding if nmo < nso
    shared_ptr<Matrix> D = factory_->create_shared_matrix("Da (MO Basis)");
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
shared_ptr<Matrix> Prop::Db_mo()
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, asking for Db makes no sense");

    // MO D are nso x nso, zeros padding if nmo < nso
    shared_ptr<Matrix> D = factory_->create_shared_matrix("Db (MO Basis)");
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
void Prop::set_Da_so(shared_ptr<Matrix> D) 
{
    Da_so_ = D;
}
void Prop::set_Db_so(shared_ptr<Matrix> D) 
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Db makes no sense");

    Db_so_ = D;
}
void Prop::set_Ca_so(shared_ptr<Matrix> C) 
{
    Ca_so_ = C;
}
void Prop::set_Cb_so(shared_ptr<Matrix> C) 
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Cb makes no sense");

    Cb_so_ = C;
}
void Prop::set_Da_ao(shared_ptr<Matrix> D) 
{
    throw FeatureNotImplemented("Prop", "Advanced set methods not implemented", __FILE__, __LINE__);    
}
void Prop::set_Db_ao(shared_ptr<Matrix> D) 
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Db makes no sense");

    throw FeatureNotImplemented("Prop", "Advanced set methods not implemented", __FILE__, __LINE__);    
}
void Prop::set_Ca_ao(shared_ptr<Matrix> C) 
{
    throw FeatureNotImplemented("Prop", "Advanced set methods not implemented", __FILE__, __LINE__);    
}
void Prop::set_Cb_ao(shared_ptr<Matrix> C) 
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Cb makes no sense");

    throw FeatureNotImplemented("Prop", "Advanced set methods not implemented", __FILE__, __LINE__);    
}
void Prop::set_Da_mo(shared_ptr<Matrix> D) 
{
    throw FeatureNotImplemented("Prop", "Advanced set methods not implemented", __FILE__, __LINE__);    
}
void Prop::set_Db_mo(shared_ptr<Matrix> D) 
{
    if (restricted_)
        throw PSIEXCEPTION("Wavefunction is restricted, setting Db makes no sense");

    throw FeatureNotImplemented("Prop", "Advanced set methods not implemented", __FILE__, __LINE__);    
}

OEProp::OEProp(shared_ptr<Wavefunction> wfn) : Prop(wfn_) 
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
}
void OEProp::compute_dipole()
{
    throw FeatureNotImplemented("OEProp::compute_dipole", "Dipole expectation value not implemented", __FILE__, __LINE__);    

    fprintf(outfile, " DIPOLE ANALYSIS [a.u.]:\n\n");

    // Awesome code goes here. 

    fflush(outfile);
}
void OEProp::compute_quadrupole()
{
    throw FeatureNotImplemented("OEProp::compute_quadrupole", "Quadrupole expectation value not implemented", __FILE__, __LINE__);    

    fprintf(outfile, " QUADRUPOLE ANALYSIS [a.u.]:\n\n");

    // Awesome code goes here. 

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

    // Awesome code goes here. 

    fflush(outfile);
}
void OEProp::compute_mulliken_charges()
{
    fprintf(outfile, " Mulliken Charges [a.u.]:\n\n");

    // OK, the world's worst Mulliken charge analysis (only works for RHF)
    // but it's 4 AM, and it gives you the idea
    // The finished version should handle R/U, do spin density, and perform a pop analysis in addition to charges

    shared_ptr<Molecule> mol = basisset_->molecule();

    double sum = 0.0; 
    double* Q = new double[mol->natom()];
    double* PS = new double[basisset_->nbf()];
    for (int A = 0; A < mol->natom(); A++) {
        Q[A] = (double) mol->Z(A);
        sum += (double) mol->Z(A);
    }
    
    shared_ptr<Matrix> D = Da_ao();
    shared_ptr<OneBodyAOInt> overlap (integral_->ao_overlap());
    shared_ptr<Matrix> S(new Matrix("Overlap Matrix", basisset_->nbf(), basisset_->nbf()));
    overlap->compute(S);

    double** Sp = S->pointer();
    double** Dp = D->pointer();
    
    for (int mu = 0; mu < basisset_->nbf(); mu++) {
        PS[mu] = 2.0*C_DDOT(basisset_->nbf(), Dp[mu], 1, Sp[mu], 1);
        Q[basisset_->shell(basisset_->function_to_shell(mu))->ncenter()] -= PS[mu]; 
        sum -= PS[mu];
    }
       
    fprintf(outfile, "   Center  Symbol  Charge \n");
    fprintf(outfile, "  ---------------------------\n");
    for (int A = 0; A < mol->natom(); A++)
        fprintf(outfile,"   %6d    %2s    %9.6f\n", A+1,mol->symbol(A).c_str(),Q[A]);
    fprintf(outfile,"\n   Total Charge: %9.6f\n\n",sum);
    fflush(outfile);

    delete[] Q;
    delete[] PS;
}
void OEProp::compute_lowdin_charges()
{
    throw FeatureNotImplemented("OEProp::compute_lowdin_charges", "Lowdin charges not implemented", __FILE__, __LINE__);    

    fprintf(outfile, " LOWDIN CHARGES [a.u.]:\n\n");

    // Awesome code goes here. 

    fflush(outfile);
}

GridProp::GridProp(shared_ptr<Wavefunction> wfn) : filename_("out.grid"), Prop(wfn)
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
    shared_ptr<Molecule> mol = basisset_->molecule();
    
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
    points_ = shared_ptr<BasisPoints>(new BasisPoints(basisset_, block_size_));    
    if (tasks_.count("GAMMA_AA") || tasks_.count("GAMMA_BB") || tasks_.count("GAMMA_AB") \
        || tasks_.count("TAU_A") || tasks_.count("TAU_B"))
        points_->setToComputeGradients(true);

    // Grid block traversal object
    shared_ptr<GridBlock> gridblock(new GridBlock());
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
    points_.reset();

    // Free the grids
    for (std::map<std::string, double***>::iterator it = grid_.begin(); it != grid_.end(); ++it) {
        free_grid((*it).second);
    }
}
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
void GridProp::compute_ESP()
{
    throw FeatureNotImplemented("GridProp", "This property not implemented", __FILE__, __LINE__);    
}

}
