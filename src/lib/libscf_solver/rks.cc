#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <cstring>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/vector3.h>
#include <libmints/gridblock.h>
#include <libmints/basisset.h>
#include <libmints/twobody.h>
#include <libmints/basispoints.h>
#include <libmints/properties.h>
#include <libmints/molecule.h>
//#include "integrator.h"
#include <libfunctional/functional.h>
#include <libmints/matrix.h>
#include "rks.h"


using namespace std;
using namespace psi;

namespace psi { namespace scf {

RKS::RKS(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt) : RHF(options, psio, chkpt)
{
    //x_functional_ = Functional::createFunctional(options.get_str("X_FUNCTIONAL"),options.get_int("N_BLOCK"));
    //c_functional_ = Functional::createFunctional(options.get_str("C_FUNCTIONAL"),options.get_int("N_BLOCK"));
/**
    integrator_ = Integrator::createIntegrator(molecule_,options);
    V_ = SharedMatrix (factory_.create_matrix("V"));
    properties_ = SharedProperties(Properties::constructProperties(basisset_,options.get_int("N_BLOCK")));
**/
    /**

    if (x_functional_->isGGA() || c_functional_-> isGGA() )
    {
        fprintf(outfile,"\n  Computing Point Gradients, X is %s, C is %s\n",(x_functional_->isGGA())?"GGA":"Not GGA",(c_functional_->isGGA())?"GGA":"Not GGA");
        properties_->setToComputeDensityGradient(true); //Need this to be able to get to \nabla \rho
    }

    fprintf(outfile,"  \n");
    fprintf(outfile,"  Exchange Functional Name: %s\n",x_functional_->getName().c_str());
    fprintf(outfile,"  Exchange Functional Description:\n  %s\n",x_functional_->getDescription().c_str());
    fprintf(outfile,"  Exchange Functional Citation:\n  %s\n",x_functional_->getCitation().c_str());
    fprintf(outfile,"  Exchange Functional Parameters:\n%s\n",x_functional_->getParametersString().c_str());
    fprintf(outfile,"  \n");

    if (options.get_bool("TEST_FUNCTIONAL") ) {
        int n_test = properties_->get_RKS_GGA_Testbed_Size();
        properties_->get_RKS_GGA_Testbed(); //Please set block size bigger than 9 or risk a segfault!
        x_functional_->computeRKSFunctional(properties_);

        const double* rhoa = properties_->getDensity();
        const double* sigmaa = properties_->getDensityGradientSquared();

        double *zk = x_functional_->getFunctional();
        double *vrhoa = x_functional_->getV_RhoA();
        double *vsigmaaa = x_functional_->getV_GammaAA();

        fprintf(outfile, "  Testing Exchange Functional:\n");
        for (int k = 0; k<n_test; k++) {
            fprintf(outfile, "   rhoa = %8.2E, sigmaa = %8.2E\n", rhoa[k], sigmaa[k]);
            fprintf(outfile, "   zk = %15.11E\n",2.0*zk[k]);
            fprintf(outfile, "   vrhoa = %15.11E\n",vrhoa[k]);
            if (x_functional_->isGGA())
                fprintf(outfile, "   vsimgaaa = %15.11E\n",vsigmaaa[k]);
        }
        fprintf(outfile, "  \n");
    }

    fprintf(outfile,"  Correlation Functional Name: %s\n",c_functional_->getName().c_str());
    fprintf(outfile,"  Correlation Functional Description:\n  %s\n",c_functional_->getDescription().c_str());
    fprintf(outfile,"  Correlation Functional Citation:\n  %s\n",c_functional_->getCitation().c_str());
    fprintf(outfile,"  Correlation Functional Parameters:\n%s\n",c_functional_->getParametersString().c_str());
    fprintf(outfile,"  \n");
    fprintf(outfile,"%s\n",(integrator_->getString()).c_str());

    if (options.get_bool("TEST_FUNCTIONAL") ) {
        int n_test = properties_->get_RKS_GGA_Testbed_Size();
        properties_->get_RKS_GGA_Testbed(); //Please set block size bigger than 9 or risk a segfault!
        c_functional_->computeRKSFunctional(properties_);

        const double* rhoa = properties_->getDensity();
        const double* sigmaa = properties_->getDensityGradientSquared();

        double *zk = c_functional_->getFunctional();
        double *vrhoa = c_functional_->getV_RhoA();
        double *vsigmaaa = c_functional_->getV_GammaAA();

        fprintf(outfile, "  Testing Correlation Functional:\n");
        for (int k = 0; k<n_test; k++) {
            fprintf(outfile, "   rhoa = %8.2E, sigmaa = %8.2E\n", rhoa[k], sigmaa[k]);
            fprintf(outfile, "   zk = %15.11E\n",2.0*zk[k]);
            fprintf(outfile, "   vrhoa = %15.11E\n",vrhoa[k]);
            if (c_functional_->isGGA())
                fprintf(outfile, "   vsimgaaa = %15.11E\n",vsigmaaa[k]);
        }
        fprintf(outfile, "  \n");
    }
    **/
}

RKS::RKS(Options& options, shared_ptr<PSIO> psio) : RHF(options, psio)
{
    //x_functional_ = Functional::createFunctional(options.get_str("X_FUNCTIONAL"),options.get_int("N_BLOCK"));
    //c_functional_ = Functional::createFunctional(options.get_str("C_FUNCTIONAL"),options.get_int("N_BLOCK"));
/**
    integrator_ = Integrator::createIntegrator(molecule_,options);
    V_ = SharedMatrix (factory_.create_matrix("V"));
    properties_ = SharedProperties(Properties::constructProperties(basisset_,options.get_int("N_BLOCK")));
**/
    /**

    if (x_functional_->isGGA() || c_functional_-> isGGA() )
    {
        fprintf(outfile,"\n  Computing Point Gradients, X is %s, C is %s\n",(x_functional_->isGGA())?"GGA":"Not GGA",(c_functional_->isGGA())?"GGA":"Not GGA");
        properties_->setToComputeDensityGradient(true); //Need this to be able to get to \nabla \rho
    }

    fprintf(outfile,"  \n");
    fprintf(outfile,"  Exchange Functional Name: %s\n",x_functional_->getName().c_str());
    fprintf(outfile,"  Exchange Functional Description:\n  %s\n",x_functional_->getDescription().c_str());
    fprintf(outfile,"  Exchange Functional Citation:\n  %s\n",x_functional_->getCitation().c_str());
    fprintf(outfile,"  Exchange Functional Parameters:\n%s\n",x_functional_->getParametersString().c_str());
    fprintf(outfile,"  \n");

    if (options.get_bool("TEST_FUNCTIONAL") ) {
        int n_test = properties_->get_RKS_GGA_Testbed_Size();
        properties_->get_RKS_GGA_Testbed(); //Please set block size bigger than 9 or risk a segfault!
        x_functional_->computeRKSFunctional(properties_);

        const double* rhoa = properties_->getDensity();
        const double* sigmaa = properties_->getDensityGradientSquared();

        double *zk = x_functional_->getFunctional();
        double *vrhoa = x_functional_->getV_RhoA();
        double *vsigmaaa = x_functional_->getV_GammaAA();

        fprintf(outfile, "  Testing Exchange Functional:\n");
        for (int k = 0; k<n_test; k++) {
            fprintf(outfile, "   rhoa = %8.2E, sigmaa = %8.2E\n", rhoa[k], sigmaa[k]);
            fprintf(outfile, "   zk = %15.11E\n",2.0*zk[k]);
            fprintf(outfile, "   vrhoa = %15.11E\n",vrhoa[k]);
            if (x_functional_->isGGA())
                fprintf(outfile, "   vsimgaaa = %15.11E\n",vsigmaaa[k]);
        }
        fprintf(outfile, "  \n");
    }

    fprintf(outfile,"  Correlation Functional Name: %s\n",c_functional_->getName().c_str());
    fprintf(outfile,"  Correlation Functional Description:\n  %s\n",c_functional_->getDescription().c_str());
    fprintf(outfile,"  Correlation Functional Citation:\n  %s\n",c_functional_->getCitation().c_str());
    fprintf(outfile,"  Correlation Functional Parameters:\n%s\n",c_functional_->getParametersString().c_str());
    fprintf(outfile,"  \n");
    fprintf(outfile,"%s\n",(integrator_->getString()).c_str());

    if (options.get_bool("TEST_FUNCTIONAL") ) {
        int n_test = properties_->get_RKS_GGA_Testbed_Size();
        properties_->get_RKS_GGA_Testbed(); //Please set block size bigger than 9 or risk a segfault!
        c_functional_->computeRKSFunctional(properties_);

        const double* rhoa = properties_->getDensity();
        const double* sigmaa = properties_->getDensityGradientSquared();

        double *zk = c_functional_->getFunctional();
        double *vrhoa = c_functional_->getV_RhoA();
        double *vsigmaaa = c_functional_->getV_GammaAA();

        fprintf(outfile, "  Testing Correlation Functional:\n");
        for (int k = 0; k<n_test; k++) {
            fprintf(outfile, "   rhoa = %8.2E, sigmaa = %8.2E\n", rhoa[k], sigmaa[k]);
            fprintf(outfile, "   zk = %15.11E\n",2.0*zk[k]);
            fprintf(outfile, "   vrhoa = %15.11E\n",vrhoa[k]);
            if (c_functional_->isGGA())
                fprintf(outfile, "   vsimgaaa = %15.11E\n",vsigmaaa[k]);
        }
        fprintf(outfile, "  \n");
    }
    **/
}

RKS::~RKS()
{

}

double RKS::compute_energy()
{
/**
    bool converged = false, diis_iter = false;
    int iteration = 0;

    // Do the initial work to get the iterations started.
    //form_multipole_integrals();  // handled by HF class
    form_H();

    H_->print(outfile);

    form_Shalf();

    S_->print(outfile);

    load_or_compute_initial_C();

    C_->print(outfile);
    D_->print(outfile);

    if (scf_type_ == "PK")
        form_PK();
    else if (scf_type_ == "DF" || scf_type_ == "CD" || scf_type_ == "1C_CD" )
        form_B();

    fprintf(outfile, "                                  Total Energy            Delta E              Density RMS\n\n");
    // SCF iterations
    do {
        iteration++;

        Dold_->copy(D_);  // Save previous density
        Eold_ = E_;       // Save previous energy

        form_J(); //J is always needed, and sometimes you get
                  //K for free-ish with that
        J_->print(outfile);
        if (functional_->isHybrid()) {
            form_K(); //Sometimes, you really need K too
            K_->print(outfile);
        }

        form_V(); //Whoa, there's that Kohn-Sham stuff
        V_->print(outfile);

        form_F();
        Fa_->print(outfile);

        if (diis_enabled_)
            save_fock();

        E_ = compute_E();

        if (diis_enabled_ == true && iteration >= min_diis_vectors_) {
            diis();
            diis_iter = true;
        } else {
            diis_iter = false;
        }

        form_C();
        form_D();

        C_->print(outfile);
        D_->print(outfile);

        converged = test_convergency();
        fprintf(outfile, "  @RKS iteration %3d energy: %20.14f    %20.14f %20.14f %s\n", iteration, E_, E_ - Eold_, Drms_, diis_iter == false ? " " : "DIIS");
        fflush(outfile);
    } while (!converged && iteration < maxiter_);

    if (converged) {
        fprintf(outfile, "\n  Energy converged.\n");
        fprintf(outfile, "\n  @RKS Final Energy: %20.14f", E_);
        if (perturb_h_) {
            fprintf(outfile, " with %f perturbation", lambda_);
        }
        fprintf(outfile, "\n");
        fprintf(outfile,"\n  @RKS Numerical Density: %14.10f\n",densityCheck_);
        fprintf(outfile,"  @RKS Numerical Dipole: <%14.10f,%14.10f,%14.10f>\n",dipoleCheckX_,dipoleCheckY_,dipoleCheckZ_);

        save_information();
    } else {
        fprintf(outfile, "\n  Failed to converged.\n");
        fprintf(outfile,"\n  @RKS Numerical Density: %14.10f\n",densityCheck_);
        fprintf(outfile,"  @RKS Numerical Dipole: <%14.10f,%14.10f,%14.10f>\n",dipoleCheckX_,dipoleCheckY_,dipoleCheckZ_);
        E_ = 0.0;
    }
    if (options_.get_bool("SAVE_NUMERICAL_GRID"))
    {
        fprintf(outfile,"  Saving DFT Numerical Grid\n");
        save_DFT_grid();
    }
    if (save_grid_)
    {
        //DOWN FOR MAINTENANCE
        //fprintf(outfile,"  Saving Cartesian Grid\n");
        //save_RHF_grid(options_, basisset_, D_, C_);
    }
    if (scf_type_ == "DF" || scf_type_ == "CD" || scf_type_ == "1C_CD")
    {
        free_B();
    }
    // Compute the final dipole.
    compute_multipole();

        //fprintf(outfile,"\nComputation Completed\n");
        //fflush(outfile);
    return E_;
**/
    return 0.0;
}

void RKS::form_J()
{
    if (scf_type_ == "PK") {
        //form_J_from_PK();
    }
    else if (scf_type_ == "DIRECT") {
        form_J_and_K_from_direct_integrals();
    }
    else if (scf_type_ == "DF" || scf_type_ == "CD" || scf_type_ == "1C_CD") {
        form_J_from_RI();
    }
    else if (scf_type_ == "OUT_OF_CORE" ){
        //form_J_and_K_disk();
    }
}
void RKS::form_K()
{
    if (scf_type_ == "PK") {
        //form_K_from_PK();
    }
    else if (scf_type_ == "OUT_OF_CORE" || scf_type_ == "DIRECT") {
        //Already formed
    }
    else if (scf_type_ == "DF" || scf_type_ == "CD" || scf_type_ == "1C_CD") {
        form_K_from_RI();
    }
    else {
        //Already formed
    }
}
double RKS::compute_E()
{
    double one_electron_energy = 2.0*D_->vector_dot(H_);
    double J_energy = D_->vector_dot(J_);
    double K_energy = 0.0;
    if (functional_->isHybrid()) {
        K_energy = functional_->getExactExchange()*D_->vector_dot(K_);
    }
    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_energy;
    Etotal += J_energy;
    Etotal += K_energy;
    Etotal += functional_energy_;

    //fprintf(outfile,"  One Electron Energy: %14.10f\n",one_electron_energy);
    //fprintf(outfile,"  Classical Coulomb Energy: %14.10f\n",J_energy);
    //fprintf(outfile,"  Functional Energy: %14.10f\n",functional_energy_);

    return Etotal;

}
void RKS::form_F()
{
    Fa_->copy(H_);
    J_->scale(2.0);
        Fa_->add(J_);
    //J_->scale(0.5);
    if (functional_->isHybrid()) {
        K_->scale(-functional_->getExactExchange());
        Fa_->add(K_);
    }
    //V_->scale(1.0);
    Fa_->add(V_);
    #ifdef _DEBUG
    if (debug_){
        Fa_->print(outfile);
    }
    #endif
    Fa_->scale(2.0);
}
void RKS::form_V()
{
/**
    V_->zero();

    bool GGA = functional_->isGGA();

    functional_energy_ = 0.0;
    densityCheck_ = 0.0;
    dipoleCheckX_ = 0.0;
    dipoleCheckY_ = 0.0;
    dipoleCheckZ_ = 0.0;

    //LSDA Setup
    double val;

    double *x_fun = functional_->getFunctional(0)->getFunctional();
    double *x_funGrad = functional_->getFunctional(0)->getV_RhoA();

    const double* rhog = properties_->getDensity();

    double **basis_points = properties_->getPoints();
    //End LSDA Setup

    //GGA Setup (Pointers will only be accessed if functional is GGA)
    double GGA_val;
    double fun_x, fun_y, fun_z;

    double *x_funGradAA = functional_->getFunctional(0)->getV_GammaAA();
    double *x_funGradAB = functional_->getFunctional(0)->getV_GammaAB();

    const double *rho_x = properties_->getDensityX();
    const double *rho_y = properties_->getDensityY();
    const double *rho_z = properties_->getDensityZ();

    double **basis_points_x = properties_->getGradientsX();
    double **basis_points_y = properties_->getGradientsY();
    double **basis_points_z = properties_->getGradientsZ();
    //End of GGA setup

    int nirreps = V_->nirrep();
    int* opi = V_->rowspi();
    double check; //numerical density contribution

    for (integrator_->reset(); !integrator_->isDone(); ) {
        SharedGridBlock q = integrator_->getNextBlock();
        int ntrue = q->getTruePoints();
        double* xg = q->getX();
        double* yg = q->getY();
        double* zg = q->getZ();
        double* wg = q->getWeights();

        properties_->computeProperties(q,D_,C_);
        functional_->getFunctional(0)->computeRKSFunctional(properties_);

        for (int grid_index = 0; grid_index<ntrue; grid_index++) {
        //BEGIN LOOP OVER GRID POINTS

            functional_energy_ += 2.0*x_fun[grid_index]*wg[grid_index];

            //GGA contribution (2 \diff f / \diff s_aa \nabla \rho + \diff f / \diff s_ab \nabla \rho)
            fun_x = 0.0;
            fun_y = 0.0;
            fun_z = 0.0;
            if (GGA) {
                fun_x += (2.0*x_funGradAA[grid_index]+x_funGradAB[grid_index])*rho_x[grid_index];
                fun_y += (2.0*x_funGradAA[grid_index]+x_funGradAB[grid_index])*rho_y[grid_index];
                fun_z += (2.0*x_funGradAA[grid_index]+x_funGradAB[grid_index])*rho_z[grid_index];
            }

            //fprintf(outfile,"  Point: <%14.10f,%14.10f,%14.10f>, w = %14.10f, rho = %14.10f\n",xg[grid_index],yg[grid_index],zg[grid_index],wg[grid_index],rhog[grid_index]);
            //fprintf(outfile,"    x_fun = %14.10f, c_fun = %14.10f, d_x_fun = %14.10f, d_c_fun = %14.10f\n",x_fun[grid_index],c_fun[grid_index],x_funGrad[grid_index],c_funGrad[grid_index]);
            if (GGA) {
                fprintf(outfile,"    del_rho: <%14.10f,%14.10f,%14.10f>\n",rho_x[grid_index],rho_y[grid_index],rho_z[grid_index]);
                fprintf(outfile,"    del_fun: <%14.10f,%14.10f,%14.10f>\n",fun_x,fun_y,fun_z);
            }
            int h_offset = 0;
            for (int h = 0; h<nirreps; h_offset+=opi[h],h++) {
                for (int i = 0; i<opi[h];i++) {
                    for (int j = 0; j<=i; j++) {
                        //fprintf(outfile,"   (%4d,%4d), phi_a = %14.10f, phi_b = %14.10f\n",i,j,basis_points[i],basis_points[j]);
                        //LSDA Contribution
                        val = wg[grid_index]*basis_points[grid_index][i+h_offset]*(x_funGrad[grid_index])*basis_points[grid_index][j+h_offset];
                        V_->add(h,i,j,val);

                        if (i!=j)
                            V_->add(h,j,i,val);

                        //GGA Contribution
                        if (GGA) {
                            GGA_val = 0.0;
                            GGA_val += fun_x*(basis_points_x[grid_index][i+h_offset]*basis_points[grid_index][j+h_offset]+basis_points[grid_index][i+h_offset]*basis_points_x[grid_index][j+h_offset]);
                            GGA_val += fun_y*(basis_points_y[grid_index][i+h_offset]*basis_points[grid_index][j+h_offset]+basis_points[grid_index][i+h_offset]*basis_points_y[grid_index][j+h_offset]);
                            GGA_val += fun_z*(basis_points_z[grid_index][i+h_offset]*basis_points[grid_index][j+h_offset]+basis_points[grid_index][i+h_offset]*basis_points_z[grid_index][j+h_offset]);
                            GGA_val *= wg[grid_index];
                            V_->add(h,i,j,GGA_val);
                            if (i!=j)
                                V_->add(h,j,i,GGA_val);
                        }

                    }
                }
            }
            //Check Numerical Density (Gives an idea of grid error)
            check = wg[grid_index]*rhog[grid_index];
            densityCheck_+=2.0*check;
            dipoleCheckX_+=-2.0*check*xg[grid_index];
            dipoleCheckY_+=-2.0*check*yg[grid_index];
            dipoleCheckZ_+=-2.0*check*zg[grid_index];
        //END LOOP OVER GRID POINTS
        }
    }
    //V_->print(outfile);
    //fprintf(outfile, "  Density Check: %14.10f\n",densityCheck_);
    //fprintf(outfile, "  Dipole  Check: <%14.10f,%14.10f,%14.10f>\n",dipoleCheckX_,dipoleCheckY_,dipoleCheckZ_);
    //fprintf(outfile, "  Functional Check: %14.10f\n",check_functional);

**/
}
void RKS::save_DFT_grid()
{
/**
 FILE * grid_file = fopen((options_.get_str("NUMERICAL_GRID_FILENAME")).c_str(),"w");
 shared_ptr<Molecule> mol = basisset_->molecule();
    fprintf(grid_file,"%d\n",mol->natom());
    fprintf(grid_file,"x,y,z,Z\n");
    for (int i=0; i<mol->natom(); i++)
        fprintf(grid_file,"%14.10f,%14.10f,%14.10f,%d\n",mol->x(i),mol->y(i), mol->z(i), mol->Z(i));

    //IntegrationPoint q;
    for (integrator_->reset(); !integrator_->isDone(); )
    {
        //q = integrator_->getNextPoint();
   //     Vector3 v = q.point;
     //   double w = q.weight;
       // fprintf(grid_file,"%14.10f, %14.10f, %14.10f, %14.10f\n",v[0],v[1], v[2], w);
    }


    fclose(grid_file);
**/
}
}}
