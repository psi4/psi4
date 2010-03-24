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
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/vector3.h>
#include <libmints/basisset.h>
#include <libmints/twobody.h>
#include <libmints/basispoints.h>
#include <libmints/properties.h>
#include "integrator.h"
#include "functional.h"
#include "functionalfactory.h"
#include <libmints/matrix.h>
#include "rks.h"


using namespace std;
using namespace psi;

namespace psi { namespace scf {
    
RKS::RKS(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt) : RHF(options, psio, chkpt)
{
	FunctionalFactory fact;
	functional_ = SharedFunctional(fact.getFunctional(options.get_str("FUNCTIONAL")));
	integrator_ = SharedIntegrator(Integrator::createIntegrator(molecule_,options));
	V_ = SharedMatrix (factory_.create_matrix("V")); 
	properties_ = SharedProperties(Properties::constructProperties(basisset_));

	fprintf(outfile,"  \n");
	fprintf(outfile,"  In RKS by Rob Parrish\n");
	fprintf(outfile,"  \n");
	fprintf(outfile,"  Functional Name: %s\n",functional_->getName().c_str());
	fprintf(outfile,"  Functional Description:\n%s\n",functional_->getDescription().c_str());
	fprintf(outfile,"  Functional Citation:\n%s\n",functional_->getCitation().c_str());
	fprintf(outfile,"  \n");
	fprintf(outfile,"%s",(integrator_->getString()).c_str());

}

RKS::~RKS()
{

}

double RKS::compute_energy()
{
    bool converged = false, diis_iter = false;
    int iteration = 0;

    // Do the initial work to get the iterations started.
    //form_multipole_integrals();  // handled by HF class
    form_H();
    
    if (ri_integrals_ == false && use_out_of_core_ == false && direct_integrals_ == false)
        form_PK();
    else if (ri_integrals_ == true)
        form_B();



    form_Shalf();
    form_initialF();
    // Check to see if there are MOs already in the checkpoint file.
    // If so, read them in instead of forming them.
    string prefix(chkpt_->build_keyword(const_cast<char*>("MO coefficients")));
    if (chkpt_->exist(const_cast<char*>(prefix.c_str()))) {
        fprintf(outfile, "  Reading previous MOs from file32.\n\n");

        // Read MOs from checkpoint and set C_ to them
        double **vectors = chkpt_->rd_scf();
        C_->set(const_cast<const double**>(vectors));
        free_block(vectors);

        form_D();

        // Read SCF energy from checkpoint file.
        E_ = chkpt_->rd_escf();
    } else {
        form_C();
        form_D();
        // Compute an initial energy using H and D
        E_ = compute_initial_E();
    }

    fprintf(outfile, "                                  Total Energy            Delta E              Density RMS\n\n");
    // SCF iterations
    do {
        iteration++;

        Dold_->copy(D_);  // Save previous density
        Eold_ = E_;       // Save previous energy

        form_J(); //J is always needed, and sometimes you get 
                  //K for free-ish with that
        J_->print(outfile);
        if (true||functional_->hasExactExchange()) {
            form_K(); //Sometimes, you really need K too
            K_->print(outfile);
        }

        form_V(); //Whoa, there's that Kohn-Sham stuff

        form_F();
        F_->print(outfile);

        if (diis_enabled_)
            save_fock();

        E_ = compute_E();

        if (diis_enabled_ == true && iteration >= num_diis_vectors_) {
            diis();
            diis_iter = true;
        } else {
            diis_iter = false;
        }

        fprintf(outfile, "  @RKS iteration %3d energy: %20.14f    %20.14f %20.14f %s\n", iteration, E_, E_ - Eold_, Drms_, diis_iter == false ? " " : "DIIS");
        fflush(outfile);

        form_C();
        form_D();

        converged = test_convergency();
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
        fprintf(outfile,"  Saving Cartesian Grid\n");
        save_RHF_grid(options_, basisset_, D_, C_);
    }
    if (ri_integrals_)
    {
        free_B();
    }
    // Compute the final dipole.
    compute_multipole();

        //fprintf(outfile,"\nComputation Completed\n");
        //fflush(outfile);
    return E_;
}

void RKS::form_J()
{
    if (ri_integrals_ == false && use_out_of_core_ == false && direct_integrals_ == false) {
        //form_J_from_PK();
    }
    else if (ri_integrals_ == false && direct_integrals_ == true) {
        form_J_and_K_from_direct_integrals();
    }
    else if (ri_integrals_ == true) {
        form_J_from_RI();
    }
    else {
        //form_J_and_K_disk();
    }
}
void RKS::form_K()
{
    if (ri_integrals_ == false && use_out_of_core_ == false && direct_integrals_ == false) {
        //form_K_from_PK();
    }
    else if (ri_integrals_ == false && direct_integrals_ == true) {
        //Already formed 
    }
    else if (ri_integrals_ == true) {
        form_K_from_RI();
    }
    else {
        //Already formed
    }
}
void RKS::form_F()
{
	F_->copy(H_);
	J_->scale(2.0*functional_->getExactCoulombCoefficient());
	F_->add(J_);
	if (functional_->hasExactExchange()) {
		K_->scale(-functional_->getExactExchangeCoefficient());
		F_->add(K_);
	}
	V_->scale(0.5);
	F_->add(V_);
	#ifdef _DEBUG
	if (debug_){
		F_->print(outfile);
	}
	#endif
}
void RKS::form_V()
{
	V_->zero();
	densityCheck_ = 0.0;
	dipoleCheckX_ = 0.0;
	dipoleCheckY_ = 0.0;
	dipoleCheckZ_ = 0.0;
	IntegrationPoint q;
	double fun,val;
	int nirreps = V_->nirreps();
	int* opi = V_->rowspi();
        double fudge = 1.0/(1.0*M_PI);
	double check;
	for (integrator_->reset(); !integrator_->isDone(); ) {
		q = integrator_->getNextPoint();
		properties_->computeProperties(q.point,D_,C_);
		const double *basis_points = properties_->getPoints();
		fun = 2.0*functional_->getValue(properties_);
		//fprintf(outfile,"  Point: <%14.10f,%14.10f,%14.10f>, w = %14.10f, functional = %14.10f\n",q.point[0],q.point[1],q.point[2],q.weight,fun);
		int h_offset = 0;
		for (int h = 0; h<nirreps; h_offset+=opi[h],h++) {
			for (int i = 0; i<opi[h];i++) {
				for (int j = 0; j<=i; j++) {
				 //fprintf(outfile,"   (%4d,%4d), phi_a = %14.10f, phi_b = %14.10f\n",i,j,basis_points[i],basis_points[j]);
					val = fudge*q.weight*basis_points[i+h_offset]*fun*basis_points[j+h_offset];
					V_->add(h,i,j,val);
					if (i!=j)
						V_->add(h,j,i,val);
				}
			}
		}
		//Check Numerical Density (Gives an idea of grid error)
		check = q.weight*properties_->getDensity();
		densityCheck_+=2.0*check;
		dipoleCheckX_+=-2.0*check*q.point[0];
		dipoleCheckY_+=-2.0*check*q.point[1];
		dipoleCheckZ_+=-2.0*check*q.point[2];
	}
	V_->print(outfile);
	
}
void RKS::save_DFT_grid()
{
 FILE * grid_file = fopen((options_.get_str("NUMERICAL_GRID_FILENAME")).c_str(),"w");
 shared_ptr<Molecule> mol = basisset_->molecule();
	fprintf(grid_file,"%d\n",mol->natom());
	fprintf(grid_file,"x,y,z,Z\n");
	for (int i=0; i<mol->natom(); i++)
		fprintf(grid_file,"%14.10f,%14.10f,%14.10f,%d\n",mol->x(i),mol->y(i), mol->z(i), mol->Z(i));

	IntegrationPoint q;	
	for (integrator_->reset(); !integrator_->isDone(); )
	{
		q = integrator_->getNextPoint();
		Vector3 v = q.point;
		double w = q.weight;
		fprintf(grid_file,"%14.10f, %14.10f, %14.10f, %14.10f\n",v[0],v[1], v[2], w);
	}


	fclose(grid_file);
}
}}
