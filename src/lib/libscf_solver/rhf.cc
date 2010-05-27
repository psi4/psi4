/*
 *  rhf.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/10/08.
 *
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>
#include <libmints/properties.h>
#include "rhf.h"
#include <psi4-dec.h>

#define TIME_SCF
#define _DEBUG

using namespace psi;
using namespace std;

namespace psi { namespace scf {

RHF::RHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : HF(options, psio, chkpt)
{
    common_init();
}

RHF::~RHF()
{
    if (pk_)
        delete[] pk_;
}

void RHF::common_init()
{
    // Defaults: attempt in core algorithm
    use_out_of_core_ = false;
    Drms_ = 0.0;

    // Allocate matrix memory
    F_    = SharedMatrix(factory_.create_matrix("F"));
//    C_    = SharedMatrix(factory_.create_matrix("MO coefficients"));
    D_    = SharedMatrix(factory_.create_matrix("D"));
    Dold_ = SharedMatrix(factory_.create_matrix("D old"));
    G_    = SharedMatrix(factory_.create_matrix("G"));
    J_    = SharedMatrix(factory_.create_matrix("J"));
    K_    = SharedMatrix(factory_.create_matrix("K"));

    if (local_K_) {
        Lref_    = SharedMatrix(factory_.create_matrix("Lref"));
        L_    = SharedMatrix(factory_.create_matrix("L"));
        I_ = block_matrix(basisset_->molecule()->natom(),doccpi_[0]);
    }

    int nao = chkpt_->rd_nao();
    chkpt_->wt_nmo(nao);

    // PK super matrix for fast G
    pk_ = NULL;

    // Check out of core
    use_out_of_core_ = options_.get_bool("OUT_OF_CORE");
    if (use_out_of_core_)
        direct_integrals_ = false;

    // Print DIIS status
    fprintf(outfile, "  DIIS %s.\n", diis_enabled_ ? "enabled" : "disabled");
    fprintf(outfile, "  Out of core %s.\n", use_out_of_core_ ? "enabled" : "disabled");
    fprintf(outfile, "  Direct %s.\n", direct_integrals_ ? "enabled": "disabled");
    fprintf(outfile, "  Density Fitting %s.\n", ri_integrals_ ? "enabled": "disabled");
    fprintf(outfile, "  Local Exchange %s.\n", local_K_ ? "enabled": "disabled");

    fflush(outfile);
    // Allocate memory for PK matrix
    if (direct_integrals_ == false && ri_integrals_ == false)
        allocate_PK();
}

double RHF::compute_energy()
{
    bool converged = false, diis_iter = false;
    iteration_ = 0;

    // Do the initial work to get the iterations started.
    //form_multipole_integrals();  // handled by HF class
    timer_on("Core Hamiltonian");
    form_H(); //Core Hamiltonian
    timer_off("Core Hamiltonian");
    timer_on("Overlap Matrix");
    form_Shalf(); //Shalf Matrix
    timer_off("Overlap Matrix");
    //Form initial MO guess by user specified method
    // Check to see if there are MOs already in the checkpoint file.
    // If so, read them in instead of forming them, unless the user disagrees.
    timer_on("Initial Guess");
    load_or_compute_initial_C();
    timer_off("Initial Guess");

    //S_->print(outfile);

    if (ri_integrals_ == false && use_out_of_core_ == false && direct_integrals_ == false)
        form_PK();
    else if (ri_integrals_ == true && local_K_ == false)
        form_B();
    else if (ri_integrals_ == true && local_K_ == true)   
        form_A();

    //if (local_K_)
    //    I_ = block_matrix(basisset_->molecule()->natom(),doccpi_[0]);

    fprintf(outfile, "                                  Total Energy            Delta E              Density RMS\n\n");
    // SCF iterations
    do {
        iteration_++;

        Dold_->copy(D_);  // Save previous density
        Eold_ = E_;       // Save previous energy

        //form_G_from_J_and_K(1.0);
        //D_->print(outfile);

        timer_on("Form G");
        if (ri_integrals_ == false && use_out_of_core_ == false && direct_integrals_ == false)
            form_G_from_PK();
        else if (ri_integrals_ == false && direct_integrals_ == true)
            form_G_from_direct_integrals();
        else if (ri_integrals_ == true && local_K_ == false)
           form_G_from_RI();
        else if (ri_integrals_ == true && local_K_ == true)
           form_G_from_RI_local_K();
        else
           form_G();
        timer_off("Form G");

        J_->print(outfile);
      
        K_->print(outfile);
        
        form_F();
        //F_->print(outfile);
        
        if (diis_enabled_)
            save_fock();

        E_ = compute_E();

        timer_on("DIIS");
        if (diis_enabled_ == true && iteration_ >= num_diis_vectors_) {
            diis();
            diis_iter = true;
        } else {
            diis_iter = false;
        }
        timer_off("DIIS");

        fprintf(outfile, "  @RHF iteration %3d energy: %20.14f    %20.14f %20.14f %s\n", iteration_, E_, E_ - Eold_, Drms_, diis_iter == false ? " " : "DIIS");
        fflush(outfile);

        timer_on("Diagonalize H");
        form_C();
        timer_off("Diagonalize H");
        form_D();

        converged = test_convergency();
    } while (!converged && iteration_ < maxiter_);

    if (converged) {
        fprintf(outfile, "\n  Energy converged.\n");
        fprintf(outfile, "\n  @RHF Final Energy: %20.14f", E_);
        if (perturb_h_) {
            fprintf(outfile, " with %f perturbation", lambda_);
        }
        fprintf(outfile, "\n");
        save_information();
    } else {
        fprintf(outfile, "\n  Failed to converged.\n");
        E_ = 0.0;
    }
    //timer_done();
    if (save_grid_) {
        // DOWN FOR MAINTENANCE
        //fprintf(outfile,"\n  Saving Cartesian Grid\n");
        //save_RHF_grid(options_, basisset_, D_, C_);
    }

    if (ri_integrals_ && !local_K_) 
        free_B();
    else if (ri_integrals_ && local_K_)
        free_A();

    //TEST LOCALIZATION
    //fully_localize_mos(); 
    //localized_Lodwin_charges();
    //propagate_local_mos(); 

    if (local_K_)
        free_block(I_);    

    // Compute the final dipole.
    compute_multipole();

    //fprintf(outfile,"\nComputation Completed\n");
    fflush(outfile);
    return E_;
}
bool RHF::load_or_compute_initial_C()
{    
    bool ret = false;
    // Check to see if there are MOs already in the checkpoint file.
    // If so, read them in instead of forming them.
    string prefix(chkpt_->build_keyword(const_cast<char*>("MO coefficients")));
    //What does the user want?
    //Options will be:
    // ""-Either READ or CORE (try READ first)
    // "READ"-try to read MOs from checkpoint
    // "BASIS2"-try to read MOs from checkpoint after basis cast up (in INPUT)
    // "CORE"-CORE Hamiltonain
    // "GWH"-Generalized Wolfsberg-Helmholtz
    // "SAD"-Superposition of Atomic Denisties TODO
    string guess_type = options_.get_str("GUESS");
    if ((guess_type == "" || guess_type == "READ"||guess_type == "BASIS2")&&chkpt_->exist(const_cast<char*>(prefix.c_str()))) {
        //Try for existing MOs already. Deuces of loaded spades
        fprintf(outfile, "  SCF Guess: Reading previous MOs.\n\n");

        // Read MOs from checkpoint and set C_ to them
        double **vectors = chkpt_->rd_scf();
        C_->set(const_cast<const double**>(vectors));
        free_block(vectors);

        form_D();

        // Read SCF energy from checkpoint file.
        E_ = chkpt_->rd_escf();
        ret = true;
    } else if (guess_type == "CORE" || guess_type == "") {
        //CORE is an old Psi standby, so we'll play this as spades
        fprintf(outfile, "  SCF Guess: Core (One-Electron) Hamiltonian.\n\n");
        F_->copy(H_); //Try the core Hamiltonian as the Fock Matrix
        form_C();
        form_D();
        // Compute an initial energy using H and D
        E_ = compute_initial_E();
    } else if (guess_type == "SAD") {
        //Superposition of Atomic Density (will be preferred when we can figure it out)
        fprintf(outfile, "  SCF Guess: Superposition of Atomic Density.\n\n");
        throw std::domain_error("SAD is not implemented for RHF/RKS yet!");
    } else if (guess_type == "GWH") {
        //Generalized Wolfsberg Helmholtz (Sounds cool, easy to code)
        fprintf(outfile, "  SCF Guess: Generalized Wolfsberg-Helmholtz.\n\n");
        F_->zero(); //Try F_{mn} = S_{mn} (H_{mm} + H_{nn})/2
        int h, i, j;
        S_->print(outfile);
        int *opi = S_->rowspi();
        int nirreps = S_->nirreps();
        for (h=0; h<nirreps; ++h) {
            for (i=0; i<opi[h]; ++i) {
                for (j=0; j<opi[h]; ++j) {
                    F_->set(h,i,j,0.5*S_->get(h,i,j)*(H_->get(h,i,i)+H_->get(h,j,j)));
                }
            }
        }
        form_C();
        form_D();
        // Compute an initial energy using H and D
        E_ = compute_initial_E();
    } else if (guess_type == "READ" || guess_type == "BASIS2") {
        throw std::invalid_argument("Checkpoint MOs requested, but do not exist!");
    }   
    return ret; 
}
void RHF::compute_multipole()
{
    // Begin dipole
    double dex, dey, dez, dx, dy, dz;
    // Convert blocked density to a full block
    SimpleMatrix D(D_->to_simple_matrix());

    dex = D.vector_dot(Dipole_[0]) * 2.0;
    dey = D.vector_dot(Dipole_[1]) * 2.0;
    dez = D.vector_dot(Dipole_[2]) * 2.0;

    dx = dex + nuclear_dipole_contribution_[0];
    dy = dey + nuclear_dipole_contribution_[1];
    dz = dez + nuclear_dipole_contribution_[2];

    double d;
    d = sqrt(dx * dx + dy * dy + dz * dz);
    // End dipole

    // Begin quadrupole
    double qexx, qexy, qexz, qeyy, qeyz, qezz;
    double mexx, mexy, mexz, meyy, meyz, mezz;
    double texx, texy, texz, teyy, teyz, tezz;

    mexx = D.vector_dot(Quadrupole_[0]) * 2.0;
    mexy = D.vector_dot(Quadrupole_[1]) * 2.0;
    mexz = D.vector_dot(Quadrupole_[2]) * 2.0;
    meyy = D.vector_dot(Quadrupole_[3]) * 2.0;
    meyz = D.vector_dot(Quadrupole_[4]) * 2.0;
    mezz = D.vector_dot(Quadrupole_[5]) * 2.0;

    texx = mexx + nuclear_quadrupole_contribution_[0];
    texy = mexy + nuclear_quadrupole_contribution_[1];
    texz = mexz + nuclear_quadrupole_contribution_[2];
    teyy = meyy + nuclear_quadrupole_contribution_[3];
    teyz = meyz + nuclear_quadrupole_contribution_[4];
    tezz = mezz + nuclear_quadrupole_contribution_[5];

    qexx = texx - (teyy+tezz)/2.0;
    qeyy = teyy - (texx+tezz)/2.0;
    qezz = tezz - (texx+teyy)/2.0;
    qexy = 1.5 * texy;
    qexz = 1.5 * texz;
    qeyz = 1.5 * teyz;

    SimpleVector evals(3);
    SimpleMatrix evecs(3, 3), temp(3, 3);

    temp.set(0, 0, qexx);
    temp.set(1, 1, qeyy);
    temp.set(2, 2, qezz);
    temp.set(0, 1, qexy);
    temp.set(0, 2, qexz);
    temp.set(1, 2, qeyz);
    temp.set(1, 0, qexy);
    temp.set(2, 0, qexz);
    temp.set(2, 1, qeyz);
    temp.diagonalize(&evecs, &evals);
    // End Quadrupole

    fprintf(outfile, "\n  Electric dipole (a.u.):\n");
    fprintf(outfile, "\t%15s\t%15s\t%15s", "X", "Y", "Z");
    fprintf(outfile, "\n    Nuclear part:\n");
    fprintf(outfile, "\t%15.10f\t%15.10f\t%15.10f\n", nuclear_dipole_contribution_[0], nuclear_dipole_contribution_[1], nuclear_dipole_contribution_[2]);

    fprintf(outfile, "\n    Electronic part:\n");
    fprintf(outfile, "\t%15.10f\t%15.10f\t%15.10f\n", dex, dey, dez);

    fprintf(outfile, "\n    Dipole moments:\n");
    fprintf(outfile, "\t%15.10f\t%15.10f\t%15.10f\n", dx, dy, dz);

    fprintf(outfile, "\n    Total dipole: %15.10f a.u.  %15.10f Debye\n", d, d*_dipmom_au2debye);
    fprintf(outfile, "    Conversion: 1.0 a.u. = %15.10f Debye\n", _dipmom_au2debye);

    // fprintf(outfile, "\n    Orbital contributions (a.u.)\n");
    // fprintf(outfile, "\t%6s %3s%15s  %15s  %15s\n", "Irrep", "MO", "X", "Y", "Z");
    // for (int h=0; h<Dipole_[0].nirreps(); ++h) {
    //   for (int i=0; i<doccpi_[h]; ++i) {
    //     fprintf(outfile, "\t%6d %3d%15.10f  %15.10f  %15.10f\n", h+1, i+1,
    //         Dipole_[0].get(h, i, i) * -2.0, Dipole_[1].get(h, i, i) * -2.0, Dipole_[2].get(h, i, i) * -2.0);
    //   }
    // }
    // fprintf(outfile, "\n  Electric quadrupole (a.u.):");
    // fprintf(outfile, "\n    Nuclear part:\n");
    // fprintf(outfile, "\txx=%15.10f\txy=%15.10f\txz=%15.10f\n", nuclear_quadrupole_contribution_[0], nuclear_quadrupole_contribution_[1], nuclear_quadrupole_contribution_[2]);
    // fprintf(outfile, "\tyy=%15.10f\tyz=%15.10f\tzz=%15.10f\n", nuclear_quadrupole_contribution_[3], nuclear_quadrupole_contribution_[4], nuclear_quadrupole_contribution_[5]);
    // fprintf(outfile, "\n    Electronic part (cross terms do not match oeprop):\n");
    // fprintf(outfile, "\txx=%15.10f\txy=%15.10f\txz=%15.10f\n", qexx, qexy, qexz);
    // fprintf(outfile, "\tyy=%15.10f\tyz=%15.10f\tzz=%15.10f\n", qeyy, qeyz, qezz);
    // fprintf(outfile, "\n    Principal values (a.u.) and axis:\n");
    // fprintf(outfile, "\tQ1=%15.10f\tV1=(%15.10f %15.10f %15.10f)\n", evals.get(0), evecs.get(0, 0), evecs.get(1, 0), evecs.get(2, 0));
    // fprintf(outfile, "\tQ2=%15.10f\tV2=(%15.10f %15.10f %15.10f)\n", evals.get(1), evecs.get(0, 1), evecs.get(1, 1), evecs.get(2, 1));
    // fprintf(outfile, "\tQ3=%15.10f\tV3=(%15.10f %15.10f %15.10f)\n", evals.get(2), evecs.get(0, 2), evecs.get(1, 2), evecs.get(2, 2));

    // Compute orbital extents
    fprintf(outfile, "\n  Orbital extents (a.u.):\n");
    fprintf(outfile, "\t%3s%15s  %15s  %15s  %15s\n", "MO", "<x^2>", "<y^2>", "<z^2>", "<r^2>");
    SimpleMatrix C(C_->to_simple_matrix());
    for (int i=0; i<C.rows(); ++i) {
        double sumx=0.0, sumy=0.0, sumz=0.0;
        for (int k=0; k<C.cols(); ++k) {
            for (int l=0; l<C.cols(); ++l) {
                double tmp = C.get(k, i) * C.get(l, i);
                sumx += Quadrupole_[0]->get(k, l) * tmp;
                sumy += Quadrupole_[3]->get(k, l) * tmp;
                sumz += Quadrupole_[5]->get(k, l) * tmp;
            }
        }
        fprintf(outfile, "\t%3d%15.10f  %15.10f  %15.10f  %15.10f\n", i+1, fabs(sumx), fabs(sumy), fabs(sumz), fabs(sumx + sumy + sumz));
        //fflush(outfile);
    }
}
/** CURRENTLY DOWN FOR MAINTENANCE
void RHF::save_RHF_grid(Options& opts, shared_ptr<BasisSet> basis, SharedMatrix D, SharedMatrix C)
{
    SharedProperties prop = SharedProperties(Properties::constructProperties(basis));
    int nmo_c;
    int* mo_inds_c;
    if (opts.get_int("N_CARTESIAN_MOS") != 0)
    {
        nmo_c = opts.get_int("N_CARTESIAN_MOS");
        mo_inds_c = init_int_array(nmo_c);
        for (int k = 0; k<nmo_c; k++)
        {
            mo_inds_c[k] = opts["CARTESIAN_MO_INDICES"][k].to_integer();
        }
        prop->setToComputeMOs(true,mo_inds_c,nmo_c);
    }
    double * extents = getCartesianGridExtents(opts, basis->molecule());
    double xmin = extents[0];
    double xmax = extents[1];
    double ymin = extents[2];
    double ymax = extents[3];
    double zmin = extents[4];
    double zmax = extents[5];
    free(extents);

    int * npoints = getCartesianGridResolution(opts);
    int npointsx = npoints[0];
    int npointsy = npoints[1];
    int npointsz = npoints[2];
    free(npoints);
    FILE * grid_file = fopen((opts.get_str("CARTESIAN_FILENAME")).c_str(),"w");
    shared_ptr<Molecule> mol = basisset_->molecule();
    fprintf(grid_file,"%d\n",mol->natom());
    fprintf(grid_file,"x,y,z,Z\n");
    for (int i=0; i<mol->natom(); i++)
        fprintf(grid_file,"%14.10f,%14.10f,%14.10f,%d\n",mol->x(i),mol->y(i), mol->z(i), mol->Z(i));
    fprintf(grid_file,"%10d,%10d,%10d\n",npointsx,npointsy,npointsz);
    fprintf(grid_file,"x,y,z,\\rho");
    if (nmo_c>0){
        for (int m=0; m<nmo_c; m++)
            fprintf(grid_file,",MO%d",mo_inds_c[m]);
    }
    fprintf(grid_file,"\n");
    double x,y,z,val;
    for (int i = 0; i<npointsx; i++) {
        if (npointsx == 1)
            x = (xmax+xmin)/2.0;
        else
            x = i*(xmax-xmin)/(npointsx-1.0)+xmin;
        for (int j = 0; j<npointsy; j++) {
            if (npointsy == 1)
                y = (ymax+ymin)/2.0;
            else
                y = j*(ymax-ymin)/(npointsy-1.0)+ymin;
            for (int k = 0; k<npointsz; k++) {
                if (npointsz == 1)
                    z = (zmax+zmin)/2.0;
                else
                    z = k*(zmax-zmin)/(npointsz-1.0)+zmin;
                Vector3 v(x,y,z);
                prop->computeProperties(v,D,C); //Needs to be updated as props changes
                val = prop->getDensity();
                fprintf(grid_file, "%14.10f,%14.10f,%14.10f,%14.10f",x,y,z,val);
                if (nmo_c>0) {
                    for (int m = 0; m<nmo_c; m++)
                        fprintf(grid_file,",%14.10f",prop->getMO(m));
                }
                fprintf(grid_file, "\n");
            }
        }
    }
    fclose(grid_file);

    if (opts.get_int("N_CARTESIAN_MOS") != 0)
        free(mo_inds_c);

}
int* RHF::getCartesianGridResolution(Options &opts)
{
    int *npoints = init_int_array(3);
    //Use individual resolutions
    if (opts["CARTESIAN_RESOLUTION_X"].has_changed()) {
     npoints[0] = opts.get_int("CARTESIAN_RESOLUTION_X");
     npoints[1] = opts.get_int("CARTESIAN_RESOLUTION_Y");
     npoints[2] = opts.get_int("CARTESIAN_RESOLUTION_Z");
    } else { //Use global or default resolutions
     npoints[0] = opts.get_int("CARTESIAN_RESOLUTION");
     npoints[1] = npoints[0];
     npoints[2] = npoints[0];
    }
    return npoints;
}
double* RHF::getCartesianGridExtents(Options &opts, shared_ptr<Molecule> mol)
{
    double *ext = init_array(6);
    if (opts["CARTESIAN_EXTENTS"].has_changed()) {
    //Use defined extents
    for (int i=0; i<6; i++)
        ext[i] = opts["CARTESIAN_EXTENTS"][i].to_double();
    } else {
        //Use molecule size plus overage
        double xmin = mol->x(0);
        double xmax = mol->x(0);
        double ymin = mol->y(0);
        double ymax = mol->y(0);
        double zmin = mol->z(0);
        double zmax = mol->z(0);
        for (int i=1;i<mol->natom(); i++) {
            if (xmin>mol->x(i))
                xmin = mol->x(i);
            if (ymin>mol->y(i))
                ymin = mol->y(i);
            if (zmin>mol->z(i))
                zmin = mol->z(i);
            if (xmax<mol->x(i))
                xmax = mol->x(i);
            if (ymax<mol->y(i))
                ymax = mol->y(i);
            if (zmax<mol->z(i))
                zmax = mol->z(i);
        }
        double overage = opts.get_double("CARTESIAN_OVERAGE");
        ext[0] = xmin - overage;
        ext[1] = xmax + overage;
        ext[2] = ymin - overage;
        ext[3] = ymax + overage;
        ext[4] = zmin - overage;
        ext[5] = zmax + overage;
    }
    return ext;
}**/
void RHF::save_information()
{
    // Print the final docc vector
    char **temp2 = chkpt_->rd_irr_labs();
    int nso = chkpt_->rd_nso();

    fprintf(outfile, "\n  Final occupation vector = (");
    for (int h=0; h<factory_.nirreps(); ++h) {
        fprintf(outfile, "%2d %3s ", doccpi_[h], temp2[h]);
    }
    fprintf(outfile, ")\n");

    // Needed for a couple of places.
    SharedMatrix eigvector(factory_.create_matrix());
    SharedVector eigvalues(factory_.create_vector());

    F_->diagonalize(eigvector, eigvalues);

    int print_mos = false;
    print_mos = options_.get_bool("PRINT_MOS");
    if (print_mos) {
        fprintf(outfile, "\n  Molecular orbitals:\n");

        C_->eivprint(eigvalues);
    }

    // Print out orbital energies.
    std::vector<std::pair<double, int> > pairs;
    for (int h=0; h<eigvalues->nirreps(); ++h) {
        for (int i=0; i<eigvalues->dimpi()[h]; ++i)
            pairs.push_back(make_pair(eigvalues->get(h, i), h));
    }
    sort(pairs.begin(),pairs.end());
    int ndocc = 0;
    for (int i=0; i<eigvalues->nirreps(); ++i)
        ndocc += doccpi_[i];

    fprintf(outfile, "\n  Orbital energies (a.u.):\n    Doubly occupied orbitals\n      ");
    for (int i=1; i<=ndocc; ++i) {
        fprintf(outfile, "%12.6f %3s  ", pairs[i-1].first, temp2[pairs[i-1].second]);
        if (i % 4 == 0)
            fprintf(outfile, "\n      ");
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "\n    Unoccupied orbitals\n      ");
    for (int i=ndocc+1; i<=nso; ++i) {
        fprintf(outfile, "%12.6f %3s  ", pairs[i-1].first, temp2[pairs[i-1].second]);
        if ((i-ndocc) % 4 == 0)
            fprintf(outfile, "\n      ");
    }
    fprintf(outfile, "\n");

    for (int i=0; i<eigvalues->nirreps(); ++i)
        free(temp2[i]);
    free(temp2);

    int *vec = new int[eigvalues->nirreps()];
    for (int i=0; i<eigvalues->nirreps(); ++i)
        vec[i] = 0;

    chkpt_->wt_nmo(nso);
    chkpt_->wt_ref(0);        // Only RHF right now
    chkpt_->wt_etot(E_);
    chkpt_->wt_escf(E_);
    chkpt_->wt_eref(E_);
    chkpt_->wt_clsdpi(doccpi_);
    chkpt_->wt_orbspi(eigvalues->dimpi());
    chkpt_->wt_openpi(vec);
    chkpt_->wt_phase_check(0);

    // Figure out frozen core orbitals
    int nfzc = chkpt_->rd_nfzc();
    int nfzv = chkpt_->rd_nfzv();
    int *frzcpi = compute_fcpi(nfzc, eigvalues);
    int *frzvpi = compute_fvpi(nfzv, eigvalues);
    chkpt_->wt_frzcpi(frzcpi);
    chkpt_->wt_frzvpi(frzvpi);
    delete[](frzcpi);
    delete[](frzvpi);

    // Save the Fock matrix
    // Need to recompute the Fock matrix as F_ is modified during the SCF interation
    form_F();
    double *ftmp = F_->to_lower_triangle();
    chkpt_->wt_fock(ftmp);
    delete[](ftmp);

    // This code currently only handles RHF
    chkpt_->wt_iopen(0);

    // Write eigenvectors and eigenvalue to checkpoint
    double *values = eigvalues->to_block_vector();
    chkpt_->wt_evals(values);
    free(values);
    double **vectors = C_->to_block_matrix();
    chkpt_->wt_scf(vectors);
    free_block(vectors);
}

void RHF::save_fock()
{
    static bool initialized_diis_manager = false;
    if (initialized_diis_manager == false) {
        diis_manager_->set_error_vector_size(1, DIISEntry::Matrix, F_.get());
        diis_manager_->set_vector_size(1, DIISEntry::Matrix, F_.get());
        initialized_diis_manager = true;
    }

    // Determine error matrix for this Fock
    SharedMatrix FDS(factory_.create_matrix()), DS(factory_.create_matrix());
    SharedMatrix SDF(factory_.create_matrix()), DF(factory_.create_matrix());

    // FDS = F_ * D_ * S_;
    DS->gemm(false, false, 1.0, D_, S_, 0.0);
    FDS->gemm(false, false, 1.0, F_, DS, 0.0);
    // SDF = S_ * D_ * F_;
    DF->gemm(false, false, 1.0, D_, F_, 0.0);
    SDF->gemm(false, false, 1.0, S_, DF, 0.0);

    Matrix FDSmSDF;
    FDSmSDF.copy(FDS);
    FDSmSDF.subtract(SDF);

    // Orthonormalize the error matrix
    FDSmSDF.transform(Shalf_);

    diis_manager_->add_entry(2, &FDSmSDF, F_.get());
}

void RHF::diis()
{
    diis_manager_->extrapolate(1, F_.get());
}

bool RHF::test_convergency()
{
    // energy difference
    double ediff = E_ - Eold_;

    // RMS of the density
    Matrix D_rms;
    D_rms.copy(D_);
    D_rms.subtract(Dold_);
    Drms_ = sqrt(D_rms.sum_of_squares());

    if (fabs(ediff) < energy_threshold_ && Drms_ < density_threshold_)
        return true;
    else
        return false;
}

void RHF::allocate_PK()
{
    if (ri_integrals_ == true || use_out_of_core_ == true || direct_integrals_ == true)
        return;

    // The size of the pk matrix is determined in HF::form_indexing
    // Allocate memory for the PK matrix (using a vector)
    if (pk_size_ < (memory_ / sizeof(double))) {
        pk_ = new double[pk_size_];

        if (pk_ == NULL) {
            fprintf(outfile, "  Insufficient free system memory for in-core PK implementation.\n");
            fprintf(outfile, "  Switching to out-of-core algorithm.\n");
            use_out_of_core_ = true;
        } else {
             // Zero out PK
            memset(pk_, 0, pk_size_*sizeof(double));
             fprintf(outfile, "  Allocated %lu elements (%lu pairs) for PK. (%5.2f MiB)\n\n", (unsigned long)pk_size_, (unsigned long)pk_pairs_, pk_size_ * 8.0 / 1048576.0);
        }
    }
    else {
        fprintf(outfile, "  Insufficient memory for in-core PK implementation.\n");
        fprintf(outfile, "  Would need %lu elements of double memory. (%5f MiB)\n", (unsigned long)pk_size_, pk_size_ * sizeof(double) / 1048576.0);
        fprintf(outfile, "  Switching to out-of-core algorithm.\n");
        use_out_of_core_ = true;
    }
}
void RHF::form_F()
{
    F_->copy(H_);
    F_->add(G_);

#ifdef _DEBUG
    if (debug_) {
        F_->print(outfile);
    }
#endif
}
void RHF::form_C()
{
    Matrix eigvec;
    Vector eigval;
    factory_.create_matrix(eigvec);
    factory_.create_vector(eigval);

    F_->transform(Shalf_);
    F_->diagonalize(eigvec, eigval);
    C_->gemm(false, false, 1.0, Shalf_, eigvec, 0.0);

    find_occupation(eigval);

    // Save C to checkpoint file.
    //double **vectors = C_->to_block_matrix();
    //chkpt_->wt_scf(vectors);
    //free_block(vectors);

#ifdef _DEBUG
    if (debug_) {
        C_->eivprint(eigval);
    }
#endif
}

void RHF::form_D()
{
    int h, i, j, m;
    int *opi = D_->rowspi();
    int nirreps = D_->nirreps();
    double val;
    for (h=0; h<nirreps; ++h) {
        for (i=0; i<opi[h]; ++i) {
            for (j=0; j<opi[h]; ++j) {
                val = 0.0;
                for (m=0; m<doccpi_[h]; ++m)
                    val += C_->get(h, i, m) * C_->get(h, j, m);
                D_->set(h, i, j, val);
            }
        }
    }

#ifdef _DEBUG
    if (debug_) {
        D_->print(outfile);
    }
#endif
}

double RHF::compute_initial_E()
{
    double Etotal = nuclearrep_ + D_->vector_dot(H_);
    fprintf(outfile, "\n  Initial RHF energy: %20.14f\n\n", Etotal);
    fflush(outfile);
    return Etotal;
}

double RHF::compute_E()
{
    Matrix HplusF;
    HplusF.copy(H_);
    HplusF.add(F_);
    double Etotal = nuclearrep_ + D_->vector_dot(HplusF);
    return Etotal;
}

void RHF::form_PK()
{
    // struct iwlbuf ERIIN;
    int ilsti, nbuf;
    int i, j, k, l;
    int ii, jj, kk, ll;
    int is, js, ks, ls;
    int fi;
    size_t bra, ket, braket=0;
    int idx;
    int counter = 0;
    int pk_counter = 0;
    bool pk_flag = false;
    double value;

    // PK zeroed out during allocation
    fprintf(outfile, "  Forming PK matrix.\n");
    fflush(outfile);

    IWL ERIIN(psio_.get(), PSIF_SO_TEI, 0.0, 1, 1);

    do {
        ilsti = ERIIN.last_buffer();
        nbuf  = ERIIN.buffer_count();

        fi = 0;
        for (idx=0; idx<nbuf; ++idx) {
            if (ERIIN.labels()[fi] >= 0) {
                i = ERIIN.labels()[fi];
                pk_flag = false;
            }
            else {
                i = -ERIIN.labels()[fi];
                pk_flag = true;
            }
            i = ERIIN.labels()[fi] >= 0 ? ERIIN.labels()[fi] : -ERIIN.labels()[fi];
            j = ERIIN.labels()[fi+1];
            k = ERIIN.labels()[fi+2];
            l = ERIIN.labels()[fi+3];
            value = ERIIN.values()[idx];
            fi += 4;

            // Get the symmetries
            is = so2symblk_[i];
            js = so2symblk_[j];
            ks = so2symblk_[k];
            ls = so2symblk_[l];

            // Get the offset of the SO index in its symblock
            ii = so2index_[i];
            jj = so2index_[j];
            kk = so2index_[k];
            ll = so2index_[l];

            // J
            if ((is == js) && (ks == ls)) {
                bra = INDEX2(ii, jj) + pk_symoffset_[is];
                ket = INDEX2(kk, ll) + pk_symoffset_[ks];
                // _pk_symoffset corrects for the symmetry offset in the _pk vector
                braket = INDEX2(bra, ket);
                pk_[braket] += value;
                // K/2 (2nd sort)
                if ((ii != jj) && (kk != ll)) {
                    if ((is == ls) && (js == ks)) {
                        bra = INDEX2(ii, ll) + pk_symoffset_[is];
                        ket = INDEX2(jj, kk) + pk_symoffset_[js];
                        braket = INDEX2(bra, ket);
                        if ((ii == ll) || (jj == kk))
                            pk_[braket] -= 0.5 * value;
                        else
                            pk_[braket] -= 0.25 * value;
                    }
                }
            }

            // K/2 (1st sort)
            if ((is == ks) && (js == ls)) {
                bra = INDEX2(ii, kk) + pk_symoffset_[is];
                ket = INDEX2(jj, ll) + pk_symoffset_[js];
                braket = INDEX2(bra, ket);
                if ((ii == kk) || (jj == ll))
                    pk_[braket] -= 0.5 * value;
                else
                    pk_[braket] -= 0.25 * value;
            }
            pk_counter++;
            counter++;

            if (pk_flag) {
                pk_counter = 0;
            }
        }

        if (!ilsti)
            ERIIN.fetch();
    } while (!ilsti);

    // Going out of scope will close the buffer
    // iwl_buf_close(&ERIIN, 1);

    // After stage two is complete, the elements of P must be halved for the case IJ=KL.
    for (size_t ij=0; ij < pk_pairs_; ++ij)
        pk_[INDEX2(ij,ij)] *= 0.5;

    fprintf(outfile, "  Processed %d two-electron integrals. form_PK", counter);
#ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "_pk:\n");
        print_array(pk_, pk_pairs_, outfile);
    }
#endif
}

void RHF::form_G_from_PK()
{
    int nirreps = factory_.nirreps();
    int *opi = factory_.rowspi();
    size_t ij;
    double *D_vector = new double[pk_pairs_];
    double *G_vector = new double[pk_pairs_];

    G_->zero();
    memset(D_vector, 0, sizeof(double) * pk_pairs_);
    memset(G_vector, 0, sizeof(double) * pk_pairs_);

    ij=0;
    for (int h=0; h<nirreps; ++h) {
        for (int p=0; p<opi[h]; ++p) {
            for (int q=0; q<=p; ++q) {
                if (p != q)
                    D_vector[ij] = 2.0 * D_->get(h, p, q);
                else
                    D_vector[ij] = D_->get(h, p, q);
                ij++;
            }
        }
    }

#ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "PK: ij = %lu\n", (unsigned long)ij);
        fflush(outfile);
        fprintf(outfile, "PK: D matrix:\n");
        D_->print(outfile);
        fprintf(outfile, "PK: D vector (appears to be OK):\n");
        for (ij=0; ij<pk_pairs_; ++ij)
            fprintf(outfile, "PK: D vector [%lu] = %20.14f\n", (unsigned long)ij, D_vector[ij]);
    }
#endif

    double G_pq,D_pq;
    double* D_rs;
    double* G_rs;
    int pq,rs;
    double* PK_block = pk_;
    int ts_pairs = pk_pairs_;
    for(pq = 0; pq < ts_pairs; ++pq){
        G_pq = 0.0;
        D_pq = D_vector[pq];
        D_rs = &D_vector[0];
        G_rs = &G_vector[0];
        for(rs = 0; rs <= pq; ++rs){
            G_pq += *PK_block * (*D_rs);
            *G_rs += *PK_block * D_pq;

            // fprintf(outfile, "pq=%d, rs=%d, G_pq=%f, PK_block=%f, D_rs=%f, G_rs=%f, D_pq=%f\n", pq, rs, G_pq, *PK_block, *D_rs, *G_rs, D_pq);

            ++D_rs;
            ++G_rs;
            ++PK_block;
        }
        G_vector[pq] += G_pq;
    }

    // Convert G to a matrix
    ij = 0;
    for(int h = 0; h < nirreps; ++h){
        for(int p = 0; p < opi[h]; ++p){
            for(int q = 0; q <= p; ++q){
                G_->set(h,p,q, 2.0 * G_vector[ij]);
                G_->set(h,q,p, 2.0 * G_vector[ij]);
                ij++;
            }
        }
    }

#ifdef _DEBUG
    if (debug_) {
        G_->print(outfile);
    }
#endif

    delete[](D_vector);
    delete[](G_vector);
}

void RHF::form_G_from_direct_integrals()
{
    timer_on("form_G_from_direct_integrals");
    double temp1, temp2, temp3, temp4, temp5, temp6;
    int itype;

    // Zero out the G matrix
    G_->zero();

    // Need to back-transform the density from SO to AO basis
    SimpleMatrix *D = D_->to_simple_matrix();

    // D->set_name("D (AO basis) pre-transform");
    // D->print();
    // D->back_transform(basisset_->uso_to_bf());
    // D->set_name("D (AO basis) post-transform");
    // D->print();

    // Need a temporary G in the AO basis
    SimpleMatrix G;
    factory_.create_simple_matrix(G, "G (AO basis)");
    G.zero();

    // Initialize an integral object
    // Begin factor ou

    IntegralFactory integral(basisset_, basisset_, basisset_, basisset_);
    if (eri_.get() == NULL) {
        eri_ = shared_ptr<TwoBodyInt>(integral.eri());
    }
    ShellCombinationsIterator iter = integral.shells_iterator();
    const double *buffer = eri_->buffer();
    // End factor out

    fprintf(outfile, "      Computing integrals..."); fflush(outfile);
    int P, Q, R, S;
    int i, j, k, l;
    int index;
    double value;

    int count=0;
    for (iter.first(); !iter.is_done(); iter.next()) {
        P = iter.p();
        Q = iter.q();
        R = iter.r();
        S = iter.s();

        // Compute quartet
        timer_on("compute_shell");
        eri_->compute_shell(P, Q, R, S);
        timer_off("compute_shell");

        // From the quartet get all the integrals
        IntegralsIterator int_iter = iter.integrals_iterator();
        for (int_iter.first(); !int_iter.is_done(); int_iter.next()) {
            i = int_iter.i();
            j = int_iter.j();
            k = int_iter.k();
            l = int_iter.l();
            index = int_iter.index();
            value = buffer[index];

            // We only care about those greater that 1.0e-14
            if (fabs(value) > 1.0e-14) {
                itype = integral_type(i, j, k, l);
                switch(itype) {
                    case 1:
                    temp1 = D->get(i, i) * value;

                    G.add(i, i, temp1);
                    break;

                    case 2:
                    temp1 = D->get(k, k) * 2.0 * value;
                    temp2 = D->get(i, k) * value;
                    temp3 = D->get(i, i) * 2.0 * value;

                    G.add(i, i, temp1);
                    G.add(k, k, temp3);
                    G.add(i, k, -temp2);
                    G.add(k, i, -temp2);
                    break;

                    case 3:
                    temp1 = D->get(i, i) * value;
                    temp2 = D->get(i, l) * value * 2.0;

                    G.add(i, l, temp1);
                    G.add(l, i, temp1);
                    G.add(i, i, temp2);
                    break;

                    case 4:
                    temp1 = D->get(j, j) * value;
                    temp2 = D->get(i, j) * value * 2.0;

                    G.add(i, j, temp1);
                    G.add(j, i, temp1);
                    G.add(j, j, temp2);
                    break;

                    case 5:
                    temp1 = D->get(i, j) * value * 3.0;
                    G.add(i, j, temp1);
                    G.add(j, i, temp1);

                    temp2 = D->get(i, i) * value;
                    temp3 = D->get(j, j) * value;
                    G.add(j, j, -temp2);
                    G.add(i, i, -temp3);
                    break;

                    case 6:
                    temp1 = D->get(k, l) * value * 4.0;
                    temp2 = D->get(i, l) * value;
                    temp3 = D->get(i, i) * value * 2.0;
                    temp4 = D->get(i, k) * value;

                    G.add(i, i, temp1);
                    G.add(i, k, -temp2);
                    G.add(k, i, -temp2);
                    G.add(k, l, temp3);
                    G.add(l, k, temp3);
                    G.add(i, l, -temp4);
                    G.add(l, i, -temp4);
                    break;

                    case 7:
                    temp1 = D->get(i, j) * value * 4.0;
                    temp2 = D->get(j, k) * value;
                    temp3 = D->get(i, k) * value;
                    temp4 = D->get(k, k) * value * 2.0;

                    G.add(k, k,  temp1);
                    G.add(i, k, -temp2);
                    G.add(k, i, -temp2);
                    G.add(j, k, -temp3);
                    G.add(k, j, -temp3);
                    G.add(i, j,  temp4);
                    G.add(j, i,  temp4);
                    break;

                    case 8:
                    temp1 = D->get(k, k) * value * 2.0;
                    temp2 = D->get(i, j) * value * 4.0;
                    temp3 = D->get(j, k) * value;
                    temp4 = D->get(i, k) * value;

                    G.add(i, j, temp1);
                    G.add(j, i, temp1);
                    G.add(k, k, temp2);
                    G.add(i, k, -temp3);
                    G.add(k, i, -temp3);
                    G.add(j, k, -temp4);
                    G.add(k, j, -temp4);
                    break;

                    case 9:
                    temp1 = D->get(i, l) * value * 3.0;
                    temp2 = D->get(i, j) * value * 3.0;
                    temp3 = D->get(j, l) * value * 2.0;
                    temp4 = D->get(i, i) * value;

                    G.add(i, j, temp1);
                    G.add(j, i, temp1);
                    G.add(i, l, temp2);
                    G.add(l, i, temp2);
                    G.add(i, i, -temp3);
                    G.add(j, l, -temp4);
                    G.add(l, j, -temp4);
                    break;

                    case 10:
                    temp1 = D->get(j, l) * value * 3.0;
                    temp2 = D->get(i, j) * value * 3.0;
                    temp3 = D->get(j, j) * value;
                    temp4 = D->get(i, l) * value * 2.0;

                    G.add(i, j, temp1);
                    G.add(j, i, temp1);
                    G.add(j, l, temp2);
                    G.add(l, j, temp2);
                    G.add(i, l, -temp3);
                    G.add(l, i, -temp3);
                    G.add(j, j, -temp4);
                    break;

                    case 11:
                    temp1 = D->get(k, j) * value * 3.0;
                    temp2 = D->get(i, j) * value * 3.0;
                    temp3 = D->get(j, j) * value;
                    temp4 = D->get(i, k) * value * 2.0;

                    G.add(i, j, temp1);
                    G.add(j, i, temp1);
                    G.add(k, j, temp2);
                    G.add(j, k, temp2);
                    G.add(i, k, -temp3);
                    G.add(k, i, -temp3);
                    G.add(j, j, -temp4);
                    break;

                    case 12:
                    case 13:
                    case 14:
                    temp1 = D->get(k, l) * value * 4.0;
                    temp2 = D->get(i, j) * value * 4.0;
                    temp3 = D->get(j, l) * value;
                    temp4 = D->get(i, k) * value;
                    temp5 = D->get(j, k) * value;
                    temp6 = D->get(i, l) * value;

                    G.add(i, j, temp1);
                    G.add(j, i, temp1);
                    G.add(k, l, temp2);
                    G.add(l, k, temp2);
                    G.add(i, k, -temp3);
                    G.add(k, i, -temp3);
                    G.add(j, l, -temp4);
                    G.add(l, j, -temp4);
                    G.add(i, l, -temp5);
                    G.add(l, i, -temp5);
                    G.add(j, k, -temp6);
                    G.add(k, j, -temp6);
                    break;
                };

                count++;
            }
        }
    }
    fprintf(outfile, "done.  %d two-electron integrals.\n", count); fflush(outfile);
//    delete eri;

    // Set RefMatrix to RefSimpleMatrix handling symmetry blocking, if needed
    // Transform G back to symmetry blocking
    // G.transform(basisset_->uso_to_bf());
    // G.print();
    G_->set(&G);
    delete D;
    // G_->print();
    timer_off("form_G_from_direct_integrals");
}

void RHF::form_G()
{
    // struct iwlbuf ERIIN;
    int ilsti, nbuf;
    int i, j, k, l;
    int ii, jj, kk, ll;
    int is, js, ks, ls;
    int fi;
    double value=0.0;
    double temp1=0.0, temp2=0.0, temp3=0.0, temp4=0.0, temp5=0.0, temp6=0.0;
    int idx;
    int itype;
    int counter = 0;

    // Zero out the G matrix
    G_->zero();

    IWL ERIIN(psio_.get(), PSIF_SO_TEI, 0.0, 1, 1);

    do {
        ilsti = ERIIN.last_buffer();
        nbuf  = ERIIN.buffer_count();

        fi = 0;
        for (idx=0; idx<nbuf; ++idx) {
            i = ERIIN.labels()[fi] > 0 ? ERIIN.labels()[fi] : -ERIIN.labels()[fi];
            j = ERIIN.labels()[fi+1];
            k = ERIIN.labels()[fi+2];
            l = ERIIN.labels()[fi+3];
            value = ERIIN.values()[idx];
            fi += 4;

            itype = integral_type(i, j, k, l);
            switch(itype) {
                case 1:
                    ii = so2index_[i];
                    is = so2symblk_[i];
                    temp1 = D_->get(is, ii, ii) * value;

                    G_->add(is, ii, ii, temp1);

// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 1:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ii, G_->get(is, ii,ii));
//                     }
// #endif
                    break;

                    case 2:
                    ii = so2index_[i];
                    kk = so2index_[k];
                    is = so2symblk_[i];
                    ks = so2symblk_[k];

                    temp1 = D_->get(ks, kk, kk) * 2.0 * value;
                    temp2 = 0.0;
                    temp3 = D_->get(is, ii, ii) * 2.0 * value;

                    G_->add(is, ii, ii, temp1);
                    G_->add(ks, kk, kk, temp3);

                    if (is == ks) {
                        temp2 = D_->get(is, ii, kk) * value;
                        G_->add(is, ii, kk, -temp2);
                        G_->add(is, kk, ii, -temp2);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 2:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ii, G_->get(is,ii,ii));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, kk, G_->get(is,ii,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", ks, kk, kk, G_->get(ks,kk,kk));
//                     }
// #endif
                    break;

                    case 3:
                    ii = so2index_[i];
                    ll = so2index_[l];
                    is = so2symblk_[i];
                    ls = so2symblk_[l];

                    temp1 = temp2 = 0.0;
                    if (is == ls) {
                        temp1 = D_->get(is, ii, ii) * value;
                        temp2 = D_->get(is, ii, ll) * value * 2.0;

                        G_->add(is, ii, ll, temp1);
                        G_->add(is, ll, ii, temp1);
                        G_->add(is, ii, ii, temp2);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 3:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ll, G_->get(is,ii,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ii, G_->get(is,ii,ii));
//                     }
// #endif
                    break;


                    case 4:
                    ii = so2index_[i];
                    jj = so2index_[j];
                    is = so2symblk_[i];
                    js = so2symblk_[j];

                    temp1 = temp2 = 0.0;
                    if (is == js) {
                        temp1 = D_->get(js, jj, jj) * value;
                        temp2 = D_->get(is, ii, jj) * value * 2.0;

                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                        G_->add(js, jj, jj, temp2);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 4:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, jj, G_->get(js,jj,jj));
//                     }
// #endif
                    break;

                    case 5:
                    ii = so2index_[i];
                    jj = so2index_[j];
                    is = so2symblk_[i];
                    js = so2symblk_[j];

                    if (is == js) {
                        temp1 = D_->get(is, ii, jj) * value * 3.0;
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }

                    temp2 = D_->get(is, ii, ii) * value;
                    temp3 = D_->get(js, jj, jj) * value;
                    G_->add(js, jj, jj, -temp2);
                    G_->add(is, ii, ii, -temp3);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 5:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, jj, G_->get(js,jj,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ii, G_->get(is,ii,ii));
//                     }
// #endif
                    break;

                    case 6:
                    ii = so2index_[i];
                    kk = so2index_[k];
                    ll = so2index_[l];
                    is = so2symblk_[i];
                    ks = so2symblk_[k];
                    ls = so2symblk_[l];

                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (ks == ls)
                        temp1 = D_->get(ks, kk, ll) * value * 4.0;
                    if (is == ls)
                        temp2 = D_->get(is, ii, ll) * value;
                    temp3 = D_->get(is, ii, ii) * value * 2.0;
                    if (is == ks)
                        temp4 = D_->get(is, ii, kk) * value;

                    G_->add(is, ii, ii, temp1);
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp2);
                        G_->add(is, kk, ii, -temp2);
                    }
                    if (ks == ls) {
                        G_->add(ks, kk, ll, temp3);
                        G_->add(ks, ll, kk, temp3);
                    }
                    if (is == ls) {
                        G_->add(is, ii, ll, -temp4);
                        G_->add(is, ll, ii, -temp4);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 6:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ii, G_->get(is,ii,ii));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, kk, G_->get(is,ii,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", ks, kk, ll, G_->get(ks,kk,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ll, G_->get(is,ii,ll));
//                     }
// #endif
                    break;

                    case 7:
                    kk = so2index_[k];
                    ii = so2index_[i];
                    jj = so2index_[j];
                    ks = so2symblk_[k];
                    is = so2symblk_[i];
                    js = so2symblk_[j];

                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (is == js)
                        temp1 = D_->get(is, ii, jj) * value * 4.0;
                    if (js == ks)
                        temp2 = D_->get(js, jj, kk) * value;
                    if (is == ks)
                        temp3 = D_->get(is, ii, kk) * value;
                    temp4 = D_->get(ks, kk, kk) * value * 2.0;

                    G_->add(ks, kk, kk, temp1);
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp2);
                        G_->add(is, kk, ii, -temp2);
                    }
                    if (js == ks) {
                        G_->add(js, jj, kk, -temp3);
                        G_->add(js, kk, jj, -temp3);
                    }
                    if (is == js) {
                        G_->add(is, ii, jj, temp4);
                        G_->add(is, jj, ii, temp4);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 7:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", ks, kk, kk, G_->get(ks,kk,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, kk, G_->get(is,ii,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, kk, G_->get(js,jj,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                     }
// #endif
                    break;

                    case 8:
                    kk = so2index_[k];
                    ii = so2index_[i];
                    jj = so2index_[j];
                    ks = so2symblk_[k];
                    is = so2symblk_[i];
                    js = so2symblk_[j];

                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    temp1 = D_->get(ks, kk, kk) * value * 2.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 4.0;
                    if (js == ks)
                        temp3 = D_->get(js, jj, kk) * value;
                    if (is == ks)
                        temp4 = D_->get(is, ii, kk) * value;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    G_->add(ks, kk, kk, temp2);
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp3);
                        G_->add(is, kk, ii, -temp3);
                    }
                    if (js == ks) {
                        G_->add(js, jj, kk, -temp4);
                        G_->add(js, kk, jj, -temp4);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 8:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", ks, kk, kk, G_->get(ks,kk,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, kk, G_->get(is,ii,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, kk, G_->get(js,jj,kk));
//                     }
// #endif
                    break;

                    case 9:
                    ii = so2index_[i];
                    jj = so2index_[j];
                    ll = so2index_[l];
                    is = so2symblk_[i];
                    js = so2symblk_[j];
                    ls = so2symblk_[l];

                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (is == ls)
                        temp1 = D_->get(is, ii, ll) * value * 3.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 3.0;
                    if (js == ls)
                        temp3 = D_->get(js, jj, ll) * value * 2.0;
                    temp4 = D_->get(is, ii, ii) * value;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    if (is == ls) {
                        G_->add(is, ii, ll, temp2);
                        G_->add(is, ll, ii, temp2);
                    }
                    G_->add(is, ii, ii, -temp3);
                    if (js == ls) {
                        G_->add(js, jj, ll, -temp4);
                        G_->add(js, ll, jj, -temp4);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 9:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ll, G_->get(is,ii,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ii, G_->get(is,ii,ii));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, ll, G_->get(js,jj,ll));
//                     }
// #endif
                    break;

                    case 10:
                    ii = so2index_[i];
                    jj = so2index_[j];
                    ll = so2index_[l];
                    is = so2symblk_[i];
                    js = so2symblk_[j];
                    ls = so2symblk_[l];

                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (js == ls)
                        temp1 = D_->get(js, jj, ll) * value * 3.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 3.0;
                    temp3 = D_->get(js, jj, jj) * value;
                    if (is == ls)
                        temp4 = D_->get(is, ii, ll) * value * 2.0;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    if (js == ls) {
                        G_->add(js, jj, ll, temp2);
                        G_->add(js, ll, jj, temp2);
                    }
                    if (is == ls) {
                        G_->add(is, ii, ll, -temp3);
                        G_->add(is, ll, ii, -temp3);
                    }
                    G_->add(js, jj, jj, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 10:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, ll, G_->get(js,jj,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ll, G_->get(is,ii,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, jj, G_->get(js,jj,jj));
//                     }
// #endif
                    break;

                    case 11:
                    ii = so2index_[i];
                    kk = so2index_[k];
                    jj = so2index_[j];
                    is = so2symblk_[i];
                    ks = so2symblk_[k];
                    js = so2symblk_[j];

                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (ks == js)
                        temp1 = D_->get(ks, kk, jj) * value * 3.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 3.0;
                    temp3 = D_->get(js, jj, jj) * value;
                    if (is == ks)
                        temp4 = D_->get(is, ii, kk) * value * 2.0;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    if (ks == js) {
                        G_->add(ks, kk, jj, temp2);
                        G_->add(ks, jj, kk, temp2);
                    }
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp3);
                        G_->add(is, kk, ii, -temp3);
                    }
                    G_->add(js, jj, jj, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 11:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", ks, kk, jj, G_->get(ks,kk,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, kk, G_->get(is,ii,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, jj, G_->get(js,jj,jj));
//                     }
// #endif
                    break;

                    case 12:
                    case 13:
                    case 14:
                    ii = so2index_[i];
                    jj = so2index_[j];
                    kk = so2index_[k];
                    ll = so2index_[l];
                    is = so2symblk_[i];
                    js = so2symblk_[j];
                    ks = so2symblk_[k];
                    ls = so2symblk_[l];

                    temp1 = temp2 = temp3 = temp4 = temp5 = temp6 = 0.0;
                    if (ks == ls)
                        temp1 = D_->get(ks, kk, ll) * value * 4.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 4.0;
                    if (js == ls)
                        temp3 = D_->get(js, jj, ll) * value;
                    if (is == ks)
                        temp4 = D_->get(is, ii, kk) * value;
                    if (js == ks)
                        temp5 = D_->get(js, jj, kk) * value;
                    if (is == ls)
                        temp6 = D_->get(is, ii, ll) * value;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    if (ks == ls) {
                        G_->add(ks, kk, ll, temp2);
                        G_->add(ks, ll, kk, temp2);
                    }
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp3);
                        G_->add(is, kk, ii, -temp3);
                    }
                    if (js == ls) {
                        G_->add(js, jj, ll, -temp4);
                        G_->add(js, ll, jj, -temp4);
                    }
                    if (is == ls) {
                        G_->add(is, ii, ll, -temp5);
                        G_->add(is, ll, ii, -temp5);
                    }
                    if (js == ks) {
                        G_->add(js, jj, kk, -temp6);
                        G_->add(js, kk, jj, -temp6);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 12,13,14:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", ks, kk, ll, G_->get(ks,kk,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, kk, G_->get(is,ii,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, ll, G_->get(js,jj,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ll, G_->get(is,ii,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, kk, G_->get(js,jj,kk));
//                     }
// #endif
                    break;
            };
            counter++;
        }

        if (!ilsti)
            ERIIN.fetch();
    } while (!ilsti);

    // Going out of scope will close the buffer
    // iwl_buf_close(&ERIIN, 1);
    fprintf(outfile, "  Processed %6d two-electron integrals.\n", counter);
}

void RHF::form_G_from_J_and_K(double scale_K_by)
{
    form_J_and_K();
    //J_->print();
    //K_->print();
    J_->scale(2.0);
    K_->scale(-1.0);
    G_->copy(K_);
    G_->scale(scale_K_by);
    G_->add(J_);
    //G_->print();
}

void RHF::form_J_and_K()
{
    if (ri_integrals_ == false && use_out_of_core_ == false && direct_integrals_ == false) {
        //form_J_from_PK();
        //form_K_from_PK();
    }
    else if (ri_integrals_ == false && direct_integrals_ == true) {
        form_J_and_K_from_direct_integrals();
    }
    else if (ri_integrals_ == true) {
        form_J_from_RI();
        form_K_from_RI();
    }
    else {
        //form_J();
        //form_K();
    }

}

void RHF::form_J_and_K_from_direct_integrals()
{
    double temp1, temp2, temp3, temp4, temp5, temp6;
    int itype;

    // Zero out the J and K matrices
    J_->zero();
    K_->zero();


    // Need to back-transform the density from SO to AO basis
    SimpleMatrix *D = D_->to_simple_matrix();

    // D->set_name("D (AO basis) pre-transform");
    // D->print();
    // D->back_transform(basisset_->uso_to_bf());
    // D->set_name("D (AO basis) post-transform");
    // D->print();

    // Need a temporary J in the AO basis
    SimpleMatrix J;
    factory_.create_simple_matrix(J, "J (AO basis)");
    J.zero();
    // Need a temporary K in the AO basis
    SimpleMatrix K;
    factory_.create_simple_matrix(K, "K (AO basis)");
    K.zero();

    // Initialize an integral object
    // Begin factor out
    IntegralFactory integral(basisset_, basisset_, basisset_, basisset_);
    if (eri_.get() == NULL) {
        eri_ = shared_ptr<TwoBodyInt>(integral.eri());
    }
    ShellCombinationsIterator iter = integral.shells_iterator();
    const double *buffer = eri_->buffer();
    // End factor out

    //fprintf(outfile, "\n      Computing integrals..."); fflush(outfile);
    int P, Q, R, S;
    int i, j, k, l;
    int index;
    double value;

    int count=0;
    for (iter.first(); !iter.is_done(); iter.next()) {
        P = iter.p();
        Q = iter.q();
        R = iter.r();
        S = iter.s();

        // Compute quartet
        eri_->compute_shell(P, Q, R, S);

        // fprintf(outfile, "Doing shell ( %d %d | %d %d )\n", P, Q, R, S); fflush(outfile);

        // From the quartet get all the integrals
        IntegralsIterator int_iter = iter.integrals_iterator();
        for (int_iter.first(); !int_iter.is_done(); int_iter.next()) {
            i = int_iter.i();
            j = int_iter.j();
            k = int_iter.k();
            l = int_iter.l();
            index = int_iter.index();
            value = buffer[index];

            //fprintf(outfile, "\tDoing integral ( %d %d | %d %d )\n", i, j, k, l); fflush(outfile);
            //fprintf(outfile, "\n (%d, %d| %d, %d) = %20.10f", i, j, k, l, value); fflush(outfile);
            // We only care about those greater that 1.0e-14
            if (fabs(value) > 1.0e-14) {
// #ifdef _DEBUG
//                 if (debug_)

// #endif
                itype = integral_type(i, j, k, l);
                switch(itype) {
                    case 1:
                    temp1 = D->get(i, i) * value;

                    J.add(i, i, 1.0*temp1);
                    K.add(i, i, -1.0*temp1);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 1:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", i, i, G.get(i,i), D->get(i, i));
//                     }
// #endif
                    break;

                    case 2:
                    temp1 = D->get(k, k) * 1.0 * value;
                    temp2 = D->get(i, k) * value;
                    temp3 = D->get(i, i) * 1.0 * value;

                    J.add(i, i, temp1);
                    J.add(k, k, temp3);
                    K.add(i, k, -temp2);
                    K.add(k, i, -temp2);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 2:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", i, i, G.get(i,i), D->get(k, k));
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", i, k, G.get(i,k), D->get(i, k));
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", k, k, G.get(k,k), D->get(i, i));
//                     }
// #endif
                    break;

                    case 3:
                    temp1 = D->get(i, i) * value;
                    temp2 = D->get(i, l) * value * 2.0;

                    J.add(i, l, 1.0*temp1);
                    J.add(l, i, 1.0*temp1);
                    K.add(i, l, -temp1);
                    K.add(l, i, -temp1);
                    J.add(i, i, 1.0*temp2);
                    K.add(i, i, -temp2);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 3:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", i, l, G.get(i,l), D->get(i, i));
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", i, i, G.get(i,i), D->get(i, l));
//                     }
// #endif
                    break;

                    case 4:
                    temp1 = D->get(j, j) * value;
                    temp2 = D->get(i, j) * value * 2.0;

                    J.add(i, j, 1.0*temp1);
                    J.add(j, i, 1.0*temp1);
                    K.add(i, j, -temp1);
                    K.add(j, i, -temp1);
                    J.add(j, j, 1.0*temp2);
                    K.add(j, j, -temp2);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 4:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", i, j, G.get(i,j), D->get(j, j));
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", j, j, G.get(j,j), D->get(i, j));
//                     }
// #endif
                    break;

                    case 5:
                    temp1 = D->get(i, j) * value;
                    J.add(i, j, 2.0*temp1);
                    J.add(j, i, 2.0*temp1);
                    K.add(i, j, -temp1);
                    K.add(j, i, -temp1);

                    temp2 = D->get(i, i) * value;
                    temp3 = D->get(j, j) * value;
                    K.add(j, j, -temp2);
                    K.add(i, i, -temp3);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 5:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, j, G.get(j,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, i, G.get(i,i));
//                     }
// #endif
                    break;

                    case 6:
                    temp1 = D->get(k, l) * value * 2.0;
                    temp2 = D->get(i, l) * value;
                    temp3 = D->get(i, i) * value * 1.0;
                    temp4 = D->get(i, k) * value;

                    J.add(i, i, temp1);
                    K.add(i, k, -temp2);
                    K.add(k, i, -temp2);
                    J.add(k, l, temp3);
                    J.add(l, k, temp3);
                    K.add(i, l, -temp4);
                    K.add(l, i, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 6:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, i, G.get(i,i));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, k, G.get(i,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", k, l, G.get(k,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, l, G.get(i,l));
//                     }
// #endif
                    break;

                    case 7:
                    temp1 = D->get(i, j) * value * 2.0;
                    temp2 = D->get(j, k) * value;
                    temp3 = D->get(i, k) * value;
                    temp4 = D->get(k, k) * value * 1.0;

                    J.add(k, k,  temp1);
                    K.add(i, k, -temp2);
                    K.add(k, i, -temp2);
                    K.add(j, k, -temp3);
                    K.add(k, j, -temp3);
                    J.add(i, j,  temp4);
                    J.add(j, i,  temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 7:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", k, k, G.get(k,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, k, G.get(i,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, k, G.get(j,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                     }
// #endif
                    break;

                    case 8:
                    temp1 = D->get(k, k) * value * 1.0;
                    temp2 = D->get(i, j) * value * 2.0;
                    temp3 = D->get(j, k) * value;
                    temp4 = D->get(i, k) * value;

                    J.add(i, j, temp1);
                    J.add(j, i, temp1);
                    J.add(k, k, temp2);
                    K.add(i, k, -temp3);
                    K.add(k, i, -temp3);
                    K.add(j, k, -temp4);
                    K.add(k, j, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 8:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", k, k, G.get(k,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, k, G.get(i,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, k, G.get(j,k));
//                     }
// #endif
                    break;

                    case 9:
                    temp1 = D->get(i, l) * value;
                    temp2 = D->get(i, j) * value;
                    temp3 = D->get(j, l) * value * 2.0;
                    temp4 = D->get(i, i) * value;

                    J.add(i, j, 2.0*temp1);
                    J.add(j, i, 2.0*temp1);
                    J.add(i, l, 2.0*temp2);
                    J.add(l, i, 2.0*temp2);
                    K.add(i, j, -temp1);
                    K.add(j, i, -temp1);
                    K.add(i, l, -temp2);
                    K.add(l, i, -temp2);
                    K.add(i, i, -temp3);
                    K.add(j, l, -temp4);
                    K.add(l, j, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 9:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, l, G.get(i,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, i, G.get(i,i));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, l, G.get(j,l));
//                     }
// #endif
                    break;

                    case 10:
                    temp1 = D->get(j, l) * value;
                    temp2 = D->get(i, j) * value;
                    temp3 = D->get(j, j) * value;
                    temp4 = D->get(i, l) * value * 2.0;

                    J.add(i, j, 2.0*temp1);
                    J.add(j, i, 2.0*temp1);
                    J.add(j, l, 2.0*temp2);
                    J.add(l, j, 2.0*temp2);
                    K.add(i, j, -temp1);
                    K.add(j, i, -temp1);
                    K.add(j, l, -temp2);
                    K.add(l, j, -temp2);
                    K.add(i, l, -temp3);
                    K.add(l, i, -temp3);
                    K.add(j, j, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 10:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, l, G.get(j,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, l, G.get(i,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, j, G.get(j,j));
//                     }
// #endif
                    break;

                    case 11:
                    temp1 = D->get(k, j) * value;
                    temp2 = D->get(i, j) * value;
                    temp3 = D->get(j, j) * value;
                    temp4 = D->get(i, k) * value * 2.0;

                    J.add(i, j, 2.0*temp1);
                    J.add(j, i, 2.0*temp1);
                    J.add(k, j, 2.0*temp2);
                    J.add(j, k, 2.0*temp2);
                    K.add(i, j, -temp1);
                    K.add(j, i, -temp1);
                    K.add(k, j, -temp2);
                    K.add(j, k, -temp2);
                    K.add(i, k, -temp3);
                    K.add(k, i, -temp3);
                    K.add(j, j, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 11:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", k, j, G.get(k,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, k, G.get(i,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, j, G.get(j,j));
//                     }
// #endif
                    break;

                    case 12:
                    case 13:
                    case 14:
                    temp1 = D->get(k, l) * value * 2.0;
                    temp2 = D->get(i, j) * value * 2.0;
                    temp3 = D->get(j, l) * value;
                    temp4 = D->get(i, k) * value;
                    temp5 = D->get(j, k) * value;
                    temp6 = D->get(i, l) * value;

                    J.add(i, j, temp1);
                    J.add(j, i, temp1);
                    J.add(k, l, temp2);
                    J.add(l, k, temp2);
                    K.add(i, k, -temp3);
                    K.add(k, i, -temp3);
                    K.add(j, l, -temp4);
                    K.add(l, j, -temp4);
                    K.add(i, l, -temp5);
                    K.add(l, i, -temp5);
                    K.add(j, k, -temp6);
                    K.add(k, j, -temp6);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 12,13,14:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", k, l, G.get(k,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, k, G.get(i,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, l, G.get(j,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, l, G.get(i,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, k, G.get(j,k));
//                     }
// #endif
                    break;
                };

                count++;
            }
        }
    }
    //fprintf(outfile, "done. %d two-electron integrals.\n", count); fflush(outfile);
//    delete eri;

    // Set RefMatrix to RefSimpleMatrix handling symmetry blocking, if needed
    // Transform G back to symmetry blocking
    // G.transform(basisset_->uso_to_bf());
    // G.print();
    J_->set(&J);
    K_->set(&K);
    J_->scale(1.0);
    K_->scale(-1.0);
    delete D;
    // G_->print();
}

}}
