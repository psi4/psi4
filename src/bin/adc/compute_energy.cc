#include "psi4-dec.h"
#include <libtrans/integraltransform.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>
#include <libqt/qt.h>
#include "physconst.h"
#include "adc.h"

namespace psi{ namespace adc {
    
struct lambda{
    double value;
    int dpdstate;
};

double 
ADC::compute_energy()
{   
    int nprint;
    char lbl[32];
    struct lambda lmax;
    double corr_energy, trace, *tracepi, **D, *pop, **C, **Dso;
    double *omega, omega_o, omega_diff, theta;
    dpdfile2 B, V; 
    
    omega_guess_ = SharedVector(new Vector(nirrep_, rpi_));

    if(options_.get_str("REFERENCE") == "RHF"){
        corr_energy = rhf_init_tensors();
        rhf_prepare_tensors();
    }
 
    // Now we got all ingredients to set up Newton-Raphson and block-Davidson procedures.
    // The secular equation to solve is written as A^{eff}_{SS}(\omega)V_S = V_S\Omega,
    // where subscript S stands for the singly excited manifold.
    
    if(!options_.get_bool("PR"))
        fprintf(outfile, "\t==> ADC(2) Computation <==\n\n");
    else
        fprintf(outfile, "\t==> PR-ADC(2) Computation <==\n\n");

    bool first;
    int iter = 0;
    double denom;
    std::string state_top = "ADC state ";
    char **irrep_      = Process::environment.molecule()->irrep_labels();
    
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    psio_->open(PSIF_ADC_SEM,      PSIO_OPEN_OLD);
    psio_->open(PSIF_ADC,          PSIO_OPEN_OLD);

    for(int irrep = 0;irrep < nirrep_;irrep++){
        if(rpi_[irrep]){
            omega = init_array(rpi_[irrep]);
            for(int root = 0;root < rpi_[irrep];root++){
                omega_o = omega_guess_->get(irrep, root);
                first = true;
                std::ostringstream oss;
                
                for(int iter = 1;iter <= pole_max_;iter++){
                    rhf_diagonalize(irrep, root+1, first, omega_o, omega);
                    first = false;
                    denom = 1 - rhf_differentiate_omega(irrep, root);
                    //denom = 1;
                    omega_diff = (omega_o-omega[root]) / denom;
                    if(DEBUG_)  printf("%e, %10.7f\n", omega_diff, 1/denom);
                    if(fabs(omega_diff) < pow(10.0,-conv_)){
                        if(DEBUG_){
                            printf("\tpole(%d)[%d] in %d iteration: %10.7lf\n", root, irrep, iter, omega[root]);
                            printf("\tpseudo-perturbative value: %10.7lf\n", poles_[irrep][root].ps_value);
                        }
                        
                        poles_[irrep][root].iter          = iter; 
                        poles_[irrep][root].iter_value    = omega[root];
                        poles_[irrep][root].renorm_factor = 1/denom;
                        
                        sprintf(lbl, "V^(%d)_[%d]12", root, irrep);
                        dpd_file2_init(&V, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                        fprintf(outfile, "->\t%d%3s state   : %10.7f (a.u.), %10.7f (eV)\n", root+1, irrep_[irrep], omega[root], omega[root]*_hartree2ev);
                        fprintf(outfile, "\tNon-iterative: %10.7f (a.u.), %10.7f (eV)\n", poles_[irrep][root].ps_value, poles_[irrep][root].ps_value*_hartree2ev);
                        fprintf(outfile, "\t         Occ Vir        Coefficient\n");
                        fprintf(outfile, "\t---------------------------------------------\n");
                        int nprint;
                        if(nxspi_[irrep] < num_amps_) nprint = nxspi_[irrep];
                        else nprint = num_amps_;
                        amps_write(&V, nprint, outfile);
                        fprintf(outfile, "\n");
                        fprintf(outfile, "\tConverged in %3d iteration.\n", iter);
                        fprintf(outfile, "\tSquared norm of the S component: %10.7f\n", poles_[irrep][root].renorm_factor);
                        
                        sprintf(lbl, "B^(%d)_[%d]12", root, irrep);
                        dpd_file2_init(&B, PSIF_ADC, irrep, ID('O'), ID('V'), lbl);
                        theta = acos(dpd_file2_dot(&B, &V)) * 180.0 / _pi;
                        if((180.0-fabs(theta)) < theta) theta = 180.0 - fabs(theta); 
                        dpd_file2_close(&B);
                        fprintf(outfile, "\tThe S vector is rotated up to %6.3f (deg.)\n", theta);
                        if(theta > ANGL_TOL_)
                            fprintf(outfile, "\t#WARNING: Strongly rotated from the CIS state!\n");
                        fprintf(outfile, "\n");

                        // Detachment / Attachment analysis Reference: JPC 99 (1995) 14261
                        // Extistence of D vector is not considered here, just akin to CIS(D_inf) manner.
                        dpd_file2_mat_init(&V);
                        dpd_file2_mat_rd(&V);
/*                        
                        trace = 0;
                        tracepi = init_array(nirrep_);
                        fprintf(outfile, "\t==> Detachment Density Analysis <==\n");
                        for(int sub_irrep = 0;sub_irrep < nirrep_;sub_irrep++){
                            int row = V.params->rowtot[sub_irrep];
                            int col = V.params->coltot[sub_irrep^irrep];
                            if(row && col){
                                D = block_matrix(row, row);
                                pop = init_array(row);
                                int lwork = 3 * row;
                                double *work = init_array(lwork);
                                // Form D_i^j = -\sum_{a} V_i^aV_j^a
                                C_DGEMM('n', 't', row, row, col,  -1.0, &(V.matrix[sub_irrep][0][0]), col, &(V.matrix[sub_irrep][0][0]), col, 0.0, &(D[0][0]), row);
                                C_DSYEV('v', 'u', row, &(D[0][0]), row, pop, work, lwork);
                                for(int i = 0;i < row;i++){
                                    if(fabs(pop[i]) > CUTOFF_DENS_){
                                        lmax.value = D[0][i];
                                        lmax.dpdstate = 0;
                                        for(int j = 0;j < row;j++){
                                            if(fabs(D[j][i]) > fabs(lmax.value)){
                                                lmax.value = D[j][i];
                                                lmax.dpdstate = j;
                                            }
                                        }
                                       fprintf(outfile, "\tIrrep %3s, State %3d, Value = %10.7f \n", irrep_[sub_irrep], lmax.dpdstate, pop[i]);
                                    }
                                    tracepi[sub_irrep] += pop[i];
                                }
                                if(DEBUG_){
                                    printf("\nDetachment population for %s...\n", irrep_[sub_irrep]);
                                    for(int i = 0;i < row;i++) printf("%10.7f\n", pop[i]);
                                    printf("\n");
                                }
                                free_block(D);
                                free(pop);
                                free(work);
                            }
                            trace += tracepi[sub_irrep];                               
                        }
                        fprintf(outfile, "\n");
                        fprintf(outfile, "\tContribution from each irrep... \n");
                        fprintf(outfile, "\t[");
                        for(int i = 0;i < nirrep_;i++) fprintf(outfile, " %6.4e ", tracepi[i]);
                        fprintf(outfile, "]\n");
                        fprintf(outfile, "\tTrace of zeroth-order detachment density-matrix: %7.5e\n", trace);
                        fprintf(outfile, "\n");
                        free(tracepi);

                        trace = 0;
                        tracepi = init_array(nirrep_);
                        fprintf(outfile, "\t==> Attachment Density Analysis <==\n");
                        for(int sub_irrep = 0;sub_irrep < nirrep_;sub_irrep++){
                            int row = V.params->coltot[sub_irrep];
                            int col = V.params->rowtot[sub_irrep^irrep];
                            if(row && col){
                                D   = block_matrix(row, row);               // Detachment density-matrix in MO basis
                                Dso = block_matrix(nsopi_[sub_irrep], row); // Detachment density-matrix in SO basis
                                C   = Ca_->pointer(sub_irrep);
                                pop = init_array(row);
                                int lwork = 3 * row;
                                double *work = init_array(lwork);
                                // Form D_a^b = \sum_{i} V_a^iV_i^b
                                C_DGEMM('t', 'n', row, row, col,  1.0, &(V.matrix[sub_irrep^irrep][0][0]), row, &(V.matrix[sub_irrep^irrep][0][0]), row, 0.0, &(D[0][0]), row);
                                C_DSYEV('v', 'u', row, &(D[0][0]), row, pop, work, lwork);
//                                C_DGEMM('n', 'n', nsopi_[sub_irrep], row, row, 1.0, &(C[0][0]), row, &(D[0][0]), row, 0.0, &(Dso[0][0]), row);
                                for(int i = 0;i < row;i++){
                                    if(fabs(pop[i]) > CUTOFF_DENS_){
                                        lmax.value = D[0][i];
                                        lmax.dpdstate = 0;
                                        for(int j = 0;j < row;j++){
                                            if(fabs(D[j][i]) > fabs(lmax.value)){
                                                lmax.value = D[j][i];
                                                lmax.dpdstate = j;
                                            }
                                        }
                                        fprintf(outfile, "\tIrrep %3s, State %3d, Value = %10.7f \n", irrep_[sub_irrep], lmax.dpdstate, pop[i]);
//                                        for(int j = 0;j < nsopi_[sub_irrep];j++){
//                                            fprintf(outfile, "\t%10.7f\n", Dso[j][i]);
//                                        }
//                                        fprintf(outfile, "\n");
                                    }
                                    tracepi[sub_irrep] += pop[i];
                                }
                                if(DEBUG_){
                                    printf("\nAttachment population for %s...\n", irrep_[sub_irrep]);
                                    for(int i = 0;i < row;i++) printf("%10.7f\n", pop[i]);
                                    printf("\n");
                                }
                                free_block(D);
                                free_block(Dso);
                                free(pop);
                                free(work);
                            }
                            trace += tracepi[sub_irrep];
                        }
                        fprintf(outfile, "\n");
                        fprintf(outfile, "\tContribution from each irrep... \n");
                        fprintf(outfile, "\t[");
                        for(int i = 0;i < nirrep_;i++) fprintf(outfile, " %6.4e ", tracepi[i]);
                        fprintf(outfile, "]\n");
                        fprintf(outfile, "\tTrace of zeroth-order attachment density-matrix: %7.5e\n", trace);
                        fprintf(outfile, "\n");
                        free(tracepi);
*/                                                
                        fflush(outfile);
                        
                        oss << state_top << root << irrep_[irrep];
                        std::string state = oss.str();
                        Process::environment.globals[state] = omega[root];
                
                        dpd_file2_close(&V);
                        
                        break;
                    }
                    else
                        omega_o -= omega_diff;
                }
            
            }
            free(omega);
        }
    }
    
    psio_->close(PSIF_ADC, 1);
    psio_->close(PSIF_ADC_SEM, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);
/*    
    fprintf(outfile, "\t------------------------------------------------------------------------------\n");
    fprintf(outfile, "\tD/A diagnosic is based on the reference:\n");
    fprintf(outfile, "\tM. Head-Gordon, A. M. Grana, D. Maurice and C. A. White, JPC 99 (1995) 14261.\n");
    fprintf(outfile, "\tN.B. Existence of D component is *NOT* considered.\n");
    fprintf(outfile, "\t------------------------------------------------------------------------------\n\n");
*/    
    energy_ += corr_energy;
    Process::environment.globals["CURRENT ENRGY"] = energy_;
    fprintf(outfile, "->\tCorreponsing GS total energy (a.u.) = %20.14f\n", energy_);

    release_mem();
    
    return energy_;
}
    
}} // End Namespaces
