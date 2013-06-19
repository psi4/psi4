/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*
 *  rhf.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/10/08.
 *
 */

#ifdef _OPENMP
#include <omp.h>
#endif

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
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/basisset_parser.h>
#include <libmints/mints.h>
#include <libfock/jk.h>
#include "libtrans/integraltransform.h"
#include "libdpd/dpd.h"

#include "rhf.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi { namespace scf {

RHF::RHF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt)
    : HF(options, psio, chkpt)
{
    common_init();
}

RHF::RHF(Options& options, boost::shared_ptr<PSIO> psio)
    : HF(options, psio)
{
    common_init();
}

RHF::~RHF()
{
}

void RHF::common_init()
{
    if (multiplicity_ != 1) throw PSIEXCEPTION("RHF: RHF reference is only for singlets.");
    Drms_ = 0.0;

    // Allocate matrix memory
    Fa_        = SharedMatrix(factory_->create_matrix("F"));
    Fb_        = Fa_;
    Ca_        = SharedMatrix(factory_->create_matrix("C"));
    Cb_        = Ca_;
    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = epsilon_a_;
    Da_        = SharedMatrix(factory_->create_matrix("SCF density"));
    Db_        = Da_;
    Lagrangian_ = SharedMatrix(factory_->create_matrix("X"));
    D_         = Da_;
    Dold_      = SharedMatrix(factory_->create_matrix("D old"));
    G_         = SharedMatrix(factory_->create_matrix("G"));
    J_         = SharedMatrix(factory_->create_matrix("J"));
    K_         = SharedMatrix(factory_->create_matrix("K"));

}

void RHF::finalize()
{
    // Form lagrangian
    for (int h=0; h<nirrep_; ++h) {
        for (int m=0; m<Lagrangian_->rowdim(h); ++m) {
            for (int n=0; n<Lagrangian_->coldim(h); ++n) {
                double sum = 0.0;
                for (int i=0; i<doccpi_[h]; ++i) {
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i);
                }
                Lagrangian_->set(h, m, n, sum);
            }
        }
    }

    Dold_.reset();
    G_.reset();
    J_.reset();
    K_.reset();

    HF::finalize();
}

SharedMatrix RHF::Da() const
{
    return D_;
}

void RHF::save_density_and_energy()
{
    Dold_->copy(D_);  // Save previous density
    Eold_ = E_;       // Save previous energy
}

void RHF::form_G()
{
    // Push the C matrix on
    std::vector<SharedMatrix> & C = jk_->C_left();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));
    
    // Run the JK object
    jk_->compute();

    // Pull the J and K matrices off
    const std::vector<SharedMatrix> & J = jk_->J();
    const std::vector<SharedMatrix> & K = jk_->K();
    J_ = J[0];
    J_->scale(2.0);
    K_ = K[0];

    G_->copy(J_);
    G_->subtract(K_);
}

void RHF::save_information()
{
}

void RHF::compute_orbital_gradient(bool save_fock)
{
    // Conventional DIIS (X'[FDS - SDF]X, where X levels things out)
    SharedMatrix gradient = form_FDSmSDF(Fa_, Da_);
    Drms_ = gradient->rms();

    if(save_fock){
        if (initialized_diis_manager_ == false) {
            if (scf_type_ == "direct")
                diis_manager_ = boost::shared_ptr<DIISManager>(new DIISManager(max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError, DIISManager::InCore));
            else
                diis_manager_ = boost::shared_ptr<DIISManager>(new DIISManager(max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError, DIISManager::OnDisk));
            diis_manager_->set_error_vector_size(1, DIISEntry::Matrix, gradient.get());
            diis_manager_->set_vector_size(1, DIISEntry::Matrix, Fa_.get());
            initialized_diis_manager_ = true;
        }
        diis_manager_->add_entry(2, gradient.get(), Fa_.get());
    }
}


bool RHF::diis()
{
    return diis_manager_->extrapolate(1, Fa_.get());
}

bool RHF::test_convergency()
{
    // energy difference
    double ediff = E_ - Eold_;

    // Drms was computed earlier
    if (fabs(ediff) < energy_threshold_ && Drms_ < density_threshold_)
        return true;
    else
        return false;
}


void RHF::form_F()
{
    Fa_->copy(H_);
    Fa_->add(G_);

    if (debug_) {
        Fa_->print(outfile);
        J_->print();
        K_->print();
        G_->print();
        
    }
}

void RHF::form_C()
{
    diagonalize_F(Fa_, Ca_, epsilon_a_);
    find_occupation();
}

void RHF::form_D()
{
    for (int h = 0; h < nirrep_; ++h) {
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int na = doccpi_[h];

        if (nso == 0 || nmo == 0) continue;

        double** Ca = Ca_->pointer(h);
        double** D = D_->pointer(h);

        if (na == 0)
            memset(static_cast<void*>(D[0]), '\0', sizeof(double)*nso*nso);

        C_DGEMM('N','T',nso,nso,na,1.0,Ca[0],nmo,Ca[0],nmo,0.0,D[0],nso);

    }

    if (debug_) {
        fprintf(outfile, "in RHF::form_D:\n");
        D_->print();
    }
}

void RHF::damp_update()
{
    for(int h = 0; h < nirrep_; ++h){
        for(int row = 0; row < D_->rowspi(h); ++row){
            for(int col = 0; col < D_->colspi(h); ++col){
                double Dold = damping_percentage_ * Dold_->get(h, row, col);
                double Dnew = (1.0 - damping_percentage_) * D_->get(h, row, col);
                D_->set(h, row, col, Dold+Dnew);
            }
        }
    }
}

double RHF::compute_initial_E()
{
    double Etotal = nuclearrep_ + D_->vector_dot(H_);
    return Etotal;
}

double RHF::compute_E()
{
    double one_electron_E = 2.0 * D_->vector_dot(H_);
    double two_electron_E = D_->vector_dot(Fa_) - 0.5 * one_electron_E;   
 
    energies_["Nuclear"] = nuclearrep_;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = two_electron_E; 
    energies_["XC"] = 0.0;
    energies_["-D"] = 0.0;

    Matrix HplusF;
    HplusF.copy(H_);
    HplusF.add(Fa_);

    double Etotal = nuclearrep_ + D_->vector_dot(HplusF);
    return Etotal;
}

void RHF::save_sapt_info()
{
    if (factory_->nirrep() != 1)
    {
        fprintf(outfile,"Must run in C1. Period.\n"); fflush(outfile);
        abort();
    }
    if (soccpi_[0] != 0)
    {
        fprintf(outfile,"Aren't we in RHF Here? Pair those electrons up cracker!\n"); fflush(outfile);
        abort();
    }

    fprintf(outfile,"\n  Saving SAPT %s file.\n",options_.get_str("SAPT").c_str());

    int fileno;
    char* body_type = (char*)malloc(400*sizeof(char));
    char* key_buffer = (char*)malloc(4000*sizeof(char));

    string type = options_.get_str("SAPT");
    if (type == "2-DIMER") {
        fileno = PSIF_SAPT_DIMER;
        sprintf(body_type,"Dimer");
    } else if (type == "2-MONOMER_A") {
        fileno = PSIF_SAPT_MONOMERA;
        sprintf(body_type,"Monomer");
    } else if (type == "2-MONOMER_B") {
        fileno = PSIF_SAPT_MONOMERB;
        sprintf(body_type,"Monomer");
    } else if (type == "3-TRIMER") {
        fileno = PSIF_3B_SAPT_TRIMER;
        sprintf(body_type,"Trimer");
    } else if (type == "3-DIMER_AB") {
        fileno = PSIF_3B_SAPT_DIMER_AB;
        sprintf(body_type,"Dimer");
    } else if (type == "3-DIMER_BC") {
        fileno = PSIF_3B_SAPT_DIMER_BC;
        sprintf(body_type,"Dimer");
    } else if (type == "3-DIMER_AC") {
        fileno = PSIF_3B_SAPT_DIMER_AC;
        sprintf(body_type,"Dimer");
    } else if (type == "3-MONOMER_A") {
        fileno = PSIF_3B_SAPT_MONOMER_A;
        sprintf(body_type,"Monomer");
    } else if (type == "3-MONOMER_B") {
        fileno = PSIF_3B_SAPT_MONOMER_B;
        sprintf(body_type,"Monomer");
    } else if (type == "3-MONOMER_C") {
        fileno = PSIF_3B_SAPT_MONOMER_C;
        sprintf(body_type,"Monomer");
    } else {
        throw std::domain_error("SAPT Output option invalid");
    }

    psio_->open(fileno,0);

    int sapt_nso = basisset_->nbf();
    int sapt_nmo = nmo_;
    int sapt_nocc = doccpi_[0];
    int sapt_nvir = sapt_nmo-sapt_nocc;
    int sapt_ne = 2*sapt_nocc;
    double sapt_E_HF = E_;
    double sapt_E_nuc = nuclearrep_;
    double *sapt_evals = epsilon_a_->to_block_vector();
    double **sapt_C = Ca_->to_block_matrix();
    //print_mat(sapt_C,sapt_nso,sapt_nso,outfile);
    double *sapt_V_ints = V_->to_lower_triangle();
    double *sapt_S_ints = S_->to_lower_triangle();

    int errcod;

    sprintf(key_buffer,"%s NSO",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_nso, sizeof(int));
    sprintf(key_buffer,"%s NMO",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_nmo, sizeof(int));
    sprintf(key_buffer,"%s NOCC",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_nocc, sizeof(int));
    sprintf(key_buffer,"%s NVIR",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_nvir, sizeof(int));
    sprintf(key_buffer,"%s Number of Electrons",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_ne,sizeof(int));
    sprintf(key_buffer,"%s HF Energy",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_E_HF,sizeof(double));
    sprintf(key_buffer,"%s Nuclear Repulsion Energy",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &sapt_E_nuc, sizeof(double));
    sprintf(key_buffer,"%s HF Eigenvalues",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &(sapt_evals[0]),sizeof(double)*sapt_nmo);
    sprintf(key_buffer,"%s Nuclear Attraction Integrals",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &(sapt_V_ints[0]), sizeof(double)*sapt_nso*(sapt_nso+1)/2);
    sprintf(key_buffer,"%s Overlap Integrals",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &(sapt_S_ints[0]), sizeof(double)*sapt_nso*(sapt_nso+1)/2);
    sprintf(key_buffer,"%s HF Coefficients",body_type);
    psio_->write_entry(fileno,key_buffer,(char *) &(sapt_C[0][0]),sizeof(double)*sapt_nmo*sapt_nso);

    psio_->close(fileno,1);

    delete[] sapt_evals;
    delete[] sapt_V_ints;
    delete[] sapt_S_ints;
    free_block(sapt_C);

    free(body_type);
    free(key_buffer);
}

void RHF::stability_analysis()
{
    if(scf_type_ == "DF" || scf_type_ == "CD"){
        throw PSIEXCEPTION("Stability analysis has not been implemented for density fitted wavefunctions yet.");
    }else{
#define ID(x) ints.DPD_ID(x)
        // Build the Fock Matrix
        SharedMatrix moF(new Matrix("MO basis fock matrix", nmopi_, nmopi_));
        moF->transform(Fa_, Ca_);

        std::vector<boost::shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::occ);
        spaces.push_back(MOSpace::vir);
        // Ref wfn is really "this"
        boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
        IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly,
                               IntegralTransform::QTOrder, IntegralTransform::None);
        ints.set_keep_dpd_so_ints(true);
        ints.transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
        ints.transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
        dpd_set_default(ints.get_dpd_id());
        dpdbuf4 Asing, Atrip,I;
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
        // Singlet A_ia_jb = 4 (ia|jb)
        global_dpd_->buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "RHF Singlet Hessian (IA|JB)", 4.0);
        // Triplet A_ia_jb = -(ib|ja)
        global_dpd_->buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, psrq,
                           ID("[O,V]"), ID("[O,V]"), "RHF Triplet Hessian (IA|JB)", -1.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
        // Triplet A_ia_jb -= (ij|ab)
        global_dpd_->buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, prqs,
                           ID("[O,V]"), ID("[O,V]"), "RHF Triplet Hessian (IA|JB)", -1.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_init(&Atrip, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "RHF Triplet Hessian (IA|JB)");
        for(int h = 0; h < Atrip.params->nirreps; ++h){
            global_dpd_->buf4_mat_irrep_init(&Atrip, h);
            global_dpd_->buf4_mat_irrep_rd(&Atrip, h);
            for(int ia = 0; ia < Atrip.params->rowtot[h]; ++ia){
                int iabs = Atrip.params->roworb[h][ia][0];
                int aabs = Atrip.params->roworb[h][ia][1];
                int isym = Atrip.params->psym[iabs];
                int asym = Atrip.params->qsym[aabs];
                int irel = iabs - Atrip.params->poff[isym];
                int arel = aabs - Atrip.params->qoff[asym] + doccpi_[asym];
                for(int jb = 0; jb < Atrip.params->coltot[h]; ++jb){
                    int jabs = Atrip.params->colorb[h][jb][0];
                    int babs = Atrip.params->colorb[h][jb][1];
                    int jsym = Atrip.params->rsym[jabs];
                    int bsym = Atrip.params->ssym[babs];
                    int jrel = jabs - Atrip.params->roff[jsym];
                    int brel = babs - Atrip.params->soff[bsym] + doccpi_[bsym];
                    // Triplet A_ia_jb += delta_ij F_ab - delta_ab F_ij
                    if((iabs == jabs) && (asym == bsym))
                        Atrip.matrix[h][ia][jb] += moF->get(asym, arel, brel);
                    if((aabs == babs) && (isym == jsym))
                        Atrip.matrix[h][ia][jb] -= moF->get(isym, irel, jrel);
                }
            }
            global_dpd_->buf4_mat_irrep_wrt(&Atrip, h);
        }
        // Singlet A += Triplet A
        global_dpd_->buf4_init(&Asing, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "RHF Singlet Hessian (IA|JB)");
        global_dpd_->buf4_axpy(&Atrip, &Asing, 1.0);
        global_dpd_->buf4_close(&Atrip);
        global_dpd_->buf4_close(&Asing);

        /*
         *  Perform the stability analysis
         */
        std::vector<std::pair<double, int> >singlet_eval_sym;
        std::vector<std::pair<double, int> >triplet_eval_sym;

        global_dpd_->buf4_init(&Asing, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "RHF Singlet Hessian (IA|JB)");
        global_dpd_->buf4_init(&Atrip, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "RHF Triplet Hessian (IA|JB)");
        for(int h = 0; h < Asing.params->nirreps; ++h) {
            int dim = Asing.params->rowtot[h];
            if(dim == 0) continue;
            double *evals = init_array(dim);
            double **evecs = block_matrix(dim, dim);

            global_dpd_->buf4_mat_irrep_init(&Asing, h);
            global_dpd_->buf4_mat_irrep_rd(&Asing, h);
            sq_rsp(dim, dim, Asing.matrix[h], evals, 1, evecs, 1e-12);
            global_dpd_->buf4_mat_irrep_close(&Asing, h);

            int mindim = dim < 5 ? dim : 5;
            for(int i = 0; i < mindim; i++)
                singlet_eval_sym.push_back(std::make_pair(evals[i], h));

            zero_arr(evals, dim);
            zero_mat(evecs, dim, dim);

            global_dpd_->buf4_mat_irrep_init(&Atrip, h);
            global_dpd_->buf4_mat_irrep_rd(&Atrip, h);
            sq_rsp(dim, dim, Atrip.matrix[h], evals, 1, evecs, 1e-12);
            global_dpd_->buf4_mat_irrep_close(&Atrip, h);

            for(int i = 0; i < mindim; i++)
                triplet_eval_sym.push_back(std::make_pair(evals[i], h));

            free_block(evecs);
            delete [] evals;
        }

        fprintf(outfile, "\tLowest singlet (RHF->RHF) stability eigenvalues:-\n");
        print_stability_analysis(singlet_eval_sym);
        fprintf(outfile, "\tLowest triplet (RHF->UHF) stability eigenvalues:-\n");
        print_stability_analysis(triplet_eval_sym);
        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }
}

}}
