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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>

#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <libmints/mints.h>
#include <libpsi4util/libpsi4util.h>
#include <libfock/jk.h>
#include <physconst.h>
#include "libtrans/integraltransform.h"
#include "libdpd/dpd.h"
#include "uhf.h"
#include "stability.h"

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

using namespace std;
using namespace psi;
using namespace boost;

namespace psi { namespace scf {

UHF::UHF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt) : HF(options, psio, chkpt)
{
    common_init();
}

UHF::UHF(Options& options, boost::shared_ptr<PSIO> psio) : HF(options, psio)
{
    common_init();
}

UHF::~UHF()
{
}

void UHF::common_init()
{

    Drms_ = 0.0;
    // TODO: Move that to the base object
    step_scale_ = options_.get_double("FOLLOW_STEP_SCALE");
    step_increment_ = options_.get_double("FOLLOW_STEP_INCREMENT");

    Fa_     = SharedMatrix(factory_->create_matrix("F alpha"));
    Fb_     = SharedMatrix(factory_->create_matrix("F beta"));
    Da_     = SharedMatrix(factory_->create_matrix("SCF alpha density"));
    Db_     = SharedMatrix(factory_->create_matrix("SCF beta density"));
    Dt_     = SharedMatrix(factory_->create_matrix("D total"));
    Dtold_  = SharedMatrix(factory_->create_matrix("D total old"));
    Lagrangian_ = SharedMatrix(factory_->create_matrix("Lagrangian"));
    Ca_     = SharedMatrix(factory_->create_matrix("C alpha"));
    Cb_     = SharedMatrix(factory_->create_matrix("C beta"));
    Ga_     = SharedMatrix(factory_->create_matrix("G alpha"));
    Gb_     = SharedMatrix(factory_->create_matrix("G beta"));
    J_      = SharedMatrix(factory_->create_matrix("J total"));
    Ka_     = SharedMatrix(factory_->create_matrix("K alpha"));
    Kb_     = SharedMatrix(factory_->create_matrix("K beta"));

    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = SharedVector(factory_->create_vector());

}

void UHF::finalize()
{
    // Form lagrangian
    for (int h=0; h<nirrep_; ++h) {
        for (int m=0; m<Lagrangian_->rowdim(h); ++m) {
            for (int n=0; n<Lagrangian_->coldim(h); ++n) {
                double sum = 0.0;
                for (int i=0; i<doccpi_[h]; ++i) {
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i)
                        +  epsilon_b_->get(h, i) * Cb_->get(h, m, i) * Cb_->get(h, n, i);
                }
                for (int i=doccpi_[h]; i<doccpi_[h]+soccpi_[h]; ++i)
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i);

                Lagrangian_->set(h, m, n, sum);
            }
        }
    }

    Dt_.reset();
    Dtold_.reset();
    Ga_.reset();
    Gb_.reset();

    compute_nos();

    HF::finalize();
}

void UHF::save_density_and_energy()
{
    Dtold_->copy(Dt_);
    Eold_ = E_;
}

void UHF::form_G()
{

    // Push the C matrix on
    std::vector<SharedMatrix> & C = jk_->C_left();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));
    C.push_back(Cb_subset("SO", "OCC"));

    // Run the JK object
    jk_->compute();

    // Pull the J and K matrices off
    const std::vector<SharedMatrix> & J = jk_->J();
    const std::vector<SharedMatrix> & K = jk_->K();
    J_->copy(J[0]);
    J_->add(J[1]);
    Ka_ = K[0];
    Kb_ = K[1];

    Ga_->copy(J_);
    Gb_->copy(Ga_);
    Ga_->subtract(Ka_);
    Gb_->subtract(Kb_);
}

void UHF::save_information()
{
}

bool UHF::test_convergency()
{
    double ediff = E_ - Eold_;

    // Drms was computed earlier
    if (fabs(ediff) < energy_threshold_ && Drms_ < density_threshold_)
        return true;
    else
        return false;
}

void UHF::form_initialF()
{
    Fa_->copy(H_);
    Fb_->copy(H_);

    if (debug_) {
        outfile->Printf( "Initial Fock alpha matrix:\n");
        Fa_->print("outfile");
        outfile->Printf( "Initial Fock beta matrix:\n");
        Fb_->print("outfile");
    }
}

void UHF::form_F()
{
    Fa_->copy(H_);
    Fa_->add(Ga_);

    Fb_->copy(H_);
    Fb_->add(Gb_);

    if (debug_) {
        Fa_->print("outfile");
        Fb_->print("outfile");
    }
}

void UHF::form_C()
{
    diagonalize_F(Fa_, Ca_, epsilon_a_);
    diagonalize_F(Fb_, Cb_, epsilon_b_);
    if (options_.get_bool("GUESS_MIX") && (iteration_ == 0)){
        if (Ca_->nirrep() == 1){
            outfile->Printf("  Mixing alpha HOMO/LUMO orbitals (%d,%d)\n\n",nalpha_,nalpha_ + 1);
            Ca_->rotate_columns(0,nalpha_ - 1,nalpha_, pc_pi * 0.25);
            Cb_->rotate_columns(0,nbeta_ - 1,nbeta_,-pc_pi * 0.25);
        }else{
            throw InputException("Warning: cannot mix alpha HOMO/LUMO orbitals. Run in C1 symmetry.", "to 'symmetry c1'", __FILE__, __LINE__);
        }
    }
    find_occupation();
    if (debug_) {
        Ca_->print("outfile");
        Cb_->print("outfile");
    }
}

void UHF::form_D()
{
    for (int h = 0; h < nirrep_; ++h) {
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];

        if (nso == 0 || nmo == 0) continue;

        double** Ca = Ca_->pointer(h);
        double** Cb = Cb_->pointer(h);
        double** Da = Da_->pointer(h);
        double** Db = Db_->pointer(h);

        if (na == 0)
            ::memset(static_cast<void*>(Da[0]), '\0', sizeof(double)*nso*nso);
        if (nb == 0)
            ::memset(static_cast<void*>(Db[0]), '\0', sizeof(double)*nso*nso);

        C_DGEMM('N','T',nso,nso,na,1.0,Ca[0],nmo,Ca[0],nmo,0.0,Da[0],nso);
        C_DGEMM('N','T',nso,nso,nb,1.0,Cb[0],nmo,Cb[0],nmo,0.0,Db[0],nso);

    }

    Dt_->copy(Da_);
    Dt_->add(Db_);

    if (debug_) {
        outfile->Printf( "in UHF::form_D:\n");
        Da_->print();
        Db_->print();
    }
}

// TODO: Once Dt_ is refactored to D_ the only difference between this and RHF::compute_initial_E is a factor of 0.5
double UHF::compute_initial_E()
{
    Dt_->copy(Da_);
    Dt_->add(Db_);
    return nuclearrep_ + 0.5 * (Dt_->vector_dot(H_));
}

double UHF::compute_E()
{
    double one_electron_E = Dt_->vector_dot(H_);
    double two_electron_E = 0.5 * (Da_->vector_dot(Fa_) + Db_->vector_dot(Fb_) - one_electron_E);

    energies_["Nuclear"] = nuclearrep_;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = two_electron_E;
    energies_["XC"] = 0.0;
    energies_["-D"] = 0.0;

    double DH  = Dt_->vector_dot(H_);
    double DFa = Da_->vector_dot(Fa_);
    double DFb = Db_->vector_dot(Fb_);
    double Eelec = 0.5 * (DH + DFa + DFb);
    // outfile->Printf( "electronic energy = %20.14f\n", Eelec);
    double Etotal = nuclearrep_ + Eelec;
    return Etotal;
}

void UHF::compute_orbital_gradient(bool save_fock)
{
    SharedMatrix gradient_a = form_FDSmSDF(Fa_, Da_);
    SharedMatrix gradient_b = form_FDSmSDF(Fb_, Db_);
    Drms_ = 0.5*(gradient_a->rms() + gradient_b->rms());

    if(save_fock){
        if (initialized_diis_manager_ == false) {
            diis_manager_ = boost::shared_ptr<DIISManager>(new DIISManager(max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError, DIISManager::OnDisk));
            diis_manager_->set_error_vector_size(2,
                                                 DIISEntry::Matrix, gradient_a.get(),
                                                 DIISEntry::Matrix, gradient_b.get());
            diis_manager_->set_vector_size(2,
                                           DIISEntry::Matrix, Fa_.get(),
                                           DIISEntry::Matrix, Fb_.get());
            initialized_diis_manager_ = true;
        }

        diis_manager_->add_entry(4, gradient_a.get(), gradient_b.get(), Fa_.get(), Fb_.get());
    }
}

bool UHF::diis()
{
    return diis_manager_->extrapolate(2, Fa_.get(), Fb_.get());
}

bool UHF::stability_analysis_pk()
{
    // ==> Old legacy stability code <==

    outfile->Printf("WARNING: PK integrals extremely slow for stability analysis.\n");
    outfile->Printf("Proceeding, but feel free to kill the computation and use other integrals.\n");
    // Build the Fock Matrix
    SharedMatrix aMoF(new Matrix("Alpha MO basis fock matrix", nmopi_, nmopi_));
    SharedMatrix bMoF(new Matrix("Beta MO basis fock matrix", nmopi_, nmopi_));
    aMoF->transform(Fa_, Ca_);
    bMoF->transform(Fb_, Cb_);

    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);
    // Ref wfn is really "this"
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
#define ID(x) ints->DPD_ID(x)
    IntegralTransform* ints = new IntegralTransform(wfn, spaces, IntegralTransform::Unrestricted, IntegralTransform::DPDOnly,
                           IntegralTransform::QTOrder, IntegralTransform::None);
    ints->set_keep_dpd_so_ints(true);
    ints->set_keep_iwl_so_ints(true);
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
    dpd_set_default(ints->get_dpd_id());
    dpdbuf4 Aaa, Aab, Abb, I;
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
    // A_IA_jb = 2 (IA|jb)
    global_dpd_->buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "UHF Hessian (IA|jb)", 2.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    // A_IA_JB = 2 (IA|JB)
    global_dpd_->buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "UHF Hessian (IA|JB)", 2.0);
    // A_IA_JB -= (IB|JA)
    global_dpd_->buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, psrq,
                       ID("[O,V]"), ID("[O,V]"), "UHF Hessian (IA|JB)", -1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
    // A_ia_jb = 2 (ia|jb)
    global_dpd_->buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "UHF Hessian (ia|jb)", 2.0);
    // A_ia_jb -= (ib|ja)
    global_dpd_->buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, psrq,
                       ID("[o,v]"), ID("[o,v]"), "UHF Hessian (ia|jb)", -1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
    // A_IA_JB -= (IJ|AB)
    global_dpd_->buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, prqs,
                       ID("[O,V]"), ID("[O,V]"), "UHF Hessian (IA|JB)", -1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>=o]+"), ID("[v>=v]+"), 0, "MO Ints (oo|vv)");
    // A_ia_jb -= (ij|ab)
    global_dpd_->buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, prqs,
                       ID("[o,v]"), ID("[o,v]"), "UHF Hessian (ia|jb)", -1.0);
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&Aaa, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "UHF Hessian (IA|JB)");
    for(int h = 0; h < Aaa.params->nirreps; ++h){
        global_dpd_->buf4_mat_irrep_init(&Aaa, h);
        global_dpd_->buf4_mat_irrep_rd(&Aaa, h);
        for(int ia = 0; ia < Aaa.params->rowtot[h]; ++ia){
            int iabs = Aaa.params->roworb[h][ia][0];
            int aabs = Aaa.params->roworb[h][ia][1];
            int isym = Aaa.params->psym[iabs];
            int asym = Aaa.params->qsym[aabs];
            int irel = iabs - Aaa.params->poff[isym];
            int arel = aabs - Aaa.params->qoff[asym] + soccpi_[asym] + doccpi_[asym];
            for(int jb = 0; jb < Aaa.params->coltot[h]; ++jb){
                int jabs = Aaa.params->colorb[h][jb][0];
                int babs = Aaa.params->colorb[h][jb][1];
                int jsym = Aaa.params->rsym[jabs];
                int bsym = Aaa.params->ssym[babs];
                int jrel = jabs - Aaa.params->roff[jsym];
                int brel = babs - Aaa.params->soff[bsym] + soccpi_[asym] + doccpi_[asym];
                // A_IA_JB += delta_IJ F_AB - delta_AB F_IJ
                if((iabs == jabs) && (asym == bsym))
                    Aaa.matrix[h][ia][jb] += aMoF->get(asym, arel, brel);
                if((aabs == babs) && (isym == jsym))
                    Aaa.matrix[h][ia][jb] -= aMoF->get(isym, irel, jrel);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Aaa, h);
    }
    global_dpd_->buf4_close(&Aaa);

    global_dpd_->buf4_init(&Abb, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "UHF Hessian (ia|jb)");

    for(int h = 0; h < Abb.params->nirreps; ++h){
        global_dpd_->buf4_mat_irrep_init(&Abb, h);
        global_dpd_->buf4_mat_irrep_rd(&Abb, h);
        for(int ia = 0; ia < Abb.params->rowtot[h]; ++ia){
            int iabs = Abb.params->roworb[h][ia][0];
            int aabs = Abb.params->roworb[h][ia][1];
            int isym = Abb.params->psym[iabs];
            int asym = Abb.params->qsym[aabs];
            int irel = iabs - Abb.params->poff[isym];
            int arel = aabs - Abb.params->qoff[asym] + doccpi_[asym];
            for(int jb = 0; jb < Abb.params->coltot[h]; ++jb){
                int jabs = Abb.params->colorb[h][jb][0];
                int babs = Abb.params->colorb[h][jb][1];
                int jsym = Abb.params->rsym[jabs];
                int bsym = Abb.params->ssym[babs];
                int jrel = jabs - Abb.params->roff[jsym];
                int brel = babs - Abb.params->soff[bsym] + doccpi_[asym];
                // A_ia_jb += delta_ij F_ab - delta_ab F_ij
                if((iabs == jabs) && (asym == bsym))
                    Abb.matrix[h][ia][jb] += bMoF->get(asym, arel, brel);
                if((aabs == babs) && (isym == jsym))
                    Abb.matrix[h][ia][jb] -= bMoF->get(isym, irel, jrel);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&Abb, h);
    }
    global_dpd_->buf4_close(&Abb);

    /*
     *  Perform the stability analysis
     */
    std::vector<std::pair<double, int> >eval_sym;

    std::string status;
    bool redo = false;
    global_dpd_->buf4_init(&Aaa, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "UHF Hessian (IA|JB)");
    global_dpd_->buf4_init(&Aab, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "UHF Hessian (IA|jb)");
    global_dpd_->buf4_init(&Abb, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "UHF Hessian (ia|jb)");
    for(int h = 0; h < Aaa.params->nirreps; ++h) {
        int aDim = Aaa.params->rowtot[h];
        int bDim = Abb.params->rowtot[h];
        int dim =  aDim + bDim;
        if(dim == 0) continue;
        double *evals = init_array(dim);
        double **evecs = block_matrix(dim, dim);
        double **A = block_matrix(dim, dim);

        // Alpha-alpha contribution to the Hessian
        global_dpd_->buf4_mat_irrep_init(&Aaa, h);
        global_dpd_->buf4_mat_irrep_rd(&Aaa, h);
        for(int ia = 0; ia < aDim; ++ia)
            for(int jb = 0; jb < aDim; ++jb)
                A[ia][jb] = Aaa.matrix[h][ia][jb];
        global_dpd_->buf4_mat_irrep_close(&Aaa, h);

        // Alpha-beta and beta-alpha contribution to the Hessian
        global_dpd_->buf4_mat_irrep_init(&Aab, h);
        global_dpd_->buf4_mat_irrep_rd(&Aab, h);
        for(int ia = 0; ia < aDim; ++ia)
            for(int jb = 0; jb < bDim; ++jb)
                A[ia][jb + aDim] = A[jb + aDim][ia] = Aab.matrix[h][ia][jb];
        global_dpd_->buf4_mat_irrep_close(&Aab, h);
        // Beta-beta contribution to the Hessian
        global_dpd_->buf4_mat_irrep_init(&Abb, h);
        global_dpd_->buf4_mat_irrep_rd(&Abb, h);
        for(int ia = 0; ia < bDim; ++ia)
            for(int jb = 0; jb < bDim; ++jb)
                A[ia + aDim][jb + aDim] = Abb.matrix[h][ia][jb];
        global_dpd_->buf4_mat_irrep_close(&Abb, h);

        sq_rsp(dim, dim, A, evals, 1, evecs, 1e-12);

        int mindim = dim < 5 ? dim : 5;
        for(int i = 0; i < mindim; i++)
            eval_sym.push_back(std::make_pair(evals[i], h));

        // Perform totally symmetric rotations, if necessary
        if(h == 0 && options_.get_str("STABILITY_ANALYSIS") == "FOLLOW"){
            if(evals[0] < 0.0){
                redo = true;
                status = "    Negative hessian eigenvalue detected: rotating orbitals.\n";
                double scale = pc_pi*options_.get_double("FOLLOW_STEP_SCALE")/2.0;
                // Rotate the alpha orbitals
//                outfile->Printf( "OLD ORBS");
//                Ca_->print();
                for(int ia = 0; ia < Aaa.params->rowtot[h]; ++ia){
                    int iabs = Aaa.params->roworb[h][ia][0];
                    int aabs = Aaa.params->roworb[h][ia][1];
                    int isym = Aaa.params->psym[iabs];
                    int asym = Aaa.params->qsym[aabs];
                    int irel = iabs - Aaa.params->poff[isym];
                    int arel = aabs - Aaa.params->qoff[asym] + doccpi_[asym] + soccpi_[asym];

                    Ca_->rotate_columns(isym, irel, arel, scale*evecs[ia][0]);
                    outfile->Printf( "Rotating %d and %d in irrep %d by %f\n",
                            irel, arel, isym, scale*evecs[ia][0]);
                }
//                outfile->Printf( "NEW ORBS");
//                Ca_->print();
                // Rotate the beta orbitals
                for(int ia = 0; ia < Abb.params->rowtot[h]; ++ia){
                    int iabs = Abb.params->roworb[h][ia][0];
                    int aabs = Abb.params->roworb[h][ia][1];
                    int isym = Abb.params->psym[iabs];
                    int asym = Abb.params->qsym[aabs];
                    int irel = iabs - Abb.params->poff[isym];
                    int arel = aabs - Abb.params->qoff[asym] + doccpi_[asym];
                    Cb_->rotate_columns(isym, irel, arel, scale*evecs[ia+aDim][0]);
                }
            }else{
                status =  "    No totally symmetric instabilities detected: "
                          "no rotation will be performed.\n";
            }
        }

        free_block(A);
        free_block(evecs);
        delete [] evals;
    }

    outfile->Printf( "    Lowest UHF->UHF stability eigenvalues:-\n");
    print_stability_analysis(eval_sym);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    delete ints;

    outfile->Printf( "%s", status.c_str());

    return redo;

}

bool UHF::stability_analysis()
{
    if(scf_type_ != "PK"){
        boost::shared_ptr<UStab> stab = boost::shared_ptr<UStab>(new UStab());
        stab->compute_energy();
        SharedMatrix eval_sym = stab->analyze();
        outfile->Printf( "    Lowest UHF->UHF stability eigenvalues: \n");
        std::vector < std::pair < double,int > >  eval_print;
        for (int h = 0; h < eval_sym->nirrep(); ++h) {
            for (int i = 0; i < eval_sym->rowdim(h); ++i) {
                eval_print.push_back(make_pair(eval_sym->get(h,i,0),h));
            }
        }
        print_stability_analysis(eval_print);

        // And now, export the eigenvalues to a PSI4 array, mainly for testing purposes

        Process::environment.arrays["SCF STABILITY EIGENVALUES"] = eval_sym;
        if (stab->is_unstable() && options_.get_str("STABILITY_ANALYSIS") == "FOLLOW") {
            if (attempt_number_ == 1 ) {
                stab_val = stab->get_eigval();
            } else if (stab_val - stab->get_eigval() < 1e-4) {
                // We probably fell on the same minimum, increase step_scale_
                outfile->Printf("    Negative eigenvalue similar to previous one, wavefunction\n");
                outfile->Printf("    likely to be in the same minimum.\n");
                step_scale_ += step_increment_;
                outfile->Printf("    Modifying FOLLOW_STEP_SCALE to %f.\n", step_scale_);
            } else {
                stab_val = stab->get_eigval();
            }
           //     outfile->Printf( "OLD ORBS");
           //     Ca_->print();
            stab->rotate_orbs(step_scale_);
           //     outfile->Printf( "NEW ORBS");
           //     Ca_->print();

           // Ask politely SCF control for a new set of iterations
           return true;
        } else {
            outfile->Printf("    Stability analysis over.\n");
            // We are done, no more iterations
            return false;
        }

    }else{
        return stability_analysis_pk();
    }
}

void UHF::compute_nos()
{
    // Compute UHF NOs and NOONs [J. Chem. Phys. 88, 4926 (1988)] -- TDC, 8/15

    // Build S^1/2
    SharedMatrix SHalf = S_->clone();
    SHalf->power(0.5);

    // Diagonalize S^1/2 Dt S^1/2
    SharedMatrix SDS = factory_->create_shared_matrix("S^1/2 Dt S^1/2");
    SDS->copy(Da_);
    SDS->add(Db_);
    SDS->transform(SHalf);

    SharedMatrix UHF_NOs = factory_->create_shared_matrix("UHF NOs");
    SharedVector UHF_NOONs(factory_->create_vector());
    SDS->diagonalize(UHF_NOs, UHF_NOONs, descending);

    // Print the NOONs -- code ripped off from OEProp::compute_no_occupations()
    int max_num;
    if(options_.get_str("PRINT_NOONS") == "ALL") max_num = nmo_;
    else max_num = to_integer(options_.get_str("PRINT_NOONS"));

    std::vector<boost::tuple<double, int, int> > metric;
    for (int h = 0; h < UHF_NOONs->nirrep(); h++)
      for (int i = 0; i < UHF_NOONs->dimpi()[h]; i++)
        metric.push_back(boost::tuple<double,int,int>(UHF_NOONs->get(h,i), h ,i));

    std::sort(metric.begin(), metric.end(), std::greater<boost::tuple<double,int,int> >());
    int offset = nalpha_;
    int start_occ = offset - max_num;
    start_occ = (start_occ < 0 ? 0 : start_occ);
    int stop_vir = offset + max_num + 1;
    stop_vir = (int)((size_t)stop_vir >= metric.size() ? metric.size() : stop_vir);
    char** labels = basisset_->molecule()->irrep_labels();
    outfile->Printf( "\n  UHF NO Occupations:\n");
    for (int index = start_occ; index < stop_vir; index++) {
      if (index < offset) {
        outfile->Printf( "  HONO-%-2d: %4d%3s %9.7f\n", offset- index - 1,
        boost::get<2>(metric[index])+1,labels[boost::get<1>(metric[index])],
        boost::get<0>(metric[index]));
      }
      else {
        outfile->Printf( "  LUNO+%-2d: %4d%3s %9.7f\n", index - offset,
        boost::get<2>(metric[index])+1,labels[boost::get<1>(metric[index])],
        boost::get<0>(metric[index]));
      }
    }
    outfile->Printf( "\n");

    if(options_.get_bool("SAVE_UHF_NOS")){
        // Save the NOs to Ca and Cb. The resulting orbitals will be restricted.

        outfile->Printf( "  Saving the UHF Natural Orbitals.\n");

        SharedMatrix SHalf_inv = S_->clone();
        SHalf_inv->power(-0.5);
        Ca_->gemm(false, false, 1.0,SHalf_inv,UHF_NOs, 0.0);

        double actv_threshold = 0.02;

        // Transform the average Fock matrix to the NO basis
        SharedMatrix F_UHF_NOs = factory_->create_shared_matrix("Fock Matrix");
        F_UHF_NOs->copy(Fa_);
        F_UHF_NOs->add(Fb_);
        F_UHF_NOs->transform(Ca_);

        // Sort orbitals according to type (core,active,virtual) and energy
        std::vector<boost::tuple<int, double, int, int> > sorted_nos;
        for (int h = 0; h < UHF_NOONs->nirrep(); h++){
            for (int i = 0; i < UHF_NOONs->dimpi()[h]; i++){
                double noon = UHF_NOONs->get(h,i);
                int type = 0; // core      NO >= 1.98
                if (noon < actv_threshold){
                    type = 2; // virtual   NO < 0.02
                }else if (noon < 2.0 - actv_threshold){
                    type = 1; // active    0.02 <= NO < 1.98
                }
                double epsilon = F_UHF_NOs->get(h,i,i);
                sorted_nos.push_back(boost::tuple<int,double,int,int>(type,epsilon,h,i));
            }
        }
        std::sort(sorted_nos.begin(), sorted_nos.end());

        // Build the final set of UHF NOs
        std::vector<int> irrep_count(nirrep_,0);

        for (size_t i = 0; i < sorted_nos.size(); i++){
            int h = boost::get<2>(sorted_nos[i]);
            int Ca_p = boost::get<3>(sorted_nos[i]);
            int Cb_p = irrep_count[h];
            for (int mu = 0; mu < Ca_->colspi(h); mu++){
                double value = Ca_->get(h,mu,Ca_p);
                Cb_->set(h,mu,Cb_p,value);
            }
            irrep_count[h] += 1;
        }

        // Copy sorted orbitals to Ca
        Ca_->copy(Cb_);

        // Suggest an active space
        Dimension corepi(nirrep_);
        Dimension actvpi(nirrep_);
        for (size_t i = 0; i < sorted_nos.size(); i++){
            int type = boost::get<0>(sorted_nos[i]);
            int h = boost::get<2>(sorted_nos[i]);
            if (type == 0) corepi[h] += 1;
            if (type == 1) actvpi[h] += 1;
        }
        outfile->Printf("\n  Active Space from UHF-NOs (NO threshold = %.4f):\n\n",actv_threshold);

        outfile->Printf("    restricted_docc = [");
        for (int h = 0; h < nirrep_; h++){
            outfile->Printf("%s%d",h ? "," : "",corepi[h]);
        }
        outfile->Printf("]\n");

        outfile->Printf("    active = [");
        for (int h = 0; h < nirrep_; h++){
            outfile->Printf("%s%d",h ? "," : "",actvpi[h]);
        }
        outfile->Printf("]\n");
    }
}

}}
