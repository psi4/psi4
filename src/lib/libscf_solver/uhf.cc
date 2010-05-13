#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/matrix.h>
#include "uhf.h"

using namespace std;
using namespace psi;

namespace psi { namespace scf {

UHF::UHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt) : HF(options, psio, chkpt)
{
    common_init();
}

UHF::~UHF()
{
    if (p_jk_)
        delete[](p_jk_);
    if (p_k_)
        delete[](p_k_);
}

void UHF::common_init()
{
    Drms_ = 0.0;
    use_out_of_core_ = false;

    Fa_     = SharedMatrix(factory_.create_matrix("F alpha"));
    Fb_     = SharedMatrix(factory_.create_matrix("F beta"));
    pertFa_ = SharedMatrix(factory_.create_matrix("Perturbed Alpha Fock Matrix"));
    pertFb_ = SharedMatrix(factory_.create_matrix("Perturbed Beta Fock Matrix"));
    Da_     = SharedMatrix(factory_.create_matrix("D alpha"));
    Db_     = SharedMatrix(factory_.create_matrix("D beta"));
    Dt_     = SharedMatrix(factory_.create_matrix("D total"));
    Dtold_  = SharedMatrix(factory_.create_matrix("D total old"));
    Ca_     = SharedMatrix(factory_.create_matrix("C alpha"));
    Cb_     = SharedMatrix(factory_.create_matrix("C beta"));
    Ga_     = SharedMatrix(factory_.create_matrix("G alpha"));
    Gb_     = SharedMatrix(factory_.create_matrix("G beta"));

    p_jk_ = NULL;
    p_k_  = NULL;
    Va_   = NULL;
    Vb_   = NULL;

    use_out_of_core_ = options_.get_bool("OUT_OF_CORE");

//    if(print_ > 1) fprintf(outfile, "  DIIS not implemented for UHF, yet.\n\n");

    fprintf(outfile, "  DIIS %s.\n", diis_enabled_ ? "enabled" : "disabled");
    fprintf(outfile, "  Out of core %s.\n", use_out_of_core_ ? "enabled" : "disabled");
    fprintf(outfile, "  Direct %s.\n", direct_integrals_ ? "enabled": "disabled");
    fprintf(outfile, "  Density Fitting %s.\n", ri_integrals_ ? "enabled": "disabled");

    if (direct_integrals_ == false && ri_integrals_ == false)
        allocate_PK();
}

double UHF::compute_energy()
{
    bool converged = false, diis_iter = false;
    int iter = 0;

    // Do the initial work to give the iterations a starting point
    form_H();
    // find_occupation(_H, _H);

    if (ri_integrals_ == false && use_out_of_core_ == false && direct_integrals_ == false)
        form_PK();
    else if (ri_integrals_ == true)
        form_B();

    form_Shalf();
    form_initialF();

    if (load_or_compute_initial_C())
        fprintf(outfile, "  Read in previous MOs from file32.\n\n");

    if (print_)
        fprintf(outfile, "                                  Total Energy            Delta E              Density RMS\n\n");
    do {
        iter++;
        iterationsNeeded_ = iter;
        Dtold_->copy(Dt_);
        Eold_ = E_;

        if (ri_integrals_ == false && use_out_of_core_ == false && direct_integrals_ == false)
            form_G_from_PK();
        else if (ri_integrals_ == false && direct_integrals_ == true)
           form_G_from_direct_integrals();
        else if (ri_integrals_ == true)
           form_G_from_RI();
        else
           form_G();

        form_F();

        if (diis_enabled_)
            save_fock();

        E_ = compute_E();

        if (diis_enabled_ == true && iter >= num_diis_vectors_) {
            diis();
            diis_iter = true;
        } else {
            diis_iter = false;
        }

        if(print_)
            fprintf(outfile, "  @UHF iteration %3d energy: %20.14f    %20.14f %20.14f %s\n", iter, E_, E_ - Eold_, Drms_, diis_iter == false ? " " : "DIIS");
        fflush(outfile);

        form_C();
        //find_occupation(Fa_, Fb_);
        form_D();

        converged = test_convergency();
    } while (!converged && iter < maxiter_);

    if (ri_integrals_)
    {
        free_B();
    }
    // Return the final RHF energy
    if (converged) {
        fprintf(outfile, "\n  Energy converged.\n");
        fprintf(outfile, "\n  @UHF Final Energy: %20.14f\n", E_);
        save_information();
    } else {
        fprintf(outfile, "\n  Failed to converge.\n");
        E_ = 0.0;
    }

    // Compute the final dipole.
    if(print_ > 2) compute_multipole();

    return E_;
}

void UHF::compute_multipole()
{
    // Begin dipole
    double dex, dey, dez, dx, dy, dz;
    // Convert blocked density to a full block
    SimpleMatrix D(Dt_->to_simple_matrix());

    dex = D.vector_dot(Dipole_[0]);
    dey = D.vector_dot(Dipole_[1]);
    dez = D.vector_dot(Dipole_[2]);

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
    fprintf(outfile, "\n  Alpha Orbital extents (a.u.):\n");
    fprintf(outfile, "\t%3s%15s  %15s  %15s  %15s\n", "MO", "<x^2>", "<y^2>", "<z^2>", "<r^2>");
    SimpleMatrix Ca(Ca_->to_simple_matrix());
    for (int i=0; i<Ca.rows(); ++i) {
        double sumx=0.0, sumy=0.0, sumz=0.0;
        for (int k=0; k<Ca.cols(); ++k) {
            for (int l=0; l<Ca.cols(); ++l) {
                double tmp = Ca.get(k, i) * Ca.get(l, i);
                sumx += Quadrupole_[0]->get(k, l) * tmp;
                sumy += Quadrupole_[3]->get(k, l) * tmp;
                sumz += Quadrupole_[5]->get(k, l) * tmp;
            }
        }
        fprintf(outfile, "\t%3d%15.10f  %15.10f  %15.10f  %15.10f\n", i+1, fabs(sumx), fabs(sumy), fabs(sumz), fabs(sumx + sumy + sumz));
    }

    // Compute orbital extents
    fprintf(outfile, "\n  Beta Orbital extents (a.u.):\n");
    fprintf(outfile, "\t%3s%15s  %15s  %15s  %15s\n", "MO", "<x^2>", "<y^2>", "<z^2>", "<r^2>");
    SimpleMatrix Cb(Cb_->to_simple_matrix());
    for (int i=0; i<Cb.rows(); ++i) {
        double sumx=0.0, sumy=0.0, sumz=0.0;
        for (int k=0; k<Cb.cols(); ++k) {
            for (int l=0; l<Cb.cols(); ++l) {
                double tmp = Cb.get(k, i) * Cb.get(l, i);
                sumx += Quadrupole_[0]->get(k, l) * tmp;
                sumy += Quadrupole_[3]->get(k, l) * tmp;
                sumz += Quadrupole_[5]->get(k, l) * tmp;
            }
        }
        fprintf(outfile, "\t%3d%15.10f  %15.10f  %15.10f  %15.10f\n", i+1, fabs(sumx), fabs(sumy), fabs(sumz), fabs(sumx + sumy + sumz));
    }
}

bool UHF::load_or_compute_initial_C()
{
    bool ret = false;
    string alpha(chkpt_->build_keyword(const_cast<char*>("Alpha MO coefficients")));
    string beta(chkpt_->build_keyword(const_cast<char*>("Beta MO coefficients")));
    if (chkpt_->exist(const_cast<char*>(alpha.c_str())) &&
        chkpt_->exist(const_cast<char*>(beta.c_str()))) {
        // Read alpha MOs from checkpoint and set C_ to them
        double **vectors = chkpt_->rd_alpha_scf();
        Ca_->set(const_cast<const double**>(vectors));
        free_block(vectors);

        // Read beta MOs from checkpoint and set C_ to them
        vectors = chkpt_->rd_beta_scf();
        Cb_->set(const_cast<const double**>(vectors));
        free_block(vectors);

        form_D();

        // Read SCF energy from checkpoint file.
        E_ = chkpt_->rd_escf();

        ret = true;
    } else {
        form_initial_C();
        form_D();
        // Compute an initial energy using H and D
        E_ = compute_initial_E();

        ret = false;
    }

    return ret;
}

void UHF::save_information()
{
    // Print the final docc vector
    char **temp2 = chkpt_->rd_irr_labs();
    int nso = chkpt_->rd_nso();
    if(print_ > 1){
        fprintf(outfile, "\n  Final doubly occupied vector = (");
        for (int h=0; h<factory_.nirreps(); ++h) {
            fprintf(outfile, "%2d %3s ", doccpi_[h], temp2[h]);
        }
        fprintf(outfile, ")\n");
        fprintf(outfile, "  Final singly occupied vector = (");
        for (int h=0; h<factory_.nirreps(); ++h) {
            fprintf(outfile, "%2d %3s ", soccpi_[h], temp2[h]);
        }
        fprintf(outfile, ")\n");
    }

    // Needed for a couple of places.
    SharedMatrix eigvectora(factory_.create_matrix());
    SharedMatrix eigvectorb(factory_.create_matrix());
    SharedVector eigvaluesa(factory_.create_vector());
    SharedVector eigvaluesb(factory_.create_vector());
    Fa_->diagonalize(eigvectora, eigvaluesa);
    Fb_->diagonalize(eigvectorb, eigvaluesb);

    bool print_mos = options_.get_bool("PRINT_MOS");
    if (print_mos) {
        fprintf(outfile, "\n  Alpha Molecular orbitals:\n");
        Ca_->eivprint(eigvaluesa);

        fprintf(outfile, "\n  Beta Molecular orbitals:\n");
        Cb_->eivprint(eigvaluesb);
    }

    // Print out orbital energies.
    std::vector<std::pair<double, int> > pairsa, pairsb;
    for (int h=0; h<eigvaluesa->nirreps(); ++h) {
        for (int i=0; i<eigvaluesa->dimpi()[h]; ++i) {
            pairsa.push_back(make_pair(eigvaluesa->get(h, i), h));
            pairsb.push_back(make_pair(eigvaluesb->get(h, i), h));
        }
    }
    sort(pairsa.begin(),pairsa.end());
    sort(pairsb.begin(),pairsb.end());
    if(print_ > 1){
        fprintf(outfile, "\n  Orbital energies (a.u.):\n    Alpha occupied\n      ");
        for (int i=1; i<=nalpha_; ++i) {
            fprintf(outfile, "%12.6f %3s  ", pairsa[i-1].first, temp2[pairsa[i-1].second]);
            if (i % 4 == 0)
                fprintf(outfile, "\n      ");
        }
        fprintf(outfile, "\n");
        fprintf(outfile, "\n    Alpha unoccupied\n      ");
        for (int i=nalpha_+1; i<=nso; ++i) {
            fprintf(outfile, "%12.6f %3s  ", pairsa[i-1].first, temp2[pairsa[i-1].second]);
            if ((i-nalpha_) % 4 == 0)
                fprintf(outfile, "\n      ");
        }
        fprintf(outfile, "\n");

        fprintf(outfile, "\n    Beta occupied\n      ");
        for (int i=1; i<=nbeta_; ++i) {
            fprintf(outfile, "%12.6f %3s  ", pairsb[i-1].first, temp2[pairsb[i-1].second]);
            if (i % 4 == 0)
                fprintf(outfile, "\n      ");
        }
        fprintf(outfile, "\n");
        fprintf(outfile, "\n    Beta unoccupied\n      ");
        for (int i=nalpha_+1; i<=nso; ++i) {
            fprintf(outfile, "%12.6f %3s  ", pairsb[i-1].first, temp2[pairsb[i-1].second]);
            if ((i-nbeta_) % 4 == 0)
                fprintf(outfile, "\n      ");
        }
        fprintf(outfile, "\n");
    }
    for (int i=0; i<eigvaluesa->nirreps(); ++i)
        free(temp2[i]);
    free(temp2);

    int *vec = new int[eigvaluesa->nirreps()];
    for (int i=0; i<eigvaluesa->nirreps(); ++i)
        vec[i] = 0;

    chkpt_->wt_nmo(nso);
    chkpt_->wt_ref(1);        // UHF
    chkpt_->wt_etot(E_);
    chkpt_->wt_escf(E_);
    chkpt_->wt_eref(E_);
    chkpt_->wt_clsdpi(doccpi_);
    chkpt_->wt_orbspi(eigvaluesa->dimpi());
    chkpt_->wt_openpi(vec);
    chkpt_->wt_phase_check(0);

    // Figure out frozen core orbitals
    int nfzc = chkpt_->rd_nfzc();
    int nfzv = chkpt_->rd_nfzv();
    int *frzcpi = compute_fcpi(nfzc, eigvaluesa);
    int *frzvpi = compute_fvpi(nfzv, eigvaluesa);
    chkpt_->wt_frzcpi(frzcpi);
    chkpt_->wt_frzvpi(frzvpi);
    delete[](frzcpi);
    delete[](frzvpi);

    // TODO: Figure out what chkpt_wt_iopen means for UHF
    chkpt_->wt_iopen(0);

    // Write eigenvectors and eigenvalue to checkpoint
    double *values = eigvaluesa->to_block_vector();
    chkpt_->wt_alpha_evals(values);
    free(values);
    double **vectors = Ca_->to_block_matrix();
    chkpt_->wt_alpha_scf(vectors);
    free_block(vectors);
    values = eigvaluesb->to_block_vector();
    chkpt_->wt_beta_evals(values);
    free(values);
    vectors = Cb_->to_block_matrix();
    chkpt_->wt_beta_scf(vectors);
    free_block(vectors);
}

bool UHF::test_convergency()
{
    double ediff = E_ - Eold_;

    // RMS of the density
    Matrix Drms;
    Drms.copy(Dt_);
    Drms.subtract(Dtold_);
    Drms_ = sqrt(Drms.sum_of_squares());

    if (fabs(ediff) < energy_threshold_ && Drms_ < density_threshold_)
        return true;
    else
        return false;
}

void UHF::allocate_PK() {
    // Figure out how many pair combinations yield A1 symmetry (done in above loop)
    //   num_pair_combinations_of_A1 = ioff[_opi[0]] + ioff[_opi[1]] + ioff[_opi[2]] + ...
    // Allocate memory for the PK matrix (using a vector)
    if (pk_size_ < (memory_ / sizeof(double) / 2)) {
        p_jk_ = new double[pk_size_];
        p_k_ = new double[pk_size_];

        if (p_jk_ == NULL || p_k_ == NULL) {
            fprintf(outfile, "  Insufficient free system memory for in-core PK implementation.\n");
            fprintf(outfile, "  Switching to out-of-core algorithm.\n");
            use_out_of_core_ = true;
        } else {
            // PK and K are zeroed out just before filling, not here
            if(print_ > 2){
                fprintf(outfile,
                "  Allocated %lu elements (%lu pairs) for PJ. (%5f MiB)\n",
                (unsigned long)pk_size_, (unsigned long)pk_pairs_, pk_size_ * 8.0 / 1048576.0);
                fprintf(outfile,
                "  Allocated %lu elements (%lu pairs) for PK. (%5f MiB)\n\n",
                (unsigned long)pk_size_, (unsigned long)pk_pairs_, pk_size_ * 8.0 / 1048576.0);
            }
        }
    } else {
        fprintf(outfile,
                        "  Insufficient memory for in-core PK implementation.\n");
        fprintf(outfile,
                        "  Would need %lu elements of double memory. (%5f MiB)\n",
                        (unsigned long)pk_size_*2, pk_size_ * 8.0 / 1048576.0 * 2.0);
        fprintf(outfile, "  Switching to out-of-core algorithm.\n");
        use_out_of_core_ = true;
    }
}

void UHF::form_initialF()
{
    Fa_->copy(H_);
    Fb_->copy(H_);

    // Transform the Focks
    // Fa_->transform(Shalf_);
    // Fb_->transform(Shalf_);

#ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "Initial Fock alpha matrix:\n");
        Fa_->print(outfile);
        fprintf(outfile, "Initial Fock beta matrix:\n");
        Fb_->print(outfile);
    }
#endif
}

void UHF::form_F() {
    Fa_->copy(H_);
    Fa_->add(Ga_);

    Fb_->copy(H_);
    Fb_->add(Gb_);

    if(addExternalPotential_){
        // If there's an external potential, its contribution has already been computed
        // in the form G routine and stored in pertF.  Now we need to add the rest of the
        // Fock matrix contribution
        pertFa_->add(Fa_);
        pertFb_->add(Fb_);
    }

#ifdef _DEBUG
    if (debug_) {
        Fa_->print(outfile);
        Fb_->print(outfile);
    }
#endif
}

void UHF::form_C()
{
    Matrix eigvec;
    Vector eigval;
        factory_.create_matrix(eigvec);
        factory_.create_vector(eigval);

        if(addExternalPotential_){
            pertFa_->transform(Shalf_);
            pertFa_->diagonalize(eigvec, eigval);
        }else{
            Fa_->transform(Shalf_);
            Fa_->diagonalize(eigvec, eigval);
        }

        find_occupation(eigval);
    // fprintf(outfile, "Fa eigenvectors/values:\n");
    // eigvec->eivprint(eigval);
    Ca_->gemm(false, false, 1.0, Shalf_, eigvec, 0.0);

        if(addExternalPotential_){
            pertFb_->transform(Shalf_);
            pertFb_->diagonalize(eigvec, eigval);
        }else{
            Fb_->transform(Shalf_);
            Fb_->diagonalize(eigvec, eigval);
        }
    // fprintf(outfile, "Fb eigenvectors/values:\n");
    // eigvec->eivprint(eigval);
    Cb_->gemm(false, false, 1.0, Shalf_, eigvec, 0.0);

#ifdef _DEBUG
    if (debug_) {
        Ca_->print(outfile);
        Cb_->print(outfile);
    }
#endif
}

void UHF::form_D()
{
    int h, i, j, m;
    int *opi = Da_->rowspi();
    int nirreps = Da_->nirreps();
    double val;
    for (h=0; h<nirreps; ++h) {
    // fprintf(outfile, "irrep %d: nalpha = %d nbeta = %d\n", h, nalphapi_[h], nbetapi_[h]);
        for (i=0; i<opi[h]; ++i) {
            for (j=0; j<opi[h]; ++j) {
                val = 0.0;
                for (m=0; m<nalphapi_[h]; ++m)
                    val += Ca_->get(h, i, m) * Ca_->get(h, j, m);
                Da_->set(h, i, j, val);

                val = 0.0;
                for (m=0; m<nbetapi_[h]; ++m)
                    val += Cb_->get(h, i, m) * Cb_->get(h, j, m);
                Db_->set(h, i, j, val);
            }
        }
    }

    // Form total density
    // TODO: Refactor Dt_ to D_ (found in HF)
    Dt_->copy(Da_);
    Dt_->add(Db_);

#ifdef _DEBUG
    if (debug_) {
        Da_->print(outfile);
        Db_->print(outfile);
    }
#endif
}

// TODO: Once Dt_ is refactored to D_ the only difference between this and RHF::compute_initial_E is a factor of 0.5
double UHF::compute_initial_E()
{
    double Etotal = nuclearrep_ + 0.5 * (Dt_->vector_dot(H_));
    if(print_ > 2) fprintf(outfile, "\n  Initial UHF energy: %20.14f\n\n", Etotal);
    fflush(outfile);
    return Etotal;
}

double UHF::compute_E()
{
    double DH  = Dt_->vector_dot(H_);
    double DFa = Da_->vector_dot(Fa_);
    double DFb = Db_->vector_dot(Fb_);
    double Eelec = 0.5 * (DH + DFa + DFb);
    // fprintf(outfile, "electronic energy = %20.14f\n", Eelec);
    double Etotal = nuclearrep_ + Eelec;
    return Etotal;
}

void UHF::form_PK()
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
    double value;

    // PK is zeroed out here, not during allocation
    // to prevent issues with multiple SCF calls
    memset(p_jk_, 0, pk_size_*sizeof(double));
    memset(p_k_, 0, pk_size_*sizeof(double));

    if(print_ > 2){
        fprintf(outfile, "  Forming PJ and PK matrices.\n");
        fflush(outfile);
    }

    IWL ERIIN(psio_.get(), PSIF_SO_TEI, 0.0, 1, 1);

    do {
        ilsti = ERIIN.last_buffer();
        nbuf  = ERIIN.buffer_count();

        fi = 0;
        for (idx=0; idx<nbuf; ++idx) {
            if (ERIIN.labels()[fi] >= 0) {
                i = ERIIN.labels()[fi];
            }
            else {
                i = -ERIIN.labels()[fi];
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
                // pk_symoffset_ corrects for the symmetry offset in the _pk vector
                braket = INDEX2(bra, ket);
                p_jk_[braket] += value;
                // K/2 (2nd sort)
                if ((ii != jj) && (kk != ll)) {
                    if ((is == ls) && (js == ks)) {
                        bra = INDEX2(ii, ll) + pk_symoffset_[is];
                        ket = INDEX2(jj, kk) + pk_symoffset_[js];
                        braket = INDEX2(bra, ket);
                        if ((ii == ll) || (jj == kk)) {
                            p_jk_[braket] -= 0.5 * value;
                            p_k_[braket]  -= 0.5 * value;
                        }
                        else {
                            p_jk_[braket] -= 0.25 * value;
                            p_k_[braket]  -= 0.25 * value;
                        }
                    }
                }
            }

            // K/2 (1st sort)
            if ((is == ks) && (js == ls)) {
                bra = INDEX2(ii, kk) + pk_symoffset_[is];
                ket = INDEX2(jj, ll) + pk_symoffset_[js];
                braket = INDEX2(bra, ket);
                if ((ii == kk) || (jj == ll)) {
                    p_jk_[braket] -= 0.5 * value;
                    p_k_[braket]  -= 0.5 * value;
                }
                else {
                    p_jk_[braket] -= 0.25 * value;
                    p_k_[braket]  -= 0.25 * value;
                }
            }
            counter++;
        }
        if (!ilsti)
            ERIIN.fetch();
    } while (!ilsti);

    // Going out of scope will close the buffer
    // iwl_buf_close(&ERIIN, 1);

    // After stage two is complete, the elements of P must be halved for the case IJ=KL.
    for (size_t ij=0; ij < pk_pairs_; ++ij) {
        p_jk_[INDEX2(ij,ij)] *= 0.5;
        p_k_[INDEX2(ij,ij)] *= 0.5;
    }

    if(print_ > 2) fprintf(outfile, "  Processed %d two-electron integrals.\n", counter);
    #ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "p_jk_:\n");
        print_array(p_jk_, pk_pairs_, outfile);
        fprintf(outfile, "p_k_:\n");
        print_array(p_k_, pk_pairs_, outfile);
    }
    #endif
}

void UHF::form_G_from_PK()
{
    int nirreps = factory_.nirreps();
    int *opi = factory_.rowspi();
    size_t ij;
    double *Da_vector = new double[pk_pairs_];
    double *Db_vector = new double[pk_pairs_];
    double *Ga_vector = new double[pk_pairs_];
    double *Gb_vector = new double[pk_pairs_];
    double *Va_vector = NULL;
    double *Vb_vector = NULL;
    double *Fa_vector = NULL;
    double *Fb_vector = NULL;
    if(addExternalPotential_){
        Va_vector = new double[pk_pairs_];
        Vb_vector = new double[pk_pairs_];
        Fa_vector = new double[pk_pairs_];
        Fb_vector = new double[pk_pairs_];
        memset(Va_vector, 0, sizeof(double) * pk_pairs_);
        memset(Vb_vector, 0, sizeof(double) * pk_pairs_);
        memset(Fa_vector, 0, sizeof(double) * pk_pairs_);
        memset(Fb_vector, 0, sizeof(double) * pk_pairs_);
    }

    Ga_->zero();
    Gb_->zero();

    memset(Da_vector, 0, sizeof(double) * pk_pairs_);
    memset(Db_vector, 0, sizeof(double) * pk_pairs_);
    memset(Ga_vector, 0, sizeof(double) * pk_pairs_);
    memset(Gb_vector, 0, sizeof(double) * pk_pairs_);

    ij=0;
    for (int h=0; h<nirreps; ++h) {
        for (int p=0; p<opi[h]; ++p) {
            for (int q=0; q<=p; ++q) {
                if (p != q) {
                    Da_vector[ij] = 2.0 * Da_->get(h, p, q);
                    Db_vector[ij] = 2.0 * Db_->get(h, p, q);
                    // Add in the external potential, if requested
                    if(addExternalPotential_) Va_vector[ij] = 2.0 * Va_[h][p][q];
                    if(addExternalPotential_) Vb_vector[ij] = 2.0 * Vb_[h][p][q];

                } else {
                    Da_vector[ij] = Da_->get(h, p, q);
                    Db_vector[ij] = Db_->get(h, p, q);
                    if(addExternalPotential_) Va_vector[ij] = Va_[h][p][q];
                    if(addExternalPotential_) Vb_vector[ij] = Vb_[h][p][q];
                }
                ij++;
            }
        }
    }

#ifdef _DEBUG
    if (debug_) {
            fprintf(outfile, "PK: ij = %lu\n", (unsigned long)ij);
            fflush(outfile);
            fprintf(outfile, "PK: Da matrix:\n");
            Da_->print(outfile);
            fprintf(outfile, "PK: Da vector (appears to be OK):\n");
            for (ij=0; ij<pk_pairs_; ++ij)
                    fprintf(outfile, "PK: Da vector [%lu] = %20.14f\n", (unsigned long)ij, Da_vector[ij]);
            fprintf(outfile, "PK: Db matrix:\n");
            Db_->print(outfile);
            fprintf(outfile, "PK: Db vector (appears to be OK):\n");
            for (ij=0; ij<pk_pairs_; ++ij)
                    fprintf(outfile, "PK: Db vector [%lu] = %20.14f\n", (unsigned long)ij, Db_vector[ij]);
    }
#endif

    /*
     * This code goes through the densities (Da_ and Db_), J, and K to form
     * two G matrices. One G matrix is for Fa_ and the other for Fb_.
     * See derivation notebook for equations.
     */
    double Ga_pq, Da_pq, Fa_pq, Va_pq;
    double Gb_pq, Db_pq, Fb_pq, Vb_pq;
    double* Da_rs;
    double* Ga_rs;
    double* Db_rs;
    double* Gb_rs;
    double* Fa_rs;
    double* Fb_rs;
    double* Va_rs;
    double* Vb_rs;
    int pq, rs;
    double* JK_block = p_jk_;
    double* K_block = p_k_;
    int ts_pairs = pk_pairs_;
    for (pq = 0; pq < ts_pairs; ++pq) {
        Ga_pq = 0.0;
        Da_pq = Da_vector[pq];
        Da_rs = &Da_vector[0];
        Ga_rs = &Ga_vector[0];
        Gb_pq = 0.0;
        Db_pq = Db_vector[pq];
        Db_rs = &Db_vector[0];
        Gb_rs = &Gb_vector[0];
        if(addExternalPotential_){
            Fa_pq = 0.0;
            Va_pq = Va_vector[pq];
            Va_rs = &Va_vector[0];
            Fa_rs = &Fa_vector[0];
            Fb_pq = 0.0;
            Vb_pq = Vb_vector[pq];
            Vb_rs = &Vb_vector[0];
            Fb_rs = &Fb_vector[0];
        }
        for (rs = 0; rs <= pq; ++rs) {
            // D_{rs}^{c} * PK_{pqrs}         Also found in RHF
            // Doing F_mn_a about to add the K term
            Ga_pq  += (*JK_block + *K_block) * (*Da_rs) + (*JK_block - *K_block) * (*Db_rs);
            *Ga_rs += (*JK_block + *K_block) * Da_pq    + (*JK_block - *K_block) * Db_pq;

            Gb_pq  += (*JK_block + *K_block) * (*Db_rs) + (*JK_block - *K_block) * (*Da_rs);
            *Gb_rs += (*JK_block + *K_block) * Db_pq    + (*JK_block - *K_block) * Da_pq;
            if(addExternalPotential_){
                Fa_pq  += (*JK_block + *K_block) * (*Va_rs) + (*JK_block - *K_block) * (*Vb_rs);
                *Fa_rs += (*JK_block + *K_block) * Va_pq    + (*JK_block - *K_block) * Vb_pq;

                Fb_pq  += (*JK_block + *K_block) * (*Vb_rs) + (*JK_block - *K_block) * (*Va_rs);
                *Fb_rs += (*JK_block + *K_block) * Vb_pq    + (*JK_block - *K_block) * Va_pq;
                ++Fa_rs;
                ++Fb_rs;
                ++Va_rs;
                ++Vb_rs;
            }
            ++Da_rs;
            ++Ga_rs;
            ++Db_rs;
            ++Gb_rs;
            ++JK_block;
            ++K_block;
        }
        Ga_vector[pq] += Ga_pq;
        Gb_vector[pq] += Gb_pq;
        if(addExternalPotential_){
            Fa_vector[pq] += Fa_pq;
            Fb_vector[pq] += Fb_pq;
        }
    }

    // Convert G to a matrix
    ij = 0;
    for (int h = 0; h < nirreps; ++h) {
        for (int p = 0; p < opi[h]; ++p) {
            for (int q = 0; q <= p; ++q) {
                Ga_->set(h, p, q, Ga_vector[ij]);
                Ga_->set(h, q, p, Ga_vector[ij]);
                Gb_->set(h, p, q, Gb_vector[ij]);
                Gb_->set(h, q, p, Gb_vector[ij]);
                if(addExternalPotential_){
                    pertFa_->set(h, p, q, Fa_vector[ij]);
                    pertFa_->set(h, q, p, Fa_vector[ij]);
                    pertFb_->set(h, p, q, Fb_vector[ij]);
                    pertFb_->set(h, q, p, Fb_vector[ij]);
                }
                ij++;
            }
        }
    }

#ifdef _DEBUG
    if (debug_) {
            Ga_->print(outfile);
            Gb_->print(outfile);
    }
#endif

    if(Va_vector != NULL) delete[](Va_vector);
    if(Vb_vector != NULL) delete[](Vb_vector);
    if(Fa_vector != NULL) delete[](Fa_vector);
    if(Fb_vector != NULL) delete[](Fb_vector);
    delete[](Da_vector);
    delete[](Db_vector);
    delete[](Ga_vector);
    delete[](Gb_vector);
}

void UHF::form_G()
{
    fprintf(stderr, "UHF out-of-core algorithm is not implemented yet!\n");
    abort();
}
void UHF::form_G_from_direct_integrals()
{
    fprintf(stderr, "UHF integral direct algorithm is not implemented yet!\n");
    abort();
}

void UHF::save_fock()
{
    static bool initialized_diis_manager = false;
    if (initialized_diis_manager == false) {
        diis_manager_->set_error_vector_size(2,
                                             DIISEntry::Matrix, Fa_.get(),
                                             DIISEntry::Matrix, Fb_.get());
        diis_manager_->set_vector_size(2,
                                       DIISEntry::Matrix, Fa_.get(),
                                       DIISEntry::Matrix, Fb_.get());
        initialized_diis_manager = true;
    }

    // Save the current Fock matrix
//    diis_F_[current_diis_fock_]->copy(F_);

    // Determine error matrix for this Fock
    SharedMatrix FaDaS(factory_.create_matrix()), DaS(factory_.create_matrix());
    SharedMatrix SDaFa(factory_.create_matrix()), DaFa(factory_.create_matrix());
    SharedMatrix FbDbS(factory_.create_matrix()), DbS(factory_.create_matrix());
    SharedMatrix SDbFb(factory_.create_matrix()), DbFb(factory_.create_matrix());

    // FDS = F_ * D_ * S_; Alpha
    DaS->gemm(false, false, 1.0, Da_, S_, 0.0);
    FaDaS->gemm(false, false, 1.0, Fa_, DaS, 0.0);
    // SDF = S_ * D_ * F_;
    DaFa->gemm(false, false, 1.0, Da_, Fa_, 0.0);
    SDaFa->gemm(false, false, 1.0, S_, DaFa, 0.0);

    // FDS = F_ * D_ * S_; Beta
    DbS->gemm(false, false, 1.0, Db_, S_, 0.0);
    FbDbS->gemm(false, false, 1.0, Fb_, DbS, 0.0);
    // SDF = S_ * D_ * F_;
    DbFb->gemm(false, false, 1.0, Db_, Fb_, 0.0);
    SDbFb->gemm(false, false, 1.0, S_, DbFb, 0.0);

    Matrix FaDaSmSDaFa;
    FaDaSmSDaFa.copy(FaDaS);
    FaDaSmSDaFa.subtract(SDaFa);
    FaDaSmSDaFa.transform(Shalf_);

    Matrix FbDbSmSDbFb;
    FbDbSmSDbFb.copy(FbDbS);
    FbDbSmSDbFb.subtract(SDbFb);
    FbDbSmSDbFb.transform(Shalf_);

    diis_manager_->add_entry(4, &FaDaSmSDaFa, &FbDbSmSDbFb, Fa_.get(), Fb_.get());
}

void UHF::diis()
{
    diis_manager_->extrapolate(2, Fa_.get(), Fb_.get());
}

}}
