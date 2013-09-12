#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libmints/view.h>
#include <libmints/local.h>
#include <libdpd/dpd.h>
#include <psifiles.h>

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

namespace psi{ namespace lmp2 {

class Lmp2 : public Wavefunction
{
public:
    Lmp2(boost::shared_ptr<Wavefunction> reference_wavefunction, Options& options);
    virtual ~Lmp2();

    double compute_energy();

private:
    Dimension nvirtpi_;
    void common_init();

    void denom(const Dimension& occOrbsPI,
               const Dimension& virOrbsPI,
               const Dimension& occOffset,
               const Dimension& virOffset);

    void pair_energies(IntegralTransform& ints, const Dimension& occOrbsPI, double** epair);
    void print_pair_energies(const Dimension& occOrbsPI, double* emp2);
};

Lmp2::Lmp2(boost::shared_ptr<Wavefunction> reference_wavefunction, Options& options)
    : Wavefunction(options, _default_psio_lib_)
{
//    Process::environment.set_wavefunction(reference_wavefunction);
    reference_wavefunction_ = reference_wavefunction;
    common_init();
}

Lmp2::~Lmp2()
{
}

/* pair_energies(): For RHF references, compute pair energies. Spin-adapt
** pair energies if SPINADAPT_ENERGIES is set to true.
**
** E(IJ) = T2(IJ,AB) * (<ij|ab> - <ij|ba>)
** E(Ij) = T2(Ij,Ab) * <ij|ab>
**
*/
void Lmp2::pair_energies(IntegralTransform& ints, const Dimension& occOrbsPI, double** epair)
{
    dpdbuf4 tau, D, E;

    int i, j, ij;
    int irrep;
    int nocc_act = 0;
    int nab;
    nocc_act = occOrbsPI.sum();
    nab = nocc_act * nocc_act;

    double* pair = new double[nab];

    for(int p=0; p<nab; p++){
        pair[p] = 0.0;
    }

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O>=O]+"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, ID("[O>=O]+"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "LMP2 tIjAb");
    global_dpd_->buf4_init(&E, PSIF_CC_TAMPS, 0, ID("[O>=O]+"), ID("[O>=O]+"), ID("[O>=O]+"), ID("[O>=O]+"), 0, "E <ij|kl>");
    global_dpd_->contract444(&D, &tau, &E, 0, 0, 1.0, 0.0);
    //global_dpd_->buf4_print(&E, outfile, 1);

    /* Extract diagonal elements (i.e. pair energies) and print them out nicely */
    for(irrep=0; irrep<nirrep_; irrep++) {
        double **block;
        dpdparams4 *Params = E.params;
        int p;
        int np = Params->rowtot[irrep];

        global_dpd_->buf4_mat_irrep_init(&E, irrep);
        global_dpd_->buf4_mat_irrep_rd(&E, irrep);
        block = E.matrix[irrep];

        for(p=0; p<np; p++) {
            int i, j, ij;

            i = Params->roworb[irrep][p][0];
            j = Params->roworb[irrep][p][1];

            pair[p] = block[p][p];
        }
        global_dpd_->buf4_mat_irrep_close(&E, irrep);
    }

    *epair = pair;

    global_dpd_->buf4_close(&tau);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&E);
}

void Lmp2::print_pair_energies(const Dimension& occOrbsPI, double* emp2)
{
    int i, j, ij;
    int irrep;
    int nocc_act = 0;
    int naa, nab;
    nocc_act = occOrbsPI.sum();
    naa = nocc_act * (nocc_act-1)/2;
    nab = nocc_act * nocc_act;

    double emp2_tot = 0.0;

    fprintf(outfile, "\tOrbital pair energies\n");
    fprintf(outfile, "\t    i       j         LMP2\n");
    fprintf(outfile, "\t  -----   -----   ------------\n");
    ij = 0;
    for(i=0; i<nocc_act; i++)
        for(j=0; j<=i; j++,ij++) {
            fprintf(outfile, "\t  %3d     %3d     %12.9lf\n", i+1, j+1, emp2[ij]);
            emp2_tot += emp2[ij];
            if (i != j)
                emp2_tot += emp2[ij];
        }
    fprintf(outfile, "\t  -------------   ------------\n");
    fprintf(outfile, "\t      Total       %12.9lf\n\n", emp2_tot);


    fprintf(outfile, "\n");
}

void Lmp2::denom(const Dimension& occOrbsPI,
                 const Dimension& virOrbsPI,
                 const Dimension& occOffset,
                 const Dimension& virOffset)
{
    Dimension openpi(nirrep_);

    int h, i, j, a, b, ij, ab;
    int I, J, A, B;
    int isym, jsym, asym, bsym;
    double fii, fjj, faa, fbb;
    dpdfile2 fIJ, fij, fAB, fab;
    dpdfile2 dIA, dia;
    dpdfile4 dIJAB, dijab, dIjAb;

    /* Grab Fock matrices from disk */
    global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_mat_init(&fIJ);
    global_dpd_->file2_mat_rd(&fIJ);

    global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
    global_dpd_->file2_mat_init(&fij);
    global_dpd_->file2_mat_rd(&fij);

    global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_mat_init(&fAB);
    global_dpd_->file2_mat_rd(&fAB);

    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
    global_dpd_->file2_mat_init(&fab);
    global_dpd_->file2_mat_rd(&fab);

    /* Alpha one-electron denominator */
    global_dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
    global_dpd_->file2_mat_init(&dIA);

    for(h=0; h < nirrep_; h++) {

        for(i=0; i < occOrbsPI[h]; i++) {
            fii = fIJ.matrix[h][i][i];

            for(a=0; a < (virOrbsPI[h] - openpi[h]); a++) {
                faa = fAB.matrix[h][a][a];

                dIA.matrix[h][i][a] = 1.0/(fii - faa);
            }
        }
    }

    global_dpd_->file2_mat_wrt(&dIA);
    global_dpd_->file2_mat_close(&dIA);
    global_dpd_->file2_close(&dIA);

    /* Beta one-electron denominator */
    global_dpd_->file2_init(&dia, PSIF_CC_OEI, 0, 0, 1, "dia");
    global_dpd_->file2_mat_init(&dia);

    for(h=0; h < nirrep_; h++) {

        for(i=0; i < (occOrbsPI[h] - openpi[h]); i++) {
            fii = fij.matrix[h][i][i];

            for(a=0; a < virOrbsPI[h]; a++) {
                faa = fab.matrix[h][a][a];

                dia.matrix[h][i][a] = 1.0/(fii - faa);
            }
        }
    }

    global_dpd_->file2_mat_wrt(&dia);
    global_dpd_->file2_mat_close(&dia);
    global_dpd_->file2_close(&dia);

    /* Alpha-alpha two-electron denominator */
    global_dpd_->file4_init(&dIJAB, PSIF_CC_DENOM, 0, 1, 6, "dIJAB");

    for(h=0; h < nirrep_; h++) {

        global_dpd_->file4_mat_irrep_init(&dIJAB, h);

        /* Loop over the rows */
        for(ij=0; ij < dIJAB.params->rowtot[h]; ij++) {
            i = dIJAB.params->roworb[h][ij][0];
            j = dIJAB.params->roworb[h][ij][1];
            isym = dIJAB.params->psym[i];
            jsym = dIJAB.params->qsym[j];

            /* Convert to relative orbital index */
            I = i - occOffset[isym];
            J = j - occOffset[jsym];

            fii = fIJ.matrix[isym][I][I];
            fjj = fIJ.matrix[jsym][J][J];

            /* Loop over the columns */
            for(ab=0; ab < dIJAB.params->coltot[h]; ab++) {
                a = dIJAB.params->colorb[h][ab][0];
                b = dIJAB.params->colorb[h][ab][1];
                asym = dIJAB.params->rsym[a];
                bsym = dIJAB.params->ssym[b];

                /* Convert to relative orbital index */
                A = a - virOffset[asym];
                B = b - virOffset[bsym];

                faa = fAB.matrix[asym][A][A];
                fbb = fAB.matrix[bsym][B][B];

                dIJAB.matrix[h][ij][ab] =
                    ((A >= (virOrbsPI[asym] - openpi[asym])) ||
                     (B >= (virOrbsPI[bsym] - openpi[bsym])) ?
                     0.0 : 1.0/(fii + fjj - faa - fbb));
            }
        }

        global_dpd_->file4_mat_irrep_wrt(&dIJAB, h);
        global_dpd_->file4_mat_irrep_close(&dIJAB, h);

    }

    global_dpd_->file4_close(&dIJAB);

    /* Beta-beta two-electron denominator */
    global_dpd_->file4_init(&dijab, PSIF_CC_DENOM, 0, 1, 6, "dijab");

    for(h=0; h < nirrep_; h++) {

        global_dpd_->file4_mat_irrep_init(&dijab, h);

        /* Loop over the rows */
        for(ij=0; ij < dijab.params->rowtot[h]; ij++) {
            i = dijab.params->roworb[h][ij][0];
            j = dijab.params->roworb[h][ij][1];
            isym = dijab.params->psym[i];
            jsym = dijab.params->qsym[j];

            /* Convert to relative orbital index */
            I = i - occOffset[isym];
            J = j - occOffset[jsym];

            fii = fij.matrix[isym][I][I];
            fjj = fij.matrix[jsym][J][J];

            /* Loop over the columns */
            for(ab=0; ab < dijab.params->coltot[h]; ab++) {
                a = dijab.params->colorb[h][ab][0];
                b = dijab.params->colorb[h][ab][1];
                asym = dijab.params->rsym[a];
                bsym = dijab.params->ssym[b];

                /* Convert to relative orbital index */
                A = a - virOffset[asym];
                B = b - virOffset[bsym];

                faa = fab.matrix[asym][A][A];
                fbb = fab.matrix[bsym][B][B];

                dijab.matrix[h][ij][ab] =
                    ((I >= (occOrbsPI[isym] - openpi[isym])) ||
                     (J >= (occOrbsPI[jsym] - openpi[jsym])) ?
                     0.0 : 1.0/(fii + fjj - faa - fbb));
            }
        }

        global_dpd_->file4_mat_irrep_wrt(&dijab, h);
        global_dpd_->file4_mat_irrep_close(&dijab, h);

    }

    global_dpd_->file4_close(&dijab);

    /* Alpha-beta two-electron denominator */
    global_dpd_->file4_init(&dIjAb, PSIF_CC_DENOM, 0, 0, 5, "dIjAb");

    for(h=0; h < nirrep_; h++) {

        global_dpd_->file4_mat_irrep_init(&dIjAb, h);

        /* Loop over the rows */
        for(ij=0; ij < dIjAb.params->rowtot[h]; ij++) {
            i = dIjAb.params->roworb[h][ij][0];
            j = dIjAb.params->roworb[h][ij][1];
            isym = dIjAb.params->psym[i];
            jsym = dIjAb.params->qsym[j];

            /* Convert to relative orbital index */
            I = i - occOffset[isym];
            J = j - occOffset[jsym];

            fii = fIJ.matrix[isym][I][I];
            fjj = fij.matrix[jsym][J][J];

            /* Loop over the columns */
            for(ab=0; ab < dIjAb.params->coltot[h]; ab++) {
                a = dIjAb.params->colorb[h][ab][0];
                b = dIjAb.params->colorb[h][ab][1];
                asym = dIjAb.params->rsym[a];
                bsym = dIjAb.params->ssym[b];

                /* Convert to relative orbital index */
                A = a - virOffset[asym];
                B = b - virOffset[bsym];

                faa = fAB.matrix[asym][A][A];
                fbb = fab.matrix[bsym][B][B];

                dIjAb.matrix[h][ij][ab] =
                    ((A >= (virOrbsPI[asym] - openpi[asym])) ||
                     (J >= (occOrbsPI[jsym] - openpi[jsym])) ?
                     0.0 : 1.0/(fii + fjj - faa - fbb));
            }
        }

        global_dpd_->file4_mat_irrep_wrt(&dIjAb, h);
        global_dpd_->file4_mat_irrep_close(&dIjAb, h);

    }

    global_dpd_->file4_close(&dIjAb);

    global_dpd_->file2_mat_close(&fIJ);
    global_dpd_->file2_mat_close(&fij);
    global_dpd_->file2_mat_close(&fAB);
    global_dpd_->file2_mat_close(&fab);
    global_dpd_->file2_close(&fIJ);
    global_dpd_->file2_close(&fij);
    global_dpd_->file2_close(&fAB);
    global_dpd_->file2_close(&fab);
}

void Lmp2::common_init()
{
    nso_     = reference_wavefunction_->nso();
    nirrep_  = reference_wavefunction_->nirrep();
    nmo_     = reference_wavefunction_->nmo();
    soccpi_  = reference_wavefunction_->soccpi();
    doccpi_  = reference_wavefunction_->doccpi();
    frzcpi_  = reference_wavefunction_->frzcpi();
    frzvpi_  = reference_wavefunction_->frzvpi();
    nmopi_   = reference_wavefunction_->nmopi();
    nsopi_   = reference_wavefunction_->nsopi();
    epsilon_a_ = reference_wavefunction_->epsilon_a();
    epsilon_b_ = reference_wavefunction_->epsilon_b();
    nvirtpi_ = nsopi_ - frzcpi_ - frzvpi_ - doccpi_;
}

double Lmp2::compute_energy()
{
    fprintf(outfile,"\n\n *******************************************************************************\n");
    fprintf(outfile,    " *                                  Local MP2                                  *\n");
    fprintf(outfile,    " *                    by Justin Turney and Brandon Magers                      *\n");
    fprintf(outfile,    " *           Parts taken from ccsort and ccenergy by Daniel Crawford           *\n");
    fprintf(outfile,    " *******************************************************************************\n\n");

    if (reference_wavefunction_->same_a_b_orbs() == false) {
        fprintf(outfile, "\tLMP2 only works for RHF reference.\n");
        throw PSIEXCEPTION("LMP2 only works for RHF reference.");
    }

    double convergence = Process::environment.options.get_double("E_CONVERGENCE");
    int maxiter = Process::environment.options.get_int("MAXITER");

    SharedMatrix C_focc = reference_wavefunction_->Ca_subset("SO", "FROZEN_OCC");
    SharedMatrix C_occ  = reference_wavefunction_->Ca_subset("SO", "ACTIVE_OCC");
    SharedMatrix C_vir  = reference_wavefunction_->Ca_subset("SO", "ACTIVE_VIR");
    SharedMatrix C_fvir = reference_wavefunction_->Ca_subset("SO", "FROZEN_VIR");

    boost::shared_ptr<Localizer> localA = Localizer::build(reference_wavefunction_->basisset(),
                                                           C_occ,
                                                           Process::environment.options);
    localA->localize();
    SharedMatrix LC_occ = localA->L();

    {
        std::vector<SharedMatrix> Cs;
        Cs.push_back(C_focc);
        Cs.push_back(LC_occ);
        Cs.push_back(C_vir);
        Cs.push_back(C_fvir);

        Ca_ = Matrix::horzcat(Cs);
    }

    // Quickly check that there are no open shell orbitals here...
    int h, i, j, a, b;
    int numAOcc, numBOcc, numAVir, numBVir;
    double eSCF;
    double *occEvals, *virEvals;

    SharedVector aEvals = reference_wavefunction_->epsilon_a();
    SharedVector bEvals = reference_wavefunction_->epsilon_b();
    char **labels       = reference_wavefunction_->molecule()->irrep_labels();
    int nmo             = reference_wavefunction_->nmo();
    Dimension mopi      = reference_wavefunction_->nmopi();
    Dimension clsdpi    = reference_wavefunction_->doccpi();
    Dimension openpi    = reference_wavefunction_->soccpi();
    Dimension frzcpi    = reference_wavefunction_->frzcpi();
    Dimension frzvpi    = reference_wavefunction_->frzvpi();
    Dimension occOrbsPI(nirrep_);
    Dimension virOrbsPI(nirrep_);
    Dimension occOffset(nirrep_);
    Dimension virOffset(nirrep_);
    numAOcc = 0; numBOcc = 0; numAVir = 0; numBVir = 0;

    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;
    for(h = 0; h < nirrep_; ++h){
        occOrbsPI[h] = clsdpi[h] - frzcpi[h];
        virOrbsPI[h] = mopi[h]   - clsdpi[h] - frzvpi[h];
        numAOcc += occOrbsPI[h];
        numAVir += virOrbsPI[h];
    }

    // form offsets
    int ocount = occOrbsPI[0], vcount = virOrbsPI[0];
    for (h=1; h<nirrep_; h++) {
        occOffset[h] = ocount;
        ocount += occOrbsPI[h];
        virOffset[h] = vcount;
        vcount += virOrbsPI[h];
    }

    occEvals = new double[numAOcc];
    virEvals = new double[numAVir];
    fprintf(outfile, "\n\n\tIrrep  Core  Docc  Socc  Occ  Vir\n");
    fprintf(outfile,     "\t===============================================\n");
    for(h = 0; h < nirrep_; ++h){
       fprintf(outfile, "\t %3s   %3d   %3d   %3d   %3d  %3d\n",
                             labels[h], frzcpi[h], clsdpi[h], openpi[h],
                             occOrbsPI[h], virOrbsPI[h]);
    }
    fprintf(outfile,     "\t===============================================\n\n");
    aOccCount = 0; aVirCount = 0;
    for(h = 0; h < nirrep_; ++h){
        for(a = frzcpi[h]; a < clsdpi[h]; ++a) occEvals[aOccCount++] = aEvals->get(h, a);
        for(a = clsdpi[h]; a < mopi[h]; ++a) virEvals[aVirCount++] = aEvals->get(h, a);
    }
    fprintf(outfile, "\tMAXITER = %d\n", maxiter);
#if 0
    if(print > 2){
        for(i = 0; i < numAOcc; ++i)
            fprintf(outfile, "\toccEvals[%2d] = %10.6f\n", i, occEvals[i]);
        fprintf(outfile, "\n");
        for(i = 0; i < numAVir; ++i)
            fprintf(outfile, "\tvirEvals[%2d] = %10.6f\n", i, virEvals[i]);
        fprintf(outfile, "\n");
    }
#endif

    // => Two-electron integral transformation <=

    // For now, we'll just transform for closed shells and generate all integrals.  For more elaborate use of the
    // LibTrans object, check out the plugin_mp2 example in the test suite.
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);

    IntegralTransform ints(C_focc, LC_occ, C_vir, C_fvir, spaces);

    {
        // Use the IntegralTransform object's DPD instance, for convenience
        dpd_set_default(ints.get_dpd_id());

        // Since we're doing multiple separate transformations, keep the SO ints around:
        ints.set_keep_dpd_so_ints(1);

        // Generate the integrals in various spaces in chemists' notation
        ints.set_dpd_int_file(PSIF_CC_AINTS);
        ints.set_aa_int_name("A (ik|jl)");
        fprintf(outfile, "\tTransforming integrals into types of A (ik|jl) ...\n\n"); fflush(outfile);
        ints.transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ, IntegralTransform::MakeAndKeep);
        fprintf(outfile, "\n");

        ints.set_dpd_int_file(PSIF_CC_CINTS);
        ints.set_aa_int_name("C (ij|ab)");
        fprintf(outfile, "\tTransforming integrals into types of C (ij|ab) ...\n\n"); fflush(outfile);
        ints.transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);
        fprintf(outfile, "\n");

        // This is the last transformation, we can delete the SO ints
        ints.set_keep_dpd_so_ints(0);
        ints.set_dpd_int_file(PSIF_CC_DINTS);
        ints.set_aa_int_name("D (ia|jb)");
        fprintf(outfile, "\tTransforming integrals into types of D (ia|jb) ...\n\n"); fflush(outfile);
        ints.transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);

        fprintf(outfile, "\n\tIntegral transformations complete.\n\n");

        // Open files we need
        psio_->open(PSIF_CC_OEI, PSIO_OPEN_NEW);
        psio_->open(PSIF_CC_AINTS, PSIO_OPEN_OLD);
        psio_->open(PSIF_CC_CINTS, PSIO_OPEN_OLD);
        psio_->open(PSIF_CC_DINTS, PSIO_OPEN_OLD);
        psio_->open(PSIF_CC_DENOM, PSIO_OPEN_OLD);
        psio_->open(PSIF_CC_TAMPS, PSIO_OPEN_OLD);

        // => Re-sort the chemists' notation integrals to physicists' notation (pq|rs) = <pr|qs> <=
        dpdbuf4 A, C, D;
        // Form A <ij|kl> from A (ik|jl)
        global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, ID("[O,O]"), ID("[O,O]"), ID("[O>=O]+"), ID("[O>=O]+"), 0, "A (ik|jl)");
        fprintf(outfile, "\tSorting A (ik|jl) to A <ij|kl> ... "); fflush(outfile);
        global_dpd_->buf4_sort(&A, PSIF_CC_AINTS, prqs, ID("[O,O]"), ID("[O,O]"), "A <ij|kl>");
        fprintf(outfile, "done\n"); fflush(outfile);
        global_dpd_->buf4_close(&A);

        // Form C <ia|jb> from C (ij|ab)
        global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>=O]+"), ID("[V>=V]+"), 0, "C (ij|ab)");
        fprintf(outfile, "\tSorting C (ij|ab) to C <ia|jb> ... "); fflush(outfile);
        global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, prqs, ID("[O,V]"), ID("[O,V]"), "C <ia|jb");
        fprintf(outfile, "done\n"); fflush(outfile);
        global_dpd_->buf4_close(&C);

        // Form D <ij|ab> from D (ia|jb)
        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "D (ia|jb)");
        fprintf(outfile, "\tSorting D (ia|jb) to D <ij|ab> ... "); fflush(outfile);
        global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, prqs, ID("[O,O]"), ID("[V,V]"), "D <ij|ab>");
        fprintf(outfile, "done\n");
        global_dpd_->buf4_close(&D);

        // Form D 2<ij|ab> - <ij|ba>
        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D <ij|ab>");
        fprintf(outfile, "\tForming D 2<ij|ab> - <ij|ba> ... "); fflush(outfile);
        global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D 2<ij|ab> - <ij|ba>");
        fprintf(outfile, "<ij|ab> done ... "); fflush(outfile);
        {
            dpdbuf4 Dba;
            global_dpd_->buf4_init(&Dba, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 1, "D <ij|ab>");
            global_dpd_->buf4_axpy(&Dba, &D, 1.0);
            fprintf(outfile, "<ij|ab> - <ij|ba> done ... "); fflush(outfile);
            global_dpd_->buf4_close(&Dba);
        }
        fprintf(outfile, "done\n");
        global_dpd_->buf4_close(&D);
    }

    // Form the Fock matrices
    {
        fprintf(outfile, "\n\tViewing Fock matrices ... "); fflush(outfile);

        // Looks like libtrans did this for us. Load it up.
        SharedMatrix f(new Matrix("f(m,n)", mopi, mopi));
        f->load(psio_, PSIF_OEI, PSIF_MO_FOCK, nmo);

        // View out our frozen orbitals
        View view_oei_ij(f, occOrbsPI, occOrbsPI, frzcpi, frzcpi);
        SharedMatrix fij = view_oei_ij();
        fij->set_name("fIJ");

        View view_oei_ab(f, virOrbsPI, virOrbsPI, frzcpi + occOrbsPI, frzcpi + occOrbsPI);
        SharedMatrix fab = view_oei_ab();
        fab->set_name("fAB");

        dpdfile2 F;
        global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, ID('O'), ID('O'), "fIJ");
        fij->write_to_dpdfile2(&F);
        global_dpd_->file2_close(&F);

        global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, ID('V'), ID('V'), "fAB");
        fab->write_to_dpdfile2(&F);
        global_dpd_->file2_close(&F);

        global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, ID('O'), ID('O'), "fij");
        fij->write_to_dpdfile2(&F);
        global_dpd_->file2_close(&F);

        global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, ID('V'), ID('V'), "fab");
        fab->write_to_dpdfile2(&F);
        global_dpd_->file2_close(&F);

        fprintf(outfile, "done\n\n");
    }

    // => Form denominators <=
    {
        fprintf(outfile, "\tForming denominators ... "); fflush(outfile);
        denom(occOrbsPI, virOrbsPI, occOffset, virOffset);
        fprintf(outfile, "done\n\n"); fflush(outfile);
    }

    // => Form initial LMP2 amplitudes <=
    {
        dpdbuf4 D;
        dpdbuf4 T2;

        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D <ij|ab>");
        global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "LMP2 tIjAb");
        global_dpd_->buf4_close(&D);

        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "LMP2 tIjAb");
        global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "dIjAb");
        global_dpd_->buf4_dirprd(&D, &T2);
        global_dpd_->buf4_close(&D);
        global_dpd_->buf4_close(&T2);
    }

    // => Compute the LMP2 energy <=
    double energy = 0.0;
    double rms = 0.0;
    int conv = 0;
    {
        dpdbuf4 D, T2, newT2;
        dpdfile2 fij, fab;

        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D 2<ij|ab> - <ij|ba>");
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "LMP2 tIjAb");
        energy = global_dpd_->buf4_dot(&D, &T2);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_close(&D);

        fprintf(outfile, "\tComputing LMP2 amplitudes:\n");
        fprintf(outfile, "\t==========================\n");
        fprintf(outfile, "\titer = %3d LMP2 Energy = %20.14lf\n", 0, energy); fflush(outfile);

        for (int iter=1; iter< maxiter; iter++) {
            global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D <ij|ab>");
            global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "New LMP2 tIjAb Increment");
            global_dpd_->buf4_close(&D);

            global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "New LMP2 tIjAb Increment");
            global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "LMP2 tIjAb");
            global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, ID('O'), ID('O'), "fIJ");
            global_dpd_->contract424(&T2, &fij, &newT2, 1, 0, 1, -1.0, 1.0);
            global_dpd_->contract244(&fij, &T2, &newT2, 0, 0, 0, -1.0, 1.0);
            global_dpd_->file2_close(&fij);

            global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, ID('V'), ID('V'), "fAB");
            global_dpd_->contract244(&fab, &T2, &newT2, 1, 2, 1, 1.0, 1.0);
            global_dpd_->contract424(&T2, &fab, &newT2, 3, 1, 0, 1.0, 1.0);
            global_dpd_->buf4_copy(&T2, PSIF_CC_TAMPS, "New LMP2 tIjAb");
            global_dpd_->buf4_close(&T2);

            global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "dIjAb");
            global_dpd_->buf4_dirprd(&D, &newT2);
            global_dpd_->buf4_close(&D);
            global_dpd_->buf4_close(&newT2);

            global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "New LMP2 tIjAb");
            global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "New LMP2 tIjAb Increment");
            global_dpd_->buf4_axpy(&T2, &newT2, 1.0);
            global_dpd_->buf4_close(&T2);

            global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D 2<ij|ab> - <ij|ba>");
            energy = global_dpd_->buf4_dot(&D, &newT2);
            global_dpd_->buf4_close(&D);

            global_dpd_->buf4_close(&newT2);

            /* Check for convergence */
            global_dpd_->buf4_init(&newT2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "New LMP2 tIjAb");
            global_dpd_->buf4_mat_irrep_init(&newT2, 0);
            global_dpd_->buf4_mat_irrep_rd(&newT2, 0);

            global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "LMP2 tIjAb");
            global_dpd_->buf4_mat_irrep_init(&T2, 0);
            global_dpd_->buf4_mat_irrep_rd(&T2, 0);

            rms = 0.0;
            for(int row=0; row < T2.params->rowtot[0]; row++)
              for(int col=0; col < T2.params->coltot[0]; col++)
                rms += (newT2.matrix[0][row][col] - T2.matrix[0][row][col]) *
                  (newT2.matrix[0][row][col] - T2.matrix[0][row][col]);

            global_dpd_->buf4_mat_irrep_close(&T2, 0);
            global_dpd_->buf4_mat_irrep_close(&newT2, 0);
            global_dpd_->buf4_close(&T2);
            global_dpd_->buf4_close(&newT2);

            rms = sqrt(rms);

            fprintf(outfile, "\titer = %3d LMP2 Energy = %20.14f   RMS = %4.3e\n", iter, energy, rms); fflush(outfile);

            if(rms < convergence) {
              conv = 1;
              fprintf(outfile, "\n\tLMP2 Iterations converged.\n");
              break;
            }
            else {
              global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "New LMP2 tIjAb");
              global_dpd_->buf4_copy(&T2, PSIF_CC_TAMPS, "LMP2 tIjAb");
              global_dpd_->buf4_close(&T2);
            }
        }

        if(!conv) {
          fprintf(outfile, "\n\tLMP2 Iterative procedure failed.\n");
          throw ConvergenceError<int>("LMP2 interative procedure failed.", maxiter, 1.0e-8, rms, __FILE__, __LINE__);
        }

        fprintf(outfile, "\tLMP2 Correlation Energy = %20.14f\n", energy);
        fprintf(outfile, "\tLMP2 Total Energy       = %20.14f\n\n", energy+Process::environment.globals["HF TOTAL ENERGY"]);

        Process::environment.globals["CURRENT_ENERGY"] = energy+Process::environment.globals["HF TOTAL ENERGY"];
        Process::environment.globals["LMP2 CORRELATION ENERGY"] = energy;

        fflush(outfile);

        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "LMP2 tIjAb");
        global_dpd_->buf4_close(&T2);

        double* epair;
        pair_energies(ints, occOrbsPI, &epair);
        print_pair_energies(occOrbsPI, epair);
    }

    // Close files we opened.
    psio_->close(PSIF_CC_OEI, 0);
    psio_->close(PSIF_CC_AINTS, 0);
    psio_->close(PSIF_CC_CINTS, 0);
    psio_->close(PSIF_CC_DINTS, 0);
    psio_->close(PSIF_CC_DENOM, 0);
    psio_->close(PSIF_CC_TAMPS, 0);

    return energy+Process::environment.globals["HF TOTAL ENERGY"];
}

PsiReturnType lmp2(Options& options)
{
    boost::shared_ptr<Wavefunction> wave(new Lmp2(Process::environment.wavefunction(), options));
    wave->compute_energy();
    Process::environment.set_wavefunction(wave);

    return Success;
}

}} // End namespaces

