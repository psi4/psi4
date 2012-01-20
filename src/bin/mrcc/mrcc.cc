#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/view.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>
#include <libdpd/dpd.h>
#include <vector>

#include <physconst.h>

#include <boost/foreach.hpp>

namespace psi{ namespace mrcc {

PsiReturnType mrcc(int level, bool pertcc)
{
    fprintf(outfile, "\n  MRCC input file generator.\n");
    fprintf(outfile, "   Will generate fort.55 and fort.56 files.\n\n");

    boost::shared_ptr<Wavefunction> wave     = Process::environment.reference_wavefunction();
    boost::shared_ptr<Molecule>     molecule = wave->molecule();
    boost::shared_ptr<BasisSet>     aobasis  = wave->basisset();
    boost::shared_ptr<SOBasisSet>   sobasis  = wave->sobasisset();

    int nbf = wave->nmopi().sum() - wave->frzcpi().sum() - wave->frzvpi().sum();
    int nirrep = wave->nirrep();
    const Dimension& docc = wave->doccpi();
    int nelectron = 2 * (docc.sum() - wave->frzcpi().sum());

    FILE* fort55 = fopen("fort.55", "w");
    fprintf(fort55, "%22d%22d\n", nbf, nelectron);

    // Print out orbital symmetries
    Dimension mopi = wave->nmopi() - wave->frzcpi() - wave->frzvpi();
    int count=0;
    for (int h=0; h<mopi.n(); ++h) {
        for (int n=0; n<mopi[h]; ++n) {
            fprintf(fort55, "%22d", h+1);  // 1 based irrep ordering

            count++;
            if (count % 3 == 0) // We only want 3 per line
                fprintf(fort55, "\n");
        }
    }

    if (count % 3) // integrals start on their own line
        fprintf(fort55, "\n");

    // Need to figure out this. I think it means closed shell
    fprintf(fort55, " 150000\n");

    fprintf(outfile, "    Beginning integral transformation.\n");
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(wave, spaces, IntegralTransform::Restricted);
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());

    fprintf(outfile, "    Transformation complete.\n\n");
    fprintf(outfile, "  Generating fort.55 integral file...");

    /*
     * Now, loop over the DPD buffer, printing the integrals
     */
    dpdbuf4 K;
    _default_psio_lib_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // We want only the permutationally unique integrals, hence [A>=A]+
    dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0,
                  ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), // In memory
                  ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), // On disk
                  0,
                  "MO Ints (AA|AA)");
    for(int h = 0; h < nirrep; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        // Need to restrict these loops
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            int p = K.params->roworb[h][pq][0];
            int q = K.params->roworb[h][pq][1];
            int psym = K.params->psym[p];
            int qsym = K.params->qsym[q];
            int prel = p - K.params->poff[psym];
            int qrel = q - K.params->qoff[qsym];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                int r = K.params->colorb[h][rs][0];
                int s = K.params->colorb[h][rs][1];
                int rsym = K.params->rsym[r];
                int ssym = K.params->ssym[s];
                int rrel = r - K.params->roff[rsym];
                int srel = s - K.params->soff[ssym];
                // Print out the absolute orbital numbers, the relative (within irrep)
                // numbers, the symmetries, and the integral itself
                if (fabs(K.matrix[h][pq][rs]) > 1.0e-12)
                    fprintf(fort55, "%28.20E%4d%4d%4d%4d\n",
                            K.matrix[h][pq][rs], p+1, q+1, r+1, s+1);
//                fprintf(outfile, "(%2d %2d | %2d %2d) = %16.10f, "
//                                 "symmetries = (%1d %1d | %1d %1d), "
//                                 "relative indices = (%2d %2d | %2d %2d)\n",
//                                 p, q, r, s, K.matrix[h][pq][rs],
//                                 psym, qsym, rsym, ssym,
//                                 prel, qrel, rrel, srel);
            }
        }
        dpd_buf4_mat_irrep_close(&K, h);
    }
    dpd_buf4_close(&K);
    _default_psio_lib_->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // Load in frozen core operator, in the event of FREEZE_CORE = FALSE this is the MO OEI
    SharedMatrix moH(new Matrix(PSIF_MO_FZC, wave->nmopi(), wave->nmopi()));
    moH->load(_default_psio_lib_, PSIF_OEI);

    View vmoH(moH, mopi, mopi, wave->frzcpi(), wave->frzcpi());
    moH = vmoH();

    // Walk through moH and save the non-zero values
    int offset = 0;
    for (int h=0; h<nirrep; ++h) {
        for (int m=0; m<mopi[h]; ++m) {
            for (int n=0; n<=m; ++n) {
                if (fabs(moH->get(h, m, n)) > 1.0e-12) {
                    fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", moH->get(h, m, n), m+offset+1, n+offset+1, 0, 0);
                }
            }
        }
        offset += mopi[h];
    }

    // Print nuclear repulsion energy.
    // Eventually needs to be changed to frozen core energy + nuclear repulsion energy
    fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", ints.get_frozen_core_energy() + molecule->nuclear_repulsion_energy(), 0, 0, 0, 0);
    fclose(fort55);

    fprintf(outfile, "done.\n  Generating fort.56 input file...");

    FILE* fort56 = fopen("fort.56", "w");
    fprintf(fort56, "%6d     1     0     0%7d      0     0      1     0     1     1     1     0      0      0     8     0     0    0.00    0%6lu\n",
            level, pertcc ? 3 : 1, Process::environment.get_memory() / 1000 / 1000);
    fprintf(fort56, "ex.lev,nsing,ntrip, rest, CC/CI,  dens,conver, symm, diag,   CS,spatial, HF,ndoub, nacto, nactv,  tol,maxex, sacc,   freq,dboc,  mem\n");

    const Dimension& active_docc = docc - wave->frzcpi();
    const Dimension& active_socc = wave->soccpi();

    for (int h=0; h < nirrep; ++h) {
        for (int nmo = 0; nmo < mopi[h]; ++nmo) {
            int val = 0;
            if (nmo < active_docc[h])
                val = 2;
            else if (nmo < (active_docc[h] + active_socc[h]))
                val = 1;

            fprintf(fort56, "%3d", val);
        }
    }
    fprintf(fort56, "\n");
    fclose(fort56);

    fprintf(outfile, "done.\n\n");
    fflush(outfile);

    return Success;
}

}} // End Namespaces

