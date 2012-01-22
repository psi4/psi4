#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/view.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>
#include <libdpd/dpd.h>
#include <libciomr/libciomr.h>
#include <vector>

#include <physconst.h>

#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/dict.hpp>

namespace psi{ namespace mrcc {

void write_oei_to_disk(FILE* fort55, SharedMatrix moH)
{
    // Walk through moH and save the non-zero values
    int offset = 0;
    for (int h=0; h<moH->nirrep(); ++h) {
        for (int m=0; m<moH->rowdim(h); ++m) {
            for (int n=0; n<=m; ++n) {
                if (fabs(moH->get(h, m, n)) > 1.0e-12) {
                    fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", moH->get(h, m, n), m+offset+1, n+offset+1, 0, 0);
                }
            }
        }
        offset += moH->rowdim(h);
    }
}

void write_tei_to_disk(FILE* fort55, int nirrep, dpdbuf4& K, double ints_tolerance)
{
    for(int h = 0; h < nirrep; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            int p = K.params->roworb[h][pq][0];
            int q = K.params->roworb[h][pq][1];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                int r = K.params->colorb[h][rs][0];
                int s = K.params->colorb[h][rs][1];

                if (fabs(K.matrix[h][pq][rs]) > ints_tolerance)
                    fprintf(fort55, "%28.20E%4d%4d%4d%4d\n",
                            K.matrix[h][pq][rs], p+1, q+1, r+1, s+1);
            }
        }
        dpd_buf4_mat_irrep_close(&K, h);
    }
}

void print_dim(const std::string& name, const Dimension& dim)
{
    fprintf(outfile, "        %15s [ ", name.c_str());
    for (int h=0; h<dim.n(); ++h) {
        fprintf(outfile, "%3d", dim[h]);
        if (h != dim.n()-1)
            fprintf(outfile, ",");
    }
    fprintf(outfile, "]\n");
}

PsiReturnType mrcc(Options& options, const boost::python::dict& level)
{
    tstart();

    fprintf(outfile, "  PSI4 interface to MRCC:\n");

    // Ensure the dict provided has everything we need.
    if (!level.has_key("method") ||
        !level.has_key("order") ||
        !level.has_key("fullname"))
        throw PSIEXCEPTION("MRCC interface: Provided dictionary is incomplete.");

    int method  = boost::python::extract<int>(level.get("method"));
    int exlevel = boost::python::extract<int>(level.get("order"));
    std::string fullname = boost::python::extract<std::string>(level.get("fullname"));
    bool pertcc = exlevel > 0 ? false : true;
    exlevel = abs(exlevel);

    fprintf(outfile, "    Generating inputs for %s.\n\n", fullname.c_str());

    fprintf(outfile, "    Automatically determined settings:\n");
    fprintf(outfile, "        method %d\n        exlevel %d\n        fullname %s\n\n",
            method, exlevel, fullname.c_str());

    boost::shared_ptr<Wavefunction> wave     = Process::environment.reference_wavefunction();
    boost::shared_ptr<Molecule>     molecule = wave->molecule();

    // Orbitals spaces
    Dimension docc        = wave->doccpi();
    Dimension frzcpi      = wave->frzcpi();
    Dimension frzvpi      = wave->frzvpi();
    Dimension active_docc = docc - frzcpi;
    Dimension active_socc = wave->soccpi();
    Dimension active_mopi = wave->nmopi() - frzcpi - frzvpi;

    int nbf = active_mopi.sum();
    int nirrep = wave->nirrep();
    int nelectron = 2 * active_docc.sum() + active_socc.sum();

    fprintf(outfile, "    Orbital Information:\n\n");
    print_dim("Frozen Core", frzcpi);
    print_dim("Active DOCC", active_docc);
    print_dim("SOCC", active_socc);
    print_dim("Frozen Virtual", frzvpi);

    fprintf(outfile, "\n");
    print_dim("Total MOs", active_mopi);

    fprintf(outfile, "\n");
    FILE* fort55 = fopen("fort.55", "w");
    fprintf(fort55, "%22d%22d\n", nbf, nelectron);

    // Print out orbital symmetries
    int count=0;
    for (int h=0; h<active_mopi.n(); ++h) {
        for (int n=0; n<active_mopi[h]; ++n) {
            fprintf(fort55, "%22d", h+1);  // 1 based irrep ordering

            count++;
            if (count % 3 == 0) // We only want 3 per line
                fprintf(fort55, "\n");
        }
    }

    if (count % 3) // integrals start on their own line
        fprintf(fort55, "\n");

    fprintf(outfile, "    Beginning integral transformation.\n");

    // Define the orbital space of the MO integrals we need.
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);

    // Check the reference.
    bool restricted = true;

    if (options.get_str("REFERENCE") == "UHF")
        restricted = false;

    if (pertcc && options.get_str("REFERENCE") == "ROHF")
        throw PSIEXCEPTION("Perturbative methods not implemented for ROHF references.");

    // Create integral transformation object
    IntegralTransform ints(wave, spaces, restricted ? IntegralTransform::Restricted : IntegralTransform::Unrestricted);

    // This transforms everything (OEI and TEI)
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);

    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());

    fprintf(outfile, "    Transformation complete.\n\n");
    fprintf(outfile, "  Generating fort.55 integral file...");

    double ints_tolerance = options.get_double("INTS_TOLERANCE");

    _default_psio_lib_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    dpdbuf4 K;

    fprintf(fort55, " 150000\n");

    // RHF
    if (restricted) {

        // We want only the permutationally unique integrals, hence [A>=A]+, see libtrans documenation for details
        dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0,
                      ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), // In memory
                      ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), // On disk
                      0,
                      "MO Ints (AA|AA)");
        write_tei_to_disk(fort55, nirrep, K, ints_tolerance);
        dpd_buf4_close(&K);

        // Load in frozen core operator, in the event of FREEZE_CORE = FALSE this is the MO OEI
        SharedMatrix moH(new Matrix(PSIF_MO_FZC, wave->nmopi(), wave->nmopi()));
        moH->load(_default_psio_lib_, PSIF_OEI);
        View vmoH(moH, active_mopi, active_mopi, wave->frzcpi(), wave->frzcpi());
        moH = vmoH();
        write_oei_to_disk(fort55, moH);

        // Print nuclear repulsion energy.
        // Eventually needs to be changed to frozen core energy + nuclear repulsion energy
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", ints.get_frozen_core_energy() + molecule->nuclear_repulsion_energy(), 0, 0, 0, 0);
    }
    else {

        // We want only the permutationally unique integrals, hence [A>=A]+, see libtrans documenation for details

        // Load up alpha alpha
        dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0,
                      ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), // In memory
                      ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), // On disk
                      0,
                      "MO Ints (AA|AA)");
        write_tei_to_disk(fort55, nirrep, K, ints_tolerance);
        dpd_buf4_close(&K);

        // Write out separator
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", 0.0, 0, 0, 0, 0);

        // Load up beta beta
        dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0,
                      ints.DPD_ID("[a>=a]+"), ints.DPD_ID("[a>=a]+"), // In memory
                      ints.DPD_ID("[a>=a]+"), ints.DPD_ID("[a>=a]+"), // On disk
                      0,
                      "MO Ints (aa|aa)");
        write_tei_to_disk(fort55, nirrep, K, ints_tolerance);
        dpd_buf4_close(&K);

        // Write out separator
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", 0.0, 0, 0, 0, 0);

        // Load up alpha beta
        dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0,
                      ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[a>=a]+"), // In memory
                      ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[a>=a]+"), // On disk
                      0,
                      "MO Ints (AA|aa)");
        write_tei_to_disk(fort55, nirrep, K, ints_tolerance);
        dpd_buf4_close(&K);

        // Write out separator
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", 0.0, 0, 0, 0, 0);

        // Load in alpha frozen core operator, in the event of FREEZE_CORE = FALSE this is the MO OEI
        SharedMatrix moH(new Matrix(PSIF_MO_A_FZC, wave->nmopi(), wave->nmopi()));
        moH->load(_default_psio_lib_, PSIF_OEI);
        View vmoH(moH, active_mopi, active_mopi, wave->frzcpi(), wave->frzcpi());
        moH = vmoH();
        write_oei_to_disk(fort55, moH);

        // Write out separator
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", 0.0, 0, 0, 0, 0);

        // Load in beta frozen core operator, in the event of FREEZE_CORE = FALSE this is the MO OEI
        SharedMatrix moHb(new Matrix(PSIF_MO_B_FZC, wave->nmopi(), wave->nmopi()));
        moHb->load(_default_psio_lib_, PSIF_OEI);
        View vmoHb(moHb, active_mopi, active_mopi, wave->frzcpi(), wave->frzcpi());
        moHb = vmoHb();
        write_oei_to_disk(fort55, moHb);

        // Write out separator
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", 0.0, 0, 0, 0, 0);

        // Print nuclear repulsion energy.
        // Eventually needs to be changed to frozen core energy + nuclear repulsion energy
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", ints.get_frozen_core_energy() + molecule->nuclear_repulsion_energy(), 0, 0, 0, 0);
    }
    _default_psio_lib_->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    fclose(fort55);

    fprintf(outfile, "done.\n  Generating fort.56 input file...");

    // Determine energy convergence to pass to MRCC
    double user_e = fabs(Process::environment.options.get_double("E_CONVERGENCE"));
    int e_conv = 0;

    if (user_e >= 1.0)
        throw InputException("Your energy convergence criterea is absurb.", "E_CONVERGENCE", user_e, __FILE__, __LINE__);

    while (user_e <= 0.1) {
        user_e *= 10.0;
        e_conv++;
    }

    if (e_conv == 0)
        throw PSIEXCEPTION("e_conv became 0. This isn't right.");

    // The following two are expert options.
    if (options["MRCC_METHOD"].has_changed())
        method = options.get_int("MRCC_METHOD");
    if (options["MRCC_LEVEL"].has_changed())
        exlevel = options.get_int("MRCC_LEVEL");

    int nsing = 1;
    int closed_shell = 1;
    int spatial_orbitals = 1;
    int ndoub = 0;
    if (pertcc && (options.get_str("REFERENCE") == "ROHF" || options.get_str("REFERENCE") == "UHF")) {
        nsing = 0;
        closed_shell = 0;
        spatial_orbitals = 0;
        ndoub = 1;
    }
    if (options["MRCC_NUM_SINGLET_ROOTS"].has_changed())
        nsing = options.get_int("MRCC_NUM_SINGLET_ROOTS");

    int symm = 0;
    for (int h=0; h<nirrep; ++h)
        for (int n=0; n<active_socc[h]; ++n)
            symm ^= h;
    symm += 1; // stupid 1 based fortran

    FILE* fort56 = fopen("fort.56", "w");
    fprintf(fort56, "%6d%6d%6d%6d%6d      0     0%6d     0%6d%6d     1%6d      0      0%6d     0     0    0.00    0%6lu\n",
            exlevel,                                         // # 1
            nsing,                                           // # 2
            options.get_int("MRCC_NUM_TRIPLET_ROOTS"),       // # 3
            options.get_int("MRCC_RESTART"),                 // # 4
            method,                                          // # 5
            symm,                                            // # 8
            closed_shell,                                    // #10
            spatial_orbitals,                                // #11
            ndoub,                                           // #13
            e_conv,                                          // #16
            Process::environment.get_memory() / 1000 / 1000  // #21
        );

    fprintf(fort56, "ex.lev,nsing,ntrip, rest,CC/CI, dens,conver, symm, diag,   CS,spatial, HF,ndoub, nacto, nactv,  tol,maxex, sacc,   freq,dboc,  mem\n");

    for (int h=0; h < nirrep; ++h) {
        for (int nmo = 0; nmo < active_mopi[h]; ++nmo) {
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

    fprintf(outfile, "done.\n");
    fflush(outfile);

    tstop();

    return Success;
}

}} // End Namespaces

