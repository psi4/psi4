#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace efp {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "EFP"|| options.read_globals()) {
        /*- The amount of information printed to the output file. -*/
        options.add_int("PRINT", 1);
        /*- Type of EFP simulation. One of single point (``SP``), gradient
        (``GRAD``), conjugate gradient geometry optimization (``CG``),
        molecular dynamics in microcanonical ensemble (``NVE``), or
        molecular dynamics in canonical ensemble (``NVT``). This
        specification will probably be moved to energy(), grad(), opt(),
        etc. eventually. -*/
        options.add_str("EFP_TYPE", "SP", "SP GRAD CG NVE NVT");
        /*- Do include electrostatics energy term in EFP computation? -*/
        options.add_bool("EFP_ELST", true);
        /*- Do include polarization energy term in EFP computation? -*/
        options.add_bool("EFP_POL", true);
        /*- Do include dispersion energy term in EFP computation? -*/
        options.add_bool("EFP_DISP", true);
        /*- Do include exchange repulsion energy term in EFP computation? -*/
        options.add_bool("EFP_EXCH", true);
        /*- Electrostatic damping type. ``SCREEN`` is a damping formula
        based on screen group in the EFP potential. ``OVERLAP`` is
        damping that computes charge penetration energy. -*/
        options.add_str("EFP_ELST_DAMPING", "SCREEN", "SCREEN OVERLAP OFF");
        /*- Dispersion damping type. ``TT`` is a damping formula by
        Tang and Toennies. ``OVERLAP`` is overlap-based dispersion damping. -*/
        options.add_str("EFP_DISP_DAMPING", "OVERLAP", "TT OVERLAP OFF");
    }

    return true;
}

extern "C" 
PsiReturnType efp(Options& options)
{
    int print = options.get_int("PRINT");

    bool elst_enabled = options.get_bool("EFP_ELST");
    bool pol_enabled  = options.get_bool("EFP_POL");
    bool disp_enabled = options.get_bool("EFP_DISP");
    bool exch_enabled = options.get_bool("EFP_EXCH");

    std::string efp_mode = options.get_str("EFP_TYPE");
    std::string elst_damping = options.get_str("EFP_ELST_DAMPING");
    std::string disp_damping = options.get_str("EFP_DISP_DAMPING");

    if (elst_damping == "SCREEN") {
    } else if (elst_damping == "OVERLAP") {
    } else if (elst_damping == "OFF") {
    }

    if (disp_damping == "TT") {
    } else if (disp_damping == "OVERLAP") {
    } else if (disp_damping == "OFF") {
    }

    if (efp_mode == "SP") {
    } else if (efp_mode == "GRAD") {
    } else if (efp_mode == "CG") {
    } else if (efp_mode == "NVE") {
    } else if (efp_mode == "NVT") {
    }

    if (elst_enabled);
    if (pol_enabled);
    if (disp_enabled);
    if (exch_enabled);

    double e_elst = 0.0;
    double e_pol = 0.0;
    double e_disp = 0.0;
    double e_exch = 0.0;
    double e_total = 0.0;

    fprintf(outfile, "\n");
    fprintf(outfile, "         --------------------------------------------------------\n");
    fprintf(outfile, "                                   EFP                           \n");
    fprintf(outfile, "                   Effective Fragment Potential Method           \n");
    fprintf(outfile, "\n");
    fprintf(outfile, "                              Ilya Kaliman                       \n");
    fprintf(outfile, "         --------------------------------------------------------\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  ==> Calculation Information <==\n\n");

    //molecule_->print();

    fprintf(outfile, "  Simulation type:        %12s\n", efp_mode.c_str());
    fprintf(outfile, "  Electrostatics damping: %12s\n", elst_damping.c_str());
    fprintf(outfile, "  Dispersion damping:     %12s\n", disp_damping.c_str());
    fprintf(outfile, "\n");

    e_total = e_elst + e_pol + e_disp + e_exch;

    fprintf(outfile, "  ==> Energetics <==\n\n");

    fprintf(outfile, "  Electrostatics Energy = %24.16f [H] %s\n", e_elst, elst_enabled ? "*" : "");
    fprintf(outfile, "  Polarization Energy =   %24.16f [H] %s\n", e_pol, pol_enabled ? "*" : "");
    fprintf(outfile, "  Dispersion Energy =     %24.16f [H] %s\n", e_disp, disp_enabled ? "*" : "");
    fprintf(outfile, "  Exchange Energy =       %24.16f [H] %s\n", e_exch, exch_enabled ? "*" : "");
    fprintf(outfile, "  Total Energy =          %24.16f [H] %s\n", e_total, "*");

    Process::environment.globals["EFP ELST ENERGY"] = e_elst;
    Process::environment.globals["EFP POL ENERGY"] = e_pol;
    Process::environment.globals["EFP DISP ENERGY"] = e_disp;
    Process::environment.globals["EFP EXCH ENERGY"] = e_exch;
    Process::environment.globals["CURRENT ENERGY"] = e_total;


    return Success;
}

}} // End namespaces

