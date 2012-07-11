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
        /*- Names of fragment files corresponding to molecule subsets.
        This is temporary until better EFP input geometry parsing is implemented. -*/
        options.add("FRAGS", new ArrayType());
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

    int nfrag = options["FRAGS"].size();
    std::vector<std::string> frag_name;

    for(int i=0; i<nfrag; i++)
        frag_name.push_back(options["FRAGS"][i].to_string());


    double e_elst = 0.0;
    double e_pol = 0.0;
    double e_disp = 0.0;
    double e_exch = 0.0;
    double e_total = 0.0;

    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

    fprintf(outfile, "\n");
    fprintf(outfile, "         --------------------------------------------------------\n");
    fprintf(outfile, "                                   EFP                           \n");
    fprintf(outfile, "                   Effective Fragment Potential Method           \n");
    fprintf(outfile, "\n");
    fprintf(outfile, "                              Ilya Kaliman                       \n");
    fprintf(outfile, "         --------------------------------------------------------\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  ==> Geometry <==\n\n");

    molecule->print();
    if (molecule->nfragments() != nfrag)
        throw InputException("Molecule doesn't have FRAGS number of fragments.", "FRAGS", nfrag, __FILE__, __LINE__);

    std::vector<int> realsA;
    realsA.push_back(0);
    std::vector<int> ghostsA;
    boost::shared_ptr<Molecule> monomerA = molecule->extract_subsets(realsA, ghostsA);
    monomerA->print();
    monomerA->print_in_bohr();

    std::vector<int> realsB;
    realsB.push_back(1);
    std::vector<int> ghostsB;
    ghostsB.push_back(0);
    boost::shared_ptr<Molecule> monomerB = molecule->extract_subsets(realsB, ghostsB);
    //monomerB->print();

    std::vector<int> realsC;
    realsC.push_back(2);
    std::vector<int> ghostsC;
    ghostsC.push_back(0);
    ghostsC.push_back(1);
    boost::shared_ptr<Molecule> monomerC = molecule->extract_subsets(realsC, ghostsC);
    //monomerC->print();

    int natomA = 0, natomB = 0, natomC = 0;
    for (int n=0; n<monomerA->natom(); n++)
      if (monomerA->Z(n)) natomA++;
    for (int n=0; n<monomerB->natom(); n++)
      if (monomerB->Z(n)) natomB++;
    for (int n=0; n<monomerC->natom(); n++)
      if (monomerC->Z(n)) natomC++;

    if (natomA != 3)
        throw InputException("Fragment doesn't have three coordinate triples.", "natomA", natomA, __FILE__, __LINE__);
    if (natomB != 3)
        throw InputException("Fragment doesn't have three coordinate triples.", "natomB", natomB, __FILE__, __LINE__);
    if (natomC != 3)
        throw InputException("Fragment doesn't have three coordinate triples.", "natomC", natomC, __FILE__, __LINE__);

    SharedMatrix xyz = SharedMatrix (new Matrix("Fragment Cartesian Coordinates(x,y,z)", monomerA->natom(), 3));
    double** xyzp = xyz->pointer();

    for (int A = 0; A < monomerA->natom(); A++) {
        xyzp[A][0] = monomerA->x(A);
        xyzp[A][1] = monomerA->y(A);
        xyzp[A][2] = monomerA->z(A);
    }

    for(int i=0; i<nfrag; i++)
        fprintf(outfile, "%s\n", frag_name[i].c_str());

    xyz->print();


    fprintf(outfile, "  ==> Calculation Information <==\n\n");

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

