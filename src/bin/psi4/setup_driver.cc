#include "psi4-dec.h"

namespace psi{
    namespace ccsort   { PsiReturnType ccsort(Options &, int argc, char *argv[]); }
    namespace ccenergy { PsiReturnType ccenergy(Options &, int argc, char *argv[]); }
    namespace cctriples { PsiReturnType cctriples(Options &, int argc, char *argv[]); }
    namespace cints    { PsiReturnType cints(Options &, int argc, char *argv[]); }
    namespace cscf     { PsiReturnType cscf(Options &, int argc, char *argv[]); }
    namespace input    { PsiReturnType input(Options &, int argc, char *argv[]); }
    namespace lmp2     { PsiReturnType lmp2(Options &, int argc, char *argv[]); }
    namespace mp2      { PsiReturnType mp2(Options &, int argc, char *argv[]); }
    namespace scf      { PsiReturnType scf(Options&, int argc, char *argv[]); }
    namespace transqt2 { PsiReturnType transqt2(Options &, int argc, char *argv[]); }
    // namespace optking  { PsiReturnType optking(Options &, int argc, char *argv[]); }
    namespace psimrcc  { PsiReturnType psimrcc(Options &, int argc, char *argv[]); }

void
setup_driver(Options &options)
{
    // The list of valid job types
    options.add_str("JOBTYPE", "SP", "SP OPT");
    // The list of valid wavefunction types
    options.add_str("WFN", "SCF", "CCSD CCSD_T MP2 MRCCSD LMP2 SCF");
    // The list of valid reference types
    options.add_str("REFERENCE", "RHF", "RHF ROHF MCSCF TCSCF UHF RKS UKS");
    // The list of valid derivative types
    options.add_str("DERTYPE", "ENERGY", "NONE ENERGY FIRST SECOND");
    options.read_ipv1();

    // make a map of function pointers to the functions
    dispatch_table["CCENERGY"]  = &(psi::ccenergy::ccenergy);
    dispatch_table["CCSORT"]    = &(psi::ccsort::ccsort);
    dispatch_table["CCTRIPLES"] = &(psi::cctriples::cctriples);
    dispatch_table["CINTS"]     = &(psi::cints::cints);
    dispatch_table["CSCF"]      = &(psi::cscf::cscf);
    dispatch_table["INPUT"]     = &(psi::input::input);
    dispatch_table["LMP2"]      = &(psi::lmp2::lmp2);
    dispatch_table["MP2"]       = &(psi::mp2::mp2);
    dispatch_table["SCF"]       = &(psi::scf::scf);
    dispatch_table["TRANSQT2"]  = &(psi::transqt2::transqt2);
    // dispatch_table["OPTKING"]   = &(psi::optking::optking);
    dispatch_table["PSIMRCC"]   = &(psi::psimrcc::psimrcc);
}

} // Namespaces
