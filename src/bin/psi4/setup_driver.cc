#include "psi4.h"
#include "psi4-dec.h"

namespace psi{
    namespace ccsort   { PsiReturnType ccsort(Options &); }
    namespace ccenergy { PsiReturnType ccenergy(Options &); }
    namespace cctriples{ PsiReturnType cctriples(Options &); }
    namespace cints    { PsiReturnType cints(Options &); }
    namespace cscf     { PsiReturnType cscf(Options &); }
    namespace dcft     { PsiReturnType dcft(Options &); }
    namespace input    { PsiReturnType input(Options &); }
    namespace lmp2     { PsiReturnType lmp2(Options &); }
    namespace mp2      { PsiReturnType mp2(Options &); }
    namespace scf      { PsiReturnType scf(Options &); }
    namespace dfmp2    { PsiReturnType dfmp2(Options &); }
    namespace transqt2 { PsiReturnType transqt2(Options &); }
    namespace psimrcc  { PsiReturnType psimrcc(Options &); }
    namespace mcscf    { PsiReturnType mcscf(Options &); }

    PsiReturnType opt_step(Options &);

void
setup_driver(Options &options)
{
    // The list of valid job types
    options.add_str("JOBTYPE", "SP", "SP OPT");
    // The list of valid wavefunction types
    options.add_str("WFN", "SCF", "DCFT CCSD CCSD_T MP2 MRCCSD LMP2 SCF MCSCF BCCD BCCD_T OOCCD");
    // The list of valid reference types
    options.add_str("REFERENCE", "RHF", "RHF ROHF MCSCF TCSCF UHF RKS UKS");
    // The list of valid derivative types
    options.add_str("DERTYPE", "ENERGY", "NONE ENERGY FIRST SECOND RESPONSE");
    options.read_ipv1();

    // make a map of function pointers to the functions
    dispatch_table["CCENERGY"]  = &(psi::ccenergy::ccenergy);
    dispatch_table["CCSORT"]    = &(psi::ccsort::ccsort);
    dispatch_table["CCTRIPLES"] = &(psi::cctriples::cctriples);
    dispatch_table["CINTS"]     = &(psi::cints::cints);
    dispatch_table["CSCF"]      = &(psi::cscf::cscf);
    dispatch_table["DCFT"]      = &(psi::dcft::dcft);
    dispatch_table["INPUT"]     = &(psi::input::input);
    dispatch_table["LMP2"]      = &(psi::lmp2::lmp2);
    dispatch_table["MP2"]       = &(psi::mp2::mp2);
    dispatch_table["SCF"]       = &(psi::scf::scf);
    dispatch_table["DFMP2"]     = &(psi::dfmp2::dfmp2);
    dispatch_table["TRANSQT2"]  = &(psi::transqt2::transqt2);
    dispatch_table["MCSCF"]     = &(psi::mcscf::mcscf);
    dispatch_table["PSIMRCC"]   = &(psi::psimrcc::psimrcc);
    dispatch_table["OPT_STEP"]  = &(psi::opt_step);
}

} // Namespaces
