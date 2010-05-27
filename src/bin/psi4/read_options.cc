/*! \file read_calculation_options
    \defgroup PSI4
*/

#include <libipv1/ip_lib.h>
#include <liboptions/liboptions.h>
#include <libparallel/parallel.h>
#include <physconst.h>
#include <psifiles.h>

namespace psi {

/**
 * This is called immediately before a module is run.  Any options
 * expected by that module must be added here
 *
 * @param name    - the name of the module.
 * @param options - the liboptions module used in the computations.
 */

int read_options(const std::string &name, Options & options) {

  ip_cwk_clear();
  ip_cwk_add(const_cast<char*>(":BASIS"));
  ip_cwk_add(const_cast<char*>(":DEFAULT"));
  ip_cwk_add(const_cast<char*>(":PSI"));
  ip_set_uppercase(1);

  options.clear();
  if (name == "PSI4") {
    options.add_str("UNITS", "ANGSTROMS", "BOHR AU ANGSTROMS ANGSTROM");
  }
  else if (name == "INPUT") {
    ip_cwk_add(const_cast<char*>(":INPUT"));
    /*- The units used for the geometry -*/
    options.add_str("UNITS", "ANGSTROMS", "BOHR AU ANGSTROMS ANGSTROM");
    /*- Whether to keep the checkpoint file. -*/
    options.add_bool("KEEP_CHKPT",false);
    /*- Read MOs from the checkpoint file and project onto new basis -*/
    options.add_bool("CHKPT_MOS",false);
  // these may be obseleted in psi4
    /*- Read geometry from checkpoint file (in findif calculations) -*/
    options.add_bool("CHKPT_GEOM",false);
    /*- Don't project MOs but simply keep them -*/
    options.add_bool("NOPROJECT",false);
    /*- Read geometry from geom.dat file (in findif calculations) -*/
    options.add_bool("GEOMDAT",false);
    /*- No center of mass shift -*/
    options.add_bool("NO_COM_SHIFT",false);
    /*- Read MOs from checkpoint file and save to a separate file -*/
    options.add_bool("SAVE_MOS",false);
    /*- Don't overwrite the output file. -*/
    options.add_bool("KEEP_OUTPUT",false);
    /*- The wavefunction desired -*/
    options.add_str("WFN", "");
    /*- Whether to skip the orientation of the molecule to make its principle
        axes cooincide with the Cartesian axes -*/
    options.add_bool("NO_REORIENT",false);
    /*- The label used when storing job information -*/
    options.add_str("LABEL","Default PSI3 Label");
    options.add_bool("SHOWNORM",false);
    options.add_bool("NORMALIZE",true);
    /*- Whether pure (spherical) angular momentum functions are to be used -*/
    options.add_bool("PUREAM",false);
    /*- Override some safety checks -*/
    options.add_bool("EXPERT",false);
    /*- The amount of information to print to the output -*/
    options.add_int("PRINT",1);
    /*- The Schoenflies symbol for the computational point group to be used -*/
    options.add_str("SUBGROUP", "", "C1 C2 CS CI C2V C2H D2");
    /*- The axis to treat as the main Cn axis in the point group used, if not
        unambiguously determined -*/
    options.add_str("UNIQUE_AXIS", "", "X Y Z");
    /*- The number of fragments to consider in the inter-fragment optimization code -*/
    options.add_int("NFRAGMENTS",1);
    options.add_bool("KEEP_REF_FRAME",false);
    options.add_bool("FRAGMENT_DISTANCE_INVERSE",false);
    // input also does ip calls to read basis, geometry
    /*- How many of the innermost orbitals to keep doubly-occupied in correlated computations -*/
    options.add_str("FREEZE_CORE","FALSE", \
      "FALSE NO TRUE YES SMALL LARGE");
    /* The number of virtual orbitals to freeze in correlated computations -*/
    options.add_int("FREEZE_VIRT",0);
  }
  else if (name == "CINTS") {
    ip_cwk_add(const_cast<char*>(":CINTS"));
    /*- The wavefunction desired -*/
    options.add_str("WFN", "");
    /*- The amount of information to pring to the output -*/
    options.add_int("PRINT",1);
    /*- The number of CPU cores to use on each node -*/
    options.add_int("NUM_THREADS",1);
    /*- -Log10 of the threshold below which integrals are assumed to be zero -*/
    options.add_int("CUTOFF",15); // cutoff on integrals
    /*- The fine-structure constant -*/
    options.add_double("FINE_STRUCTURE_ALPHA", 1.0/(_c_au));
    /*- Whether to compute two electron repusion integrals -*/
    options.add_bool("MAKE_ERI", 1);
    /*- Whether empirical dispersion terms will be used -*/
    options.add_bool("EMPIRICAL_DISPERSION", 0);
    /*- Whether this computation is a restart from a previous run -*/
    options.add_bool("RESTART", 0);
    options.add_int("RESTART_TASK", 0);
    /*- The file to which the SO overlap integral will be written -*/
    options.add_int("S_FILE", PSIF_OEI);
    /*- The file to which the SO kinetic energy integrals will be written -*/
    options.add_int("T_FILE", PSIF_OEI);
    /*- The file to which the SO potential energy integrals will be written -*/
    options.add_int("V_FILE", PSIF_OEI);
    /*- The file to which the SO two electron repulsion integrals will be written -*/
    options.add_int("ERI_FILE", PSIF_SO_TEI);

  }
  else if (name == "SCF") {
    /*- The DFT grid specification, such as SG1 -*/
    options.add_str("GRID_STRING","","SG1");
    /*- The number of radial points in the DFT grid -*/
    options.add_int("N_RADIAL",99);
    /*- The number of spherical points in the DFT grid -*/
    options.add_int("N_SPHERICAL",590);  
    /*- The number of grid points per evaluation block -*/
    options.add_int("N_BLOCK",5000);  
    /*- The spherical quadrature type for DFT, usually Lebedev-*/
    options.add_str("SPHERICAL_TYPE","LEBEDEV","LEBEDEV");
    /*- The radial quadrature type for DFT, Treutler is best -*/
    options.add_str("RADIAL_TYPE","TREUTLER","TREUTLER BECKE EULER-MACLAURIN");
    /*- The nuclear partition type for DFT, Treutler is best  -*/
    options.add_str("NUCLEAR_TYPE","TREUTLER","NAIVE BECKE TREUTLER");
    /*- The exchange functional name or alias  -*/
    options.add_str("X_FUNCTIONAL","");
    /*- The correlation functional name  -*/
    options.add_str("C_FUNCTIONAL","");
    /*- Print Functional Test Data?  -*/
    options.add_bool("TEST_FUNCTIONAL",false);
     
    /*- Save a grid or not?  -*/
    options.add_bool("SAVE_CARTESIAN_GRID",false);
    /*- Grid filename  -*/
    options.add_str("CARTESIAN_FILENAME","Grid.out");
    /*- Grid global resolution (npoints)  -*/
    options.add_int("CARTESIAN_RESOLUTION",30);
    /*- Grid x resolution (npoints)  -*/
    options.add_int("CARTESIAN_RESOLUTION_X",1);
    /*- Grid y resolution (npoints)  -*/
    options.add_int("CARTESIAN_RESOLUTION_Y",1);
    /*- Grid z resolution (npoints)  -*/
    options.add_int("CARTESIAN_RESOLUTION_Z",1);
    /*- Grid MO indices  -*/
    options.add("CARTESIAN_MO_INDICES",new ArrayType());
    /*- Number MO indices  -*/
    options.add_int("N_CARTESIAN_MOS",0);

    
    /*- Save a grid or not?  -*/
    options.add_bool("SAVE_NUMERICAL_GRID",false);
    /*- Grid filename  -*/
    options.add_str("NUMERICAL_GRID_FILENAME","ngrid.out");

    /*- Grid overage at sides  -*/
    options.add_double("CARTESIAN_OVERAGE",2.00);
    /*- Grid extents [xl, xr, yl, yr, zl, zr]  -*/
    options.add("CARTESIAN_EXTENTS",new ArrayType());
    /*- The name of the auxiliary basis to be used in RI computations -*/
    options.add_str("RI_BASIS_SCF", "");
    options.read_ipv1();

    /*- Whether to use density fitting or not -*/
    if(options.get_str("RI_BASIS_SCF").empty())
      options.add_bool("RI_HF", false);
    else
      options.add_bool("RI_HF", true);

    /*- Use Local fitting for K or not? -*/
    options.add_bool("L_HF",false);
    /*- Atomic Charge cutoff (for primary domain) -*/
    options.add_double("CHARGE_CUTOFF",0.05);
    /*- Extended domain radius, Angstrom -*/
    options.add_double("R_EXT",3.0);
    /*- Iterations per full Pipek-Mizey Localization -*/
    options.add_int("STEPS_PER_LOCALIZE",1);
    

    /*- The guess type to be used in the computation -*/
    options.add_str("GUESS", "", "CORE GWH SAD READ BASIS2");
    /*- The reference wavefunction used in the computation -*/
    options.add_str("REFERENCE", "RHF");
    /*- The maximum number of iterations -*/
    options.add_int("MAXITER", 40);
    /*- An array containing the number of doubly-occupied orbitals per irrep (in Cotton order) -*/
    options.add("DOCC", new ArrayType());
    /*- An array containing the number of singly-occupied orbitals per irrep (in Cotton order) -*/
    options.add("SOCC", new ArrayType());
    /*- The amount of information to print to the output -*/
    options.add_int("PRINT", 3);
    /*- Whether to perturb the Hamiltonian or not -*/
    options.add_bool("PERTURB_H", false);
    options.add_double("LAMBDA", 0.0);
    /*- The storage scheme for the three index tensors in density fitting -*/
    options.add_str("RI_STORAGE", "DEFAULT", "DEFAULT DOUBLE_IN_CORE IN_CORE FLIP_B_DISK K_IN_CORE DISK");
    /*- The operator used to perturb the Hamiltonian, if requested -*/
    options.add_str("PERTURB_WITH", "DIPOLE_X", "DIPOLE_X DIPOLE_Y DIPOLE_Z");
    /*- Whether the computation is to be performed without disk-based storage of integrals -*/
    options.add_bool("DIRECT", false);
    /*- The maximum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("DIIS_VECTORS", 4);
    /*- Whether DIIS extrapolation is used to accelerate convergence -*/
    options.add_bool("DIIS", true);
    /*- Whether to use hard disk based storage in the computation -*/
    options.add_bool("OUT_OF_CORE", false);
    /*- The molecular charge -*/
    options.add_int("CHARGE", 0);
    /*- (2$\times M_s+1$), e.g. 1 for a singlet state, 2 for a doublet, 3 for a triplet, etc. -*/
    options.add_int("MULTP", 0);
    /*- Whether to print the molecular orbitals -*/
    options.add_bool("PRINT_MOS", false);
    /*- The amount of debugging information to print -*/
    options.add_int("DEBUG", 0);
    /*- -Log10 of the energy convergence criterion -*/
    options.add_int("E_CONVERGE", 8);
    /*- -Log10 of the density convergence criterion -*/
    options.add_int("D_CONVERGE", 8);
    /*- -Minimum absolute TEI value for seive -*/
    options.add_double("SCHWARZ_CUTOFF", 0.0);
    /*- -Minimum absolute D matrix value for DF-SCF exchange -*/
    //options.add_double("DENSITY_CUTOFF", 0.0);
    /*- -Minimum absolute S matrix value for DF-SCF exchange -*/
    options.add_double("OVERLAP_CUTOFF", 0.0);
    /*- -Minimum absolute three-index value for DF-SCF seive -*/
    options.add_double("THREE_INDEX_CUTOFF", 0.0);
    /*- -Maximum number of rows to read/write in each DF-SCf operation */
    options.add_int("ROWS_PER_READ", 1);
  }
  else if (name == "MP2") {
    options.add_str("WFN", "");
    options.add_str("REFERENECE", "RHF");
    options.add_str("JOBTYPE", "SP");
    options.add_str("DERTYPE", "NONE");
    options.add_int("PRINT", 0);
    options.add_int("CACHELEV", 2);
    options.add_str("CACHETYPE", "LRU", "LRU LOW");
    options.add_bool("SCS","false");
    options.add_bool("SCS_N", "false");
    options.add_double("SCALE_OS", 6.0/5.0);
    options.add_double("SCALE_SS", 1.0/3.0);
  }
  else if(name == "TRANSQT2") {
    options.add_str("WFN", "");
    options.add_str("REFERENCE","RHF");
    options.add_str("DERTYPE", "NONE");
    options.add_int("PRINT", 1);
    options.add_bool("PRINT_TEI", false);
    options.add_int("TOLERANCE", 14);
    options.add_int("CACHELEV", 2);
    options.add_str("AO_BASIS", "NONE");
    options.add_bool("DELETE_TEI", true);
  }
  else if(name == "CUSP"){
    options.add("FROZEN_DOCC", new ArrayType());
    options.add("FROZEN_UOCC", new ArrayType());
  }
  else if(name == "CCSORT") {
    options.add_str("WFN", "");
    options.add_str("REFERENCE", "RHF");
    options.add_str("DERTYPE", "NONE");
    options.add_str("PROPERTY", "POLARIZABILITY");
    options.add_bool("LOCAL", false);
    options.add_double("LOCAL_CUTOFF", 0.02);
    options.add_double("LOCAL_CPHF_CUTOFF", 0.10);
    options.add_double("LOCAL_CORE_CUTOFF",0.05);
    options.add_str("LOCAL_METHOD","WERNER");
    options.add_str("LOCAL_WEAKP","NONE");
    options.add_str("FREEZE_CORE","NONE");
    options.add_str("LOCAL_PAIRDEF","BP");
    options.add_bool("LOCAL_DOMAIN_POLAR", false);
    options.add_bool("LOCAL_DOMAIN_MAG", false);
    options.add_bool("LOCAL_DOMAIN_SEP", false);
    options.add_bool("LOCAL_FILTER_SINGLES", false);
    options.add_str("AO_BASIS","NONE");
    options.add_str("EOM_REFERENCE","RHF");
    options.add_int("PRINT", 0);
    options.add_bool("KEEP_TEIFILE", false);
    options.add_bool("KEEP_OEIFILE", false);
    options.add_int("TOLERANCE", 14);
    options.add_int("CACHELEV", 2);
    options.add_bool("LOCAL", false);
    options.add("OMEGA", new ArrayType());
  }
  else if(name == "CCTRIPLES") {
    options.add_str("WFN", "");
    options.add_int("NTHREADS",1);
    options.add_str("REFERENCE","RHF");
    options.add_str("DERTYPE","NONE");
  }
else if(name == "CCDENSITY") {
    options.add_str("WFN", "SCF");
    options.add_str("REFERENCE","RHF");
    options.add_str("DERTYPE","NONE");
    options.add_int("TOLERANCE",14);
    options.add_int("CACHELEVEL",2);
    options.add_bool("AO_BASIS",false);
    options.add_bool("AEL",false);
    options.add_str("GAUGE","LENGTH");
    options.add_bool("RELAX_OPDM",false);
    options.add_bool("CONNECT_XI",false);
    options.add("STATES_PER_IRREP", new ArrayType());
    options.add_bool("PROP_ALL",false);
    options.add_int("PROP_SYM", 0);
    options.add_int("PROP_ROOT", 0);
  }
  else if(name == "CCLAMBDA") {
    options.add_str("WFN","SCF");
    options.add_int("CONVERGENCE",7);
    options.add_bool("RESTART",false);
    options.add_int("PRINT",0);
    options.add_int("CACHELEVEL",2);
    options.add_bool("SEKINO",false);
    options.add_bool("DIIS",true);
    options.add_bool("AO_BASIS",false);
    options.add_str("ABCD","NEW");
    options.add_int("NUM_AMPS",10);
    options.add_str("DERTYPE","NONE");
    options.add_str("JOBTYPE","");
    options.add_bool("LOCAL",false);
    options.add_double("LOCAL_CUTOFF",0.02);
    options.add_str("LOCAL_METHOD","WERNER");
    options.add_bool("LOCAL_FILTER_SINGLES",true);
    options.add_double("LOCAL_CPHF_CUTOFF",0.10);
    options.add_str("FREEZE_CORE","FALSE");
    options.add_str("LOCAL_PAIRDEF","");
    options.add("STATES_PER_IRREP", new ArrayType());
    options.add_int("PROP_SYM",1);
    options.add_int("PROP_ROOT",1);
    options.add_int("MAXITER",50);
  }
  else if(name == "CLAG") {
    options.add_bool("WRITE_CAS_FILES",0);
    options.add_str("DERTYPE","NONE");
    options.add_str("WFN","NONE");
    options.add_int("ROOT",1);
  }
  else if(name == "STABLE") {
    options.add_int("PRINT",1);
    options.add_int("CACHELEV",2);
    options.add_str("REFERENCE",0);
    options.add_bool("FOLLOW",false);
    options.add_int("NUM_EVECS_PRINT",0);
    options.add_int("ROTATION_METHOD",0);
    options.add_double("SCALE",0.5);
  }
  else if(name == "OEPROP") {
    options.add_int("NUM_ROOTS",0);
    options.add_int("ROOT",0);
    options.add_int("GRID",0);
    options.add_str("MO_TO_PLOT","");
    options.add_int("GRID_ORIGIN",0);
    options.add_int("GRID_UNIT_X",0);
    options.add("GRID_XY0", new ArrayType());
    options.add("GRID_XY1", new ArrayType());
    options.add("GRID_XYZ0", new ArrayType());
    options.add("GRID_XYZ1", new ArrayType());
    options.add_int("NIX",0);
    options.add_int("NIY",0);
    options.add_int("NIZ",0);
    options.add_str("GRID_FORMAT","");
    options.add_double("GRID_ZMIN",0);
    options.add_double("GRID_ZMAX",0);
    options.add_int("EDGRAD_LOGSCALE",0);
    options.add_str("WFN","");
    options.add_int("TRANSITION_DENSITY",0);
    options.add_str("REFERENCE", "RHF");
    options.add_int("READ_OPDM", 1);
    options.add_double("OPDM_FILE", 0.0);
    options.add_str("OPDM_BASIS", "MO");
    options.add_str("OPDM_FORMAT", "SQUARE");
    options.add_int("WRTNOS", 0);
    options.add_int("ASYMM_OPDM", 0);
    options.add_int("SPIN_PROP", 0);
    options.add_double("PRINT", 1);
    options.add_int("PRINT_NOS", 0);
    options.add_int("CORREL_CORR", 0);
    options.add_double("ZVEC_FILE", 0);
    options.add_int("DELETE_ZVEC", 0);
    options.add_double("MPMAX", 0.0);
    options.add("MP_REF_XYZ", new ArrayType());
    options.add_double("MP_REF", 0.0);
    options.add("LM_REF_XYZ", new ArrayType());
    options.add_int("NUC_ESP", 0);
    options.add_double("FINE_STRUCTURE_ALPHA", 1/(_c_au));
    options.add_int("QED_DARWIN", 0);
    options.add_int("FREEZE_CORE", 0);
  }
  else if(name == "CCHBAR") {
    options.add_bool("TAMPLITUDE",false);
    options.add_int("CACHELEV",2);
    options.add_int("PRINT",0);
    options.add_str("WFN", "SCF");
    options.add_str("DERTYPE", "ENERGY");
    options.add_bool("WABEI_LOWDISK", false);
  }
  else if(name == "CCRESPONSE") {
    options.add_str("WFN", "SCF");
    options.add_int("PRINT",1);
    options.add_int("CACHELEV",2);
    options.add_str("REFERENCE","RHF");
    options.add_str("DERTYPE",0);
    options.add_str("GAUGE","LENGTH");
    options.add_int("MAXITER",50);
    options.add_int("CONVERGENCE",7);
    options.add_bool("DIIS",1);
    options.add_str("PROPERTY","POLARIZABILITY");
    options.add_str("ABCD","NEW");
    options.add_bool("RESTART",1);
    options.add_bool("LOCAL",0);
    options.add_double("LOCAL_CUTOFF",0.01);
    options.add_str("LOCAL_METHOD","WERNER");
    options.add_str("LOCAL_WEAKP","NONE");
    options.add_bool("LOCAL_FILER_SINGLES", false);
    options.add_double("LOCAL_CPHF_CUTOFF",0.10);
    options.add_str("FREEZE_CORE","FALSE");
    options.add_str("LOCAL_PAIRDEF","NONE");
    options.add_bool("ANALYZE",0);
    options.add_int("NUM_AMPS",5);
    options.add_bool("SEKINO",0);
    options.add_bool("LINEAR",0);
    options.add("OMEGA",new ArrayType());
  }
  else if(name == "MVO") {
   options.add_str("WFN","CCSD");
   options.add_int("FZC_FILE", PSIF_OEI);
   options.add_bool("PRINT_MOS",false);
   options.add_int("PRINT",1);
   options.add_bool("OEI_ERASE",false);
   options.add_bool("FZC",true);
   options.add_bool("DELETE_RESTR_DOCC",true);
   options.add_bool("MP2NOS",false);
   options.add_bool("UNOS",false);
   options.add_double("FZC_FOCK_COEFF",1.0);
   options.add_double("FOCK_COEFF",0.0);
   options.add_bool("IVO",false);
   options.add_bool("CANONICAL",false);
   options.add("FROZEN_DOCC", new ArrayType());
   options.add("FROZEN_UOCC", new ArrayType());
   options.add("RESTRICTED_DOCC", new ArrayType());
   options.add("RESTRICTED_UOCC", new ArrayType());
   options.add("DOCC", new ArrayType());
   options.add("SOCC", new ArrayType());
   options.add("DOCC_VIRT", new ArrayType());
  }
  else if(name == "RESPONSE"){
    options.add_int("PRINT", 1);
    options.add_str("REFERENCE", "RHF");
    options.add("OMEGA", new ArrayType());
    options.add_str("PROPERTY","POLARIZABILITY");
  }
  else if(name == "MCSCF") {
    /*- The molecular charge -*/
    options.add_int("CHARGE", 0);
    /*- (2$\times M_s+1$), e.g. 1 for a singlet state, 2 for a doublet, 3 for a triplet, etc. -*/
    options.add_int("MULTP", 1);
    options.add_int("CONVERGENCE",9);
    options.add_int("LEVELSHIFT",0);
    /*- The amount of debugging information to print -*/
    options.add_int("DEBUG", 0);
    /*- -Log10 of the energy convergence criterion -*/
    options.add_int("E_CONVERGE", 12);
    /*- -Log10 of the density convergence criterion -*/
    options.add_int("D_CONVERGE", 12);
    options.add_int("MAXITER",100);
    options.add_int("NDIIS",7);
    options.add_int("ROOT",1);
    options.add_int("START_FAVG",5);
    options.add_int("TURN_ON_ACTV",0);
    options.add_int("ROTATE_MO_ANGLE",0);
    options.add_int("ROTATE_MO_IRREP",1);  // IRREP is one-based
    options.add_int("ROTATE_MO_P",1);      // P and Q are one-based
    options.add_int("ROTATE_MO_Q",2);
    options.add_bool("CI_DIIS",false);
    options.add_bool("USE_DIIS",true);
    options.add_bool("READ_MOS",true);
    options.add_bool("USE_FAVG",false);
    options.add_bool("CANONICALIZE_ACTIVE_FAVG",false);
    options.add_bool("CANONICALIZE_INACTIVE_FAVG",false);
    options.add_bool("INTERNAL_ROTATIONS",true);
    options.add_bool("FORCE_TWOCON",false);

    options.add_str("REFERENCE","RHF","RHF ROHF UHF TWOCON MCSCF GENERAL");
    options.add_str("WFN_SYM","1","A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");
  }
  else if(name == "EXTREMA") {
    options.add_str("COORDINATES","foo");
  }
  else if(name == "CCENERGY") {
    options.add_bool("NEWTRIPS", 1);
    options.add_str("WFN", "NONE", "CCSD CCSD_T EOM_CCSD LEOM_CCSD BCCD BCCD_T CC2 CC3 EOM_CC2 EOM_CC3 CCSD_MVD");
    options.add_str("REFERENCE", "RHF");
    options.add_bool("ANALYZE", 0);
    options.add_str("DERTYPE", "NONE", "NONE FIRST RESPONSE");
    options.add_int("PRINT", 0);
    options.add_int("MAXITER", 50);
    options.add_int("CONVERGENCE", 7);
    options.add_bool("RESTART",1);
    options.add_bool("FORCE_RESTART", 0);
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    options.add_int("CACHELEV", 2);
    options.add_str("CACHETYPE", "LOW", "LOW LRU");
    options.add_int("NTHREADS",1);
    options.add_bool("DIIS", true);
    options.add_bool("T2_COUPLED", false);
    options.add_str("PROPERTY", "POLARIZABILITY", "POLARIZABILITY ROTATION MAGNETIZABILITY ROA ALL");
    options.add_str("ABCD", "NEW", "NEW OLD");
    options.add_bool("LOCAL", 0);
    options.add_double("LOCAL_CUTOFF", 0.02);
    options.add_double("LOCAL_MOS", 0);
    options.add_str("LOCAL_METHOD", "WERNER", "WERNER AOBASIS");
    options.add_str("LOCAL_WEAKP", "NONE", "NONE NEGLECT MP2");
    //options.add_int("LOCAL_FILTER_SINGLES", 1);
    options.add_double("LOCAL_CPHF_CUTOFF", 0.10);
    options.add_bool("FREEZE_CORE", 0);
    options.add_str("LOCAL_PAIRDEF", "BP", "BP RESPONSE");
    options.add_int("NUM_AMPS", 10);
    options.add_int("BRUECKNER_CONV", 5);
    options.add_bool("PRINT_MP2_AMPS", 0);
    options.add_bool("PRINT_PAIR_ENERGIES", 0);
    options.add_bool("SPINADAPT_ENERGIES", false);
    options.add_bool("T3_WS_INCORE", 0);
    options.add_bool("SCSN_MP2", 0);
    options.add_bool("SCS_MP2", 0);
    options.add_bool("SCS_CCSD", 0);
    options.add_double("MP2_SCALE_OS",1.20);
    options.add_double("MP2_SCALE_SS",1.0/3.0);
    options.add_double("CC_SCALE_OS", 1.27);
    options.add_double("CC_SCALE_SS",1.13);
  }
  else if(name == "CIS") {
    options.add_str("WFN", "CIS", "CCSD CCSD_T EOM_CCSD CIS");
    options.add_str("REFERENCE", "RHF", "RHF ROHF UHF");
    options.add_double("LOCAL_AMP_PRINT_CUTOFF", 0.60);
    options.add_int("PRINT", 0);
    options.add_int("MAXITER", 500);
    options.add_int("CONVERGENCE", 7);
    options.add("STATES_PER_IRREP", new ArrayType());
    options.add_str("DIAG_METHOD", "DAVIDSON", "DAVIDSON FULL");
    options.add_bool("LOCAL", false);
    options.add_double("LOCAL_CUTOFF", 0.02);
    options.add_str("LOCAL_METHOD", "WERNER", "AOBASIS WERNER");
    options.add_str("LOCAL_WEAKP", "MP2", "MP2 NEGLECT NONE");
    options.add_int("LOCAL_GHOST", -1);
    options.add("DOMAINS", new ArrayType());
    options.add_bool("DOMAIN_PRINT", 0);
  }
  else if(name == "LMP2") {
    /*- The wavefunction desired -*/
    options.add_str("RI_BASIS_MP2", "NONE");
    options.read_ipv1();
    if(options.get_str("RI_BASIS_MP2") != "NONE")
      options.add_bool("RI_LMP2", true);
    else
      options.add_bool("RI_LMP2", false);
    options.add_str("WFN", "LMP2");
    options.add_str("REFERENCE", "RHF", "RHF");
    options.add_int("PRINT", 0);
    options.add_int("MAXITER", 50);
    options.add_int("ENERGY_CONV", 7);
    options.add_int("RMS_CONV", 5);
    options.add_int("FSKIP", 2);
    options.add_bool("USE_DIIS", 1);
    options.add_bool("NEGLECT_DP", 1);
    options.add_double("DISTANT_PAIR", 8.0);
    options.add_int("DIISSTART", 3);
    options.add_int("NDIIS", 6);
    options.add_double("LOCAL_CUTOFF", 0.02);
    options.add_int("MEMORY", 2000);
    options.add_bool("SCS","false");
    options.add_bool("SCS_N", "false");
    options.add_double("SCALE_OS", 6.0/5.0);
    options.add_double("SCALE_SS", 1.0/3.0);
    options.add_int("SCREENING", 7);
    options.add_bool("SCREEN_INTS", false);
    options.add_int("SCHWARTZ_TOL", 12);
   }
  else if(name=="DF-MP2") {
    options.add_str("WFN", "RI-MP2");
    options.add_str("RI_BASIS_MP2","NONE");
    options.add_str("BASIS","NONE");
    options.add_bool("SCS","false");
    options.add_bool("SCS_N", "false");
    options.add_double("SCALE_OS", 6.0/5.0);
    options.add_double("SCALE_SS", 1.0/3.0);
  }
  else if(name == "PSIMRCC") {
    options.add_int("CORR_CHARGE",0);
    options.add_int("DEBUG",0);
    options.add_int("DAMPING_FACTOR",0);
    options.add_int("MAXDIIS",7);
    options.add_int("NUM_THREADS",1);
    options.add_int("NEL",0);
    options.add_int("ROOT",1);
    options.add_int("CONVERGENCE",9);
    options.add_int("MAXITER",100);
    options.add_int("DENOMINATOR_SHIFT",0);
    options.add_int("START_DIIS",2);
    options.add_int("TIKHONOW_OMEGA",0);  // Omega = TIKHONOW_OMEGA / 1000
    options.add_int("TIKHONOW_MAX",5);

    options.add_bool("DIIS_TRIPLES",false);
    options.add_bool("LOCK_SINGLET",false);
    options.add_bool("MP2_GUESS",true);
    options.add_bool("FAVG_CCSD_T",false);
    options.add_bool("HEFF4",true);
    options.add_bool("OFFDIAGONAL_CCSD_T",true);
    options.add_bool("DIAGONAL_CCSD_T",true);
    options.add_bool("DIAGONALIZE_HEFF",false);
    options.add_bool("ONLY_CLOSED_SHELL",false);
    options.add_bool("USE_DIIS",true);
    options.add_bool("USE_SPIN_SYMMETRY",true);
    options.add_bool("ZERO_INTERNAL_AMPS",true);
    options.add_bool("COUPLING_TERMS",true);
    options.add_bool("PRINT_HEFF",false);
    options.add_bool("PERT_CBS",false);
    options.add_bool("PERT_CBS_COUPLING",true);
    options.add_bool("RESTRICTED_TRIPLES",false);
    options.add_bool("TIKHONOW_TRIPLES",false);

    options.add_str("WFN","MRCCSD","MRCCSD");
    options.add_str("PT_ENERGY","SECOND_ORDER","SECOND_ORDER SCS_SECOND_ORDER PSEUDO_SECOND_ORDER SCS_PSEUDO_SECOND_ORDER");
    options.add_str("CORR_WFN","CCSD","PT2 CCSD MP2-CCSD CCSD_T");
    options.add_str("CORR_CCSD_T","STANDARD","STANDARD PITTNER");
    options.add_str("CORR_REFERENCE","GENERAL","RHF ROHF TCSCF MCSCF GENERAL");
    options.add_str("CORR_ANSATZ","MK","SR MK BW APBW");
    options.add_str("COUPLING","CUBIC","NONE LINEAR QUADRATIC CUBIC");
    options.add_str("WFN_SYM","1","A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");
    options.add_str("TRIPLES_ALGORITHM","RESTRICTED","SPIN_ADAPTED RESTRICTED UNRESTRICTED");
    options.add_str("MP2_CCSD_METHOD","II","I IA II");
  }
  options.read_ipv1();
  if(Communicator::world->me() == 0)
    options.print();
}

} //end ::psi

