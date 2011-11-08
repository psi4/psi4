/*! \file read_calculation_options
    \defgroup PSI4
*/

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
 * @param suppress_printing - boolean to specify whether to print to output file [false]
 */
int read_options(const std::string &name, Options & options, bool suppress_printing)
{
//  options.clear();

  /*- Units used in geometry specification -*/
  options.add_str("UNITS", "ANGSTROMS", "BOHR AU A.U. ANGSTROMS ANG ANGSTROM");
  /*- The molecular charge -*/
  options.add_int("CHARGE", 0);
  /*- Spin multiplicity, (2S+1), e.g. 1 for a singlet state, 2 for a doublet, 3 for a triplet, etc. -*/
  options.add_int("MULTP", 1);

  /*- An array containing the number of doubly-occupied orbitals per irrep 
  (in Cotton order) -*/
  options.add("DOCC", new ArrayType());

  /*- An array containing the number of singly-occupied orbitals per irrep 
  (in Cotton order) -*/
  options.add("SOCC", new ArrayType());

  /*- An array containing the number of frozen doubly-occupied orbitals per 
  irrep (these are not excited in a CI, nor can they be optimized in 
  MCSCF -*/
  options.add("FROZEN_DOCC", new ArrayType());
  /*- An array containing the number of frozen unoccupied orbitals per 
  irrep (these are not populated in a CI, nor can they be optimized in 
  MCSCF -*/
  options.add("FROZEN_UOCC", new ArrayType());
  /*- The scope of core orbitals to freeze in later correlated computations
      FROZEN_DOCC trumps this option -*/
  options.add_int("NUM_FROZEN_DOCC", 0);
  /*- The scope of virtual orbitals to freeze in later correlated computations
      FROZEN_UOCC trumps this option -*/
  options.add_int("NUM_FROZEN_UOCC", 0);
  /*- The scope of core orbitals to freeze in later correlated computations
      FROZEN_DOCC or NUM_FROZEN_DOCC trumps this option  -*/
  options.add_str("FREEZE_CORE","FALSE", \
    "FALSE TRUE SMALL LARGE");

  /*- Whether to use pure angular momentum basis functions -*/
  options.add_bool("PUREAM", true);
  /*- The amount of information to print to the output file -*/
  options.add_int("PRINT", 1);
  /*- The amount of information to print to the output file -*/
  options.add_int("DEBUG", 0);
  /*- Default number of geometry optimization steps -*/
  options.add_int("GEOM_MAXITER", 20);
  /*- Wavefunction type -*/
  options.add_str("WFN", "SCF");
  /*- Derivative level -*/
  options.add_str("DERTYPE", "NONE", "NONE FIRST SECOND RESPONSE");
  /*- Number of columns to print in calls to Matrix::print_mat -*/
  options.add_int("PRINT_MAT_NCOLUMN", 5);

  // CDS-TODO: We should go through and check that the user hasn't done
  // something silly like specify frozen_docc in DETCI but not in TRANSQT.
  // That would create problems.  (This was formerly checked in DETCI
  // itself, but I don't think DETCI will have the info available to check
  // this anymore).  This problem has affected users in the past.
  // Same goes for restricted_docc, restricted_uocc, ras1, ras2, ras3,
  // frozen_uocc.

  if (name == "DETCI" || options.read_globals()) {
    /*- Reference wavefunction -*/
    options.add_str("REFERENCE","RHF", "RHF ROHF");

    /*- Convergence is achieved when the RMS of the error in the CI vector is
    less than 10**(-n).  The default is 4 for energies and 7 for gradients. -*/
    options.add_int("CONVERGENCE", 4);

    /*- -Log10 of the energy convergence criterion -*/
    options.add_int("E_CONVERGE", 6);

    /*- Wavefunction type -*/
    options.add_str("WFN", "DETCI", "DETCI CI ZAPTN DETCAS CASSCF RASSCF");

    /*- Do a full CI (FCI)? -*/
    options.add_bool("FCI",false);

    /*- The CI excitation level -*/
    options.add_int("EX_LVL", 2);

    /*- The CC excitation level -*/
    options.add_int("CC_EX_LVL", 2);

    /*- The valence excitation level -*/
    options.add_int("VAL_EX_LVL", 0);

    /*- The CC valence excitation level -*/
    options.add_int("CC_VAL_EX_LVL", 0);

    /*- number of CI roots to find -*/
    options.add_int("NUM_ROOTS", 1);

    /*- The amount of information to print to the output file -*/
    options.add_int("PRINT", 1);

    /*- An array giving the number of orbitals per irrep for RAS1 !expert -*/
    options.add("RAS1", new ArrayType());

    /*- An array giving the number of orbitals per irrep for RAS2 !expert -*/
    options.add("RAS2", new ArrayType());

    /*- An array giving the number of orbitals per irrep for RAS3 !expert -*/
    options.add("RAS3", new ArrayType());

    /*- An array giving the number of orbitals per irrep for RAS4 !expert -*/
    options.add("RAS4", new ArrayType());

    /*- An array giving the number of restricted doubly-occupied orbitals per 
    irrep (not excited in CI wavefunctions, but orbitals can be optimized
    in MCSCF) -*/
    options.add("RESTRICTED_DOCC", new ArrayType());

    /*- An array giving the number of restricted unoccupied orbitals per 
    irrep (not occupied in CI wavefunctions, but orbitals can be optimized
    in MCSCF) -*/
    options.add("RESTRICTED_UOCC", new ArrayType());

    /*- An array giving the number of active orbitals (occupied plus 
    unoccupied) per irrep (shorthand to make MCSCF easier to specify than
    using RAS keywords) -*/
    options.add("ACTIVE", new ArrayType());

    /*- If TRUE, use the Ms=0 component of the state.  Defaults to TRUE
    if closed-shell and FALSE otherwise.  Related to the S parameter. -*/
    options.add_bool("MS0",false);

    /*- If TRUE then DETCI will stop after string information is formed
    and before integrals are read. -*/
    options.add_bool("ISTOP",false);

    /*- print a summary of the CI blocks? -*/
    options.add_bool("PRINT_CIBLKS",false);

    /*- Guess vector type.  Accepted values are UNIT for a unit vector
    guess (NUM_ROOTS and NUM_INIT_VECS must both be 1); H0_BLOCK to use
    eigenvectors from the H0 BLOCK submatrix (default); DFILE to use
    NUM_ROOTS previously converged vectors in the D file; IMPORT to
    import a guess previously exported from a CI computation
    (possibly using a different CI space) !expert -*/
    options.add_str("GUESS_VECTOR", "H0_BLOCK", "UNIT H0_BLOCK DFILE IMPORT");

    /*- The number of initial vectors to use in the CI iterative procedure.
    Defaults to the number of roots. !expert -*/
    options.add_int("NUM_INIT_VECS", 0);

    /*- Irrep for CI vectors;  -1 = find automatically.
    This option allows the user to look for CI vectors of a different irrep
    than the reference.  This probably only makes sense for Full CI,
    and it would probably not work with unit vector guesses.  Numbering
    starts from zero for the totally-symmetric irrep. !expert -*/
    options.add_int("REF_SYM", -1);

    /*- This parameter specifies the size of the H0 block of the Hamiltonian
    which is solved exactly.  The n determinants with the lowest SCF
    energy are selected, and a submatrix of the Hamiltonian is formed
    using these determinants.  This submatrix is used to accelerate
    convergence of the CI iterations in the BOLSEN and MITRUSHENKOV
    iteration schemes, and also to find a good starting guess for the
    SEM method if GUESS_VECTOR = H0_BLOCK.  Defaults to 400.
    Note that the program may change the given size for Ms=0 cases
    (Ms0 = TRUE) if it determines that the H0 block includes only
    one member of a pair of determinants related by time reversal symmetry.
    For very small block sizes, this could conceivably eliminate the entire
    H0 block; the program should print warnings if this occurs. !expert -*/
    options.add_int("H0_BLOCKSIZE", 400);

    /*- size of H0 block for initial guess !expert -*/
    options.add_int("H0_GUESS_SIZE", 400);

    /*- Use coupling block in preconditioner? !expert -*/
    options.add_bool("H0_BLOCK_COUPLING",false);

    /*- Parameters which specifies the size of the coupling block
     within the generalized davidson preconditioner. !expert -*/
    options.add_int("H0_BLOCK_COUPLING_SIZE",0);

    /*- number of important determinants to print out -*/
    options.add_int("NPRINT",20);

    /*- number of important CC amps per ex lvl to print -*/
    options.add_int("CC_NPRINT",10);

    /*- How to average H diag energies over spin coupling sets.
      HD_EXACT uses the exact diagonal energies which results in expansion
      vectors which break spin symmetry. HD_KAVE averages the diagonal
      energies over a spin-coupling set yielding spin pure expansion vectors.
      ORB_ENER employs the sum of orbital energy approximation giving
      spin pure expansion vectors but usually doubles the number of Davidson
      iterations. EVANGELISTI uses the sums and differences of orbital
      energies with the SCF reference energy to produce spin pure expansion
      vectors. LEININGER approximation which subtracts the one-electron
      contribution from the orbital energies, multiplies by 0.5, and adds
      the one-electron contribution back in, producing spin pure expansion
      vectors and developed by Matt Leininger and works as well as
      EVANGELISTI. !expert -*/
    options.add_str("HD_AVE", "EVANGELISTI",
      "EVANGELISTI HD_EXACT HD_KAVE ORB_ENER LEININGER Z_KAVE");

    /*- If TRUE the diagonal elements of the Hamiltonian matrix are
      computed on-the-fly, otherwise a diagonal element vector is written
      to a separate file on disk. !expert -*/
    options.add_bool("HD_OTF",true);

    /*- If TRUE, use the last vector space in the BVEC file to write
      scratch DVEC rather than using a separate DVEC file. (Only
      possible if NUM_ROOTS = 1). !expert -*/
    options.add_bool("NODFILE",false);

    /*- Freeze core orbitals? -*/
    // CDS-TODO: Need to make DETCI compatible with normal FREEZE_CORE
    options.add_bool("DETCI_FREEZE_CORE",true);

    /*- Store strings specifically for FCI? (Default to TRUE for FCI)
        !expert -*/
    options.add_bool("FCI_STRINGS",false);

    /*- This determines whether `mixed' RAS II/RAS III excitations are
      allowed into the CI space.  This is useful for placing additional
      constraints on a RAS CI. !expert -*/
    options.add_bool("MIXED",true);

    /*- This determines whether `mixed' excitations involving RAS IV are
      allowed into the CI space.  This is useful for placing additional
      constraints on a RAS CI. !expert -*/
    options.add_bool("MIXED4",true);

    /*- Restrict strings with e- in RAS IV: i.e. if an electron is in
      RAS IV, then the holes in RAS I must equal the particles in RAS III
      + RAS IV else the string is discarded !expert -*/
    options.add_bool("R4S",false);

    /*- Tells DETCI whether or not to do string replacements on the fly.  Can
      save a gigantic amount of memory (especially for truncated CI's) but
      is somewhat flaky and hasn't been tested for a while.  It may work
      only works for certain classes of RAS calculations.  The current
      code is very slow with this option turned on. !expert -*/
    options.add_bool("REPL_OTF",false);

    /*- If TRUE, calculate the value of $<S^2>$ for each root -*/
    options.add_bool("CALC_SSQ",false);

    /*- When this option is TRUE DETCI will compute the MPn series out to
    kth order where k is determined by maxnvect.  For open-shell systems
    (REF=ROHF, WFN = ZAPTN), DETCI will compute the ZAPTn series.
    GUESS_VECTOR must be set to UNIT, HD_OTF must be set to TRUE, and
    HD_AVE must be set to orb_ener; these should happen by default for
    MPN=TRUE. -*/
    options.add_bool("MPN",false);

    /*- If TRUE, save MP(2n-1) energy; else, save MPn energy -*/
    options.add_bool("SAVE_MPN2",false);

    /*- If TRUE, an orthonormal vector space is employed rather than
      storing the kth order wfn !expert -*/
    options.add_bool("MPN_SCHMIDT",false);

    /*- Use Wigner formulas in the Empn series? !expert -*/
    options.add_bool("WIGNER",false);

    /*- z in H = H0 + z * H1 !expert -*/
    options.add_double("PERTURBATION_PARAMETER",1.0);

    /*- maximum number of alpha electrons in RAS III -*/
    options.add_int("A_RAS3_MAX",-1);

    /*- maximum number of beta electrons in RAS III -*/
    options.add_int("B_RAS3_MAX",-1);

    /*- maximum number of alpha electrons in RAS III, for CC -*/
    options.add_int("CC_A_RAS3_MAX",-1);

    /*- maximum number of beta electrons in RAS III, for CC -*/
    options.add_int("CC_B_RAS3_MAX",-1);

    /*- maximum number of electrons in RAS III -*/
    options.add_int("RAS3_MAX",-1);

    /*- maximum number of electrons in RAS III, for CC -*/
    options.add_int("CC_RAS3_MAX",-1);

    /*- maximum number of electrons in RAS IV -*/
    options.add_int("RAS4_MAX",-1);

    /*- maximum number of electrons in RAS IV, for CC -*/
    options.add_int("CC_RAS4_MAX",-1);

    /*- maximum number of electrons in RAS III + IV -*/
    options.add_int("RAS34_MAX",-1);

    /*- maximum number of electrons in RAS III + IV, for CC -*/
    options.add_int("CC_RAS34_MAX",-1);

    /*- Specifies how to handle buffering of CI vectors.  A value of 0
    makes the program perform I/O one RAS subblock at a time; 1
    uses entire CI vectors at a time; and 2 uses one irrep block
    at a time.  Values of 0 or 2 cause some inefficiency in the I/O
    (requiring multiple reads of the C vector when constructing
    H in the iterative subspace if DIAG_METHOD = SEM), but require
    less core memory. -*/
    options.add_int("ICORE", 1);

    /*- This specifies which method is to be used in diagonalizing the
    Hamiltonian.  The valid options are: RSP, to form the entire H
    matrix and diagonalize using libciomr to obtain all eigenvalues
    (n.b. requires HUGE memory); OLSEN, to use Olsen's preconditioned
    inverse subspace method (1990); MITRUSHENKOV, to use a 2x2
    Olsen/Davidson method; and DAVIDSON (or SEM) to use Liu's
    Simultaneous Expansion Method, which is identical to the Davidson method
    if only one root is to be found.  There also exists a SEM debugging mode,
    SEMTEST.  The SEM method is the most robust, but it also
    requires 2(N*M)+1 CI vectors on disk, where N is the maximum number of
    iterations and M is the number of roots. -*/
    options.add_str("DIAG_METHOD", "SEM",
      "RSP OLSEN MITRUSHENKOV DAVIDSON SEM SEMTEST");

    /*- This specifies the type of preconditioner to use in the selected
    diagonalization method.  The valid options are: DAVIDSON which
    approximates the Hamiltonian matrix by the diagonal elements;
    H0BLOCK_INV which uses an exact Hamiltonian of H0_BLOCKSIZE and
    explicitly inverts it; GEN_DAVIDSON which does a spectral
    decomposition of H0BLOCK; ITER_INV using an iterative approach
    to obtain the correction vector of H0BLOCK.  The H0BLOCK_INV, GEN_DAVIDSON,
    and ITER_INV approaches are all formally equivalent but the ITER_INV is
    less computationally expensive.  Default is DAVIDSON. -*/
    options.add_str("PRECONDITIONER", "DAVIDSON",
      "LANCZOS DAVIDSON GEN_DAVIDSON H0BLOCK H0BLOCK_INV ITER_INV H0BLOCK_COUPLING EVANGELISTI");

    /*- DAVIDSON employs the standard DAVIDSON update or correction vector
    formula, while OLSEN uses the OLSEN correction vector.  Default
    is DAVIDSON. -*/
    options.add_str("UPDATE", "DAVIDSON", "DAVIDSON OLSEN");

    /*- Gives the maximum number of Davidson subspace vectors which can
    be held on disk for the CI coefficient and sigma vectors.  (There
    is one H(diag) vector and the number of D vectors is equal to the
    number of roots).  When the number of vectors on disk reaches
    the value of MAXNVECT, the Davidson subspace will be
    collapsed to COLLAPSE_SIZE vectors for each root.  This is very
    helpful for saving disk space.  Defaults to MAXITER * NUM_ROOTS
    + NUM_INIT_VECS. -*/
    options.add_int("MAXNVECT", 0);

    /*- Gives the number of vectors to retain when the Davidson subspace is
    collapsed (see MAXNVECT below).  If greater than one, the
    collapsed subspace retains the best estimate of the CI vector for
    the previous n iterations.   Defaults to 1. -*/
    options.add_int("COLLAPSE_SIZE", 1);

    /*- Use least-squares extrapolation in iterative solution of CI
    vector? -*/
    options.add_bool("LSE",false);

    /*- Number of iterations between least-squares extrapolations -*/
    options.add_int("LSE_COLLAPSE", 3);

    /*- Energy must be converged to $10^{-n}$ for least-squares
    extrapolation to be performed -*/
    options.add_int("LSE_TOLERANCE", 3);

    /*- This option allows the user to resume a DETCI iteration that
    terminated prematurely.  It assumes that the CI and sigma vectors are on
    disk; the number of vectors specified by RESTART_VECS is collapsed
    down to one vector per root. -*/
    options.add_bool("RESTART",false);

    /*- Use some routines to calculate sigma based on the papers of Bendazzoli
    et al.  Seems to be slower and not worthwhile; may disappear
    eventually.  Works only for full CI and I don't remember if I could see
    how their clever scheme might be extended to RAS in general. !expert -*/
    options.add_bool("BENDAZZOLI", false);

    /*- Do coupled-cluster computation? -*/
    options.add_bool("CC", false);

    /*- Compute one-particle density matrix if not otherwise required? -*/
    options.add_bool("OPDM", false);

    /*- Compute two-particle density matrix if not otherwise required? -*/
    options.add_bool("TPDM", false);

    /*- Maximum number of iterations to diagonalize the Hamiltonian. -*/
    options.add_int("MAXITER", 12);

    /*- Print the one-particle density matrix for each root? -*/
    options.add_bool("OPDM_PRINT", false);

    /*- Write the natural orbitals? -*/
    options.add_bool("WRTNOS", false);

    /*- Flag for whether or not to average the OPDM over several roots in
    order to obtain a state-average one-particle density matrix.  This
    density matrix can be diagonalized to obtain the CI natural orbitals. -*/
    options.add_bool("OPDM_AVE", false);

    /*- Sets the root number for which CI natural orbitals are written
    to PSIF_CHKPT.  The default value is 1 (lowest root). -*/
    options.add_int("ORBS_ROOT", 1);

    /*- Compute the kinetic energy contribution from the correlated part of
    the one-particle density matrix !expert -*/
    options.add_bool("OPDM_KE", false);

    /*- Print the two-particle density matrix? (Warning: large tensor) -*/
    options.add_bool("TPDM_PRINT", false);

    /*- The root to write out the two-particle density matrix for
    (the one-particle density matrices are written for all roots).
    Useful for a state-specific CASSCF or CI optimization on an
    excited state. -*/
    options.add_int("ROOT", 1);

    /*- Compute the transition density? -*/
    options.add_bool("TRANSITION_DENSITY", false);

    /*- Write the transition density? -*/
    options.add_bool("TDM_WRITE", false);

    /*- Print the transition density? -*/
    options.add_bool("TDM_PRINT", false);

    /*- Compute the dipole moment? -*/
    options.add_bool("DIPMOM", false);

    /*- Number of threads -*/
    options.add_int("NTHREADS", 1);

    /*- This specifies whether to store converged vector(s) at the end of
    the run.  The vector(s) is(are) stored in a transparent format such that
    other programs can use it easily. The format is specified in
    src/lib/libqt/slaterdset.h. The default is false. -*/
    options.add_bool("EXPORT_VECTOR", false);

    /*- Number of vectors to export -*/
    options.add_int("NUM_EXPORT", 1);

    /*- Eliminate determinants not valid for spin-complete spin-flip CI's
    [see J. S. Sears et al, J. Chem. Phys. 118, 9084-9094 (2003)] !expert -*/
    options.add_bool("SF_RESTRICT", false);

    /*- Print the sigma overlap matrix?  Not generally useful.  !expert -*/
    options.add_bool("SIGMA_OVERLAP", false);

    /*- The value of the spin quantum number S is given by this option.
    The default is determined by the value of MULTP.  This is used for two
    things: (1) determining the phase of the redundant half of the CI vector
    when the Ms=0 component is used (i.e., Ms0 = TRUE), and (2) making sure
    the guess vector has the desired value of $<S^2>$ (if CALC_SSQ is TRUE
    and ICORE=1). -*/
    options.add_double("S", 0.0);

    /*- An array of length EX_LVL specifying whether each excitation type
    (S,D,T, etc.) is allowed (1 is allowed, 0 is disallowed).  Used to
    specify non-standard CI spaces such as CIST.  !expert -*/
    options.add("EX_ALLOW", new ArrayType());

    /*- The FILTER_GUESS options are used to filter out some trial
    vectors which may not have the appropriate phase convention
    between two determinants.  This is useful to remove, e.g.,
    delta states when a sigma state is desired.  The user
    inputs two determinants (by giving the absolute alpha string
    number and beta string number for each), and also the
    desired phase between these two determinants for guesses
    which are to be kept.  FILTER_GUESS = TRUE turns on the filtering
    routine (and requires additional keywords). !expert -*/
    options.add_bool("FILTER_GUESS", false);

    /*- The required phase (1 or -1) between the two determinants specified
    by FILTER_GUESS_DET1 and FILTER_GUESS_DET2 !expert -*/
    options.add_int("FILTER_GUESS_SIGN", 1);

    /*- Array specifying the absolute alpha string number and beta string
    number for the first determinant in the filter procedure. !expert -*/
    options.add("FILTER_GUESS_DET1", new ArrayType());

    /*- Array specifying the absolute alpha string number and beta string
    number for the second determinant in the filter procedure. !expert -*/
    options.add("FILTER_GUESS_DET2", new ArrayType());

    /*- If present, the code will try to filter out a particular determinant
    by setting its CI coefficient to zero.  FILTER_ZERO_DET = (alphastr
    betastr) specifies the absolute alpha and beta string numbers of the
    target determinant. This could be useful for trying to exclude states
    that have a nonzero CI coefficient for the given determinant.  However,
    this option was experimental and may not be effective.  !expert -*/
    options.add("FILTER_ZERO_DET", new ArrayType());

    /*- Array giving the root numbers of the states to average in a
    state-averaged procedure such as SA-CASSCF. Root numbering starts
    from 1. -*/
    options.add("AVERAGE_STATES", new ArrayType());

    /*- Array giving the weights for each state in a state-averaged
    procedure -*/
    // CDS:TODO - Does this work for doubles??
    options.add("AVERAGE_WEIGHTS", new ArrayType());

    /*- In following a particular root (see ROOT keyword), sometimes the
    root number changes.  To follow a root of a particular character,
    one can specify a list of determinants and their coefficients,
    and the code will follow the root with the closest overlap.  The
    user specifies arrays containing the absolute alpha string indices
    (A_i below), absolute beta indices (B_i below), and CI coefficients
    (C_i below) to form the desired vector.
    FOLLOW_VECTOR_ALPHAS specifies the alpha string indices. The
    format is FOLLOW_VECTOR = [ [[A_1, B_1], C_1], [[A_2, B_2], C_2], ...].
    !expert -*/
    options.add("FOLLOW_VECTOR", new ArrayType());

    /*- Export a CC vector to disk? -*/
    options.add_bool("CC_EXPORT", false);

    /*- Import a CC vector from disk? -*/
    options.add_bool("CC_IMPORT", false);

    /*- Fix amplitudes involving RAS I or RAS IV?  Useful in mixed
    MP2-CC methods. !expert -*/
    options.add_bool("CC_FIX_EXTERNAL", false);

    /*- Number of external indices before amplitude gets fixed by
    CC_FIX_EXTERNAL.  Experimental. !expert -*/
    options.add_int("CC_FIX_EXTERNAL_MIN", 1);

    /*- Use variational energy expression in CC computation?
    Experimental.  !expert -*/
    options.add_bool("CC_VARIATIONAL", false);

    /*- ignore block if num holes in RAS I and II is > cc_ex_lvl and if
    any indices correspond to RAS I or IV (i.e., include only all-active
    higher excitations) !expert -*/
    options.add_bool("CC_MIXED", true);

    /*- Update T amplitudes with orbital eigenvalues? (Usually would
    do this).  Not doing this is experimental.  !expert -*/
    options.add_bool("CC_UPDATE_EPS", true);

    /*- Do DIIS? -*/
    options.add_bool("DIIS", true);

    // CDS-TODO: Check difference between DIIS_START AND DIIS_MIN_VECS
    /*- how many diis vectors built up before start -*/
    options.add_int("DIIS_START", 1);

    /*- how often to do a DIIS exterpolation.  1 means do DIIS every
    iteration, 2 is every other iteration, etc. -*/
    options.add_int("DIIS_FREQ", 1);

    /*- how many vectors required before do diis? -*/
    options.add_int("DIIS_MIN_VECS", 2);

    /*- how many vectors maximum to hold? -*/
    options.add_int("DIIS_MAX_VECS", 5);

    /*- CC_MACRO = [ [ex_lvl, max_holes_I, max_parts_IV, max_I+IV],
                     [ex_lvl, max_holes_I, max_parts_IV, max_I+IV], ... ]
    Optional additional restrictions on allowed exictations in
    coupled-cluster computations, based on macroconfiguration selection.
    For each sub-array, [ex_lvl, max_holes_I, max_parts_IV, max_I+IV],
    eliminate cluster amplitudes in which: [the excitation level
    (holes in I + II) is equal to ex_lvl] AND [there are more than
    max_holes_I holes in RAS I, there are more than max_parts_IV
    particles in RAS IV, OR there are more than max_I+IV quasiparticles
    in RAS I + RAS IV].  !expert -*/
    options.add("CC_MACRO", new ArrayType());

  }

  if (name == "SAPT"|| options.read_globals()) {
    /*- The level of theory for SAPT -*/
    options.add_str("SAPT_LEVEL","SAPT0","SAPT0 SAPT2 SAPT2+ SAPT2+3 MP2C");
    /*- The ubiquitous debug flag -*/
    options.add_int("DEBUG",0);
    /*- The amount of information to print to the output file -*/
    options.add_int("PRINT",1);
    /*- How many digits after the decimal to converge the energy to -*/
    options.add_int("E_CONVERGE",10);
    /*- How many digits after the decimal to converge the density to -*/
    options.add_int("D_CONVERGE",8);
    /*- Don't solve the CPHF equations -*/
    options.add_bool("NO_RESPONSE",false);
    /*- Use asynchronous I/O in the CPHF solver -*/
    options.add_bool("AIO_CPHF",false);
    /*- Use asynchronous I/O in the DF integral formation -*/
    options.add_bool("AIO_DFINTS",false);
    /*- Max CPHF iterations -*/
    options.add_int("MAXITER",50);
    /*- The number of DIIS vectors used to extrapolate -*/
    options.add_int("DIISVECS",5);
    /*- Compute Third-order Corrections -*/
    options.add_bool("DO_THIRD_ORDER",false);
    /*- Compute Natural Orbitals -*/
    options.add_bool("NAT_ORBS",false);
    /*- Use Natural Orbitals for T2's -*/
    options.add_bool("NAT_ORBS_T2",false);
    /*- Natural Orbital Occupation Cutoff -*/
    options.add_double("OCC_CUTOFF",1.0E-6);
    /*- Schwarz cutoff -*/
    options.add_double("SCHWARZ_CUTOFF",1.0E-12);
    /*- Memory safety -*/
    options.add_double("SAPT_MEM_SAFETY",0.9);
    /*- SAPT DF Basis -*/
    options.add_str("RI_BASIS_SAPT", "");
    /*- SAPT DF Basis for Elst10 and Exch10 -*/
    options.add_str("RI_BASIS_ELST", "");
    /*- Maximum denominator error allowed (Max error norm in Delta tensor) -*/
    options.add_double("DENOMINATOR_DELTA", 1.0E-6);
    /*- Denominator algorithm for PT methods -*/
    options.add_str("DENOMINATOR_ALGORITHM", "LAPLACE", "LAPLACE CHOLESKY");
    /*- Number of Omega points for Casimir-Polder -*/
    options.add_int("OMEGA_POINTS",8);
    /*- Omega (atomic wavenumbers) to center Casimir-Polder on -*/
    options.add_double("OMEGA_CENTER", 0.4);
  }
  if(name == "DCFT"|| options.read_globals()) {
//      ip_cwk_add(":DCFT");
      /*- How to cache quantities within the DPD library -*/
      options.add_int("CACHELEV", 2);
      /*- The amount of memory available (in Mb) -*/
      options.add_int("MEMORY", 2000);
      /*- The shift applied to the denominator -*/
      options.add_double("REGULARIZER", 0.0);
      /*- The maximum number of lambda iterations per macro-iteration -*/
      options.add_int("LAMBDA_MAXITER", 50);
      /*- The maximum number of SCF iterations per cycle -*/
      options.add_int("SCF_MAXITER", 50);
      /*- The maximum number iterations allowed -*/
      options.add_int("MAXITER", 40);
      /*- Whether to compute the full two particle density matrix at the end of the computation, for properties -*/
      options.add_bool("COMPUTE_TPDM", 0);
      /*- The number of decimal digits required in the SCF density -*/
      options.add_int("SCF_CONV", 8);
      /*- The number of decimal digits required in the determination of lambda -*/
      options.add_int("CONVERGENCE", 10);
      /*- Whether to relax the orbitals or not -*/
      options.add_bool("RELAX_ORBITALS", true);
      /*- The damping factor used in the initial SCF procedure (0 - 1000) 0 means full standard SCF update
          is performed, 1000 will completely damp the iteration to the extent that no update is performed -*/
      options.add_int("DAMPING_FACTOR", 0);
      /*- Should the tau terms be included? -*/
      options.add_bool("IGNORE_TAU", false);
      /*- Whether to compute the DCFT energy with the $\tau^{2}$ correction to $\tau$ or not-*/
      options.add_bool("TAU_SQUARED", false);
      /*- An integral is considered to be zero if it's magnitude is less than $10^{-int_thresh}$ -*/
      options.add_int("INT_THRESH", 14);
      /*- DIIS starts when the  RMS lambda and SCF errors are less than $10^{diis_start}$ -*/
      options.add_int("DIIS_START", 3);
      /*- The maximum number of vectors used in DIIS extrapolation -*/
      options.add_int("MAX_DIIS", 6);
      /*- The number of DIIS vectors needed for extrapolation to start -*/
      options.add_int("DIIS_NUM_VECS", 3);
      /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
      options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
      /*- The algorithm to use for lambda and orbital updates -*/
      options.add_str("ALGORITHM", "SIMULTANEOUS", "TWOSTEP SIMULTANEOUS");
      /*- Whether to force the occupation to be that of the SCF starting point -*/
      options.add_bool("LOCK_OCCUPATION", true);
      /*- The molecular charge -*/
      options.add_int("CHARGE", 0);
      /*- (2$\times M_s+1$), e.g. 1 for a singlet state, 2 for a doublet, 3 for a triplet, etc. -*/
      options.add_int("MULTP", 0);
  }
  if (name == "MINTS"|| options.read_globals()) {
      /*- primary basis set -*/
      options.add_str("BASIS","");
  }
  if (name == "DFT" || options.read_globals()) {
    /*- The DFT grid specification, such as SG1 -*/
    options.add_str("DFT_GRID_NAME","","SG1");
    /*- Maximum order of spherical grids -*/
    options.add_int("DFT_ORDER_SPHERICAL", 15);
    /*- Number of radial points -*/
    options.add_int("DFT_N_RADIAL", 99);
    /*- Spherical Scheme -*/
    options.add_str("DFT_SPHERICAL_SCHEME", "LEBEDEV", "LEBEDEV");
    /*- Radial Scheme -*/
    options.add_str("DFT_RADIAL_SCHEME", "TREUTLER", "TREUTLER BECKE MULTIEXP EM MURA");
    /*- Nuclear Scheme -*/
    options.add_str("DFT_NUCLEAR_SCHEME", "TREUTLER", "TREUTLER BECKE NAIVE STRATMANN");
    /*- Pruning Scheme -*/
    options.add_str("DFT_PRUNING_SCHEME", "FLAT", "FLAT P_GAUSSIAN D_GAUSSIAN P_SLATER D_SLATER LOG_GAUSSIAN LOG_SLATER");
    /*- Factor for effective BS radius in radial grid -*/
    options.add_double("DFT_BS_RADIUS_ALPHA",1.0);
    /*- Spread alpha for logarithmic pruning -*/
    options.add_double("DFT_PRUNING_ALPHA",1.0);
    /*- The number of grid points per evaluation block -*/
    options.add_int("DFT_MAX_POINTS",5000);
    /*- The number of grid points per evaluation block -*/
    options.add_int("DFT_MIN_POINTS",0);
    /*- The DFT basis cutoff -*/
    options.add_double("DFT_BASIS_CUTOFF", 0.0);
    /*- The DFT combined functional name (for now) -*/
    options.add_str("DFT_FUNCTIONAL", "");
    /*- The DFT Range-separation parameter (only used if changed by the user) -*/
    options.add_double("DFT_OMEGA", 0.0);
  }
  if (name == "SCF"|| options.read_globals()) {

    /*- Are going to do SAPT? If so, what part?  -*/
    options.add_str("SAPT","FALSE","FALSE 2-DIMER 2-MONOMER_A 2-MONOMER_B 3-TRIMER 3-DIMER_AB 3-DIMER_BC 3-DIMER_AC 3-MONOMER_A 3-MONOMER_B 3-MONOMER_C");

    /*- The name of the auxiliary basis to be used in RI computations -*/
    options.add_str("RI_BASIS_SCF", "");

    /*- Atomic Charge cutoff (for primary domain) -*/
    options.add_double("CHARGE_CUTOFF",0.05);
    /*- Extended domain radius, Angstrom -*/
    options.add_double("R_EXT",3.0);
    /*- Iterations per full Pipek-Mizey Localization -*/
    options.add_int("STEPS_PER_LOCALIZE",1);
    /*- What algorithm to use for the SCF computation -*/
    options.add_str("SCF_TYPE","PK","PK OUT_OF_CORE DIRECT DF PSEUDOSPECTRAL POISSON L_DF CD 1C_CD");
    /*- Whether to run in parallel or not -*/
    options.add_bool("PARALLEL", false);

    /*- primary basis set -*/
    options.add_str("BASIS","");

    /*- The type of guess orbitals -*/
    options.add_str("GUESS", "CORE", "CORE GWH SAD READ");
    /*- The reference wavefunction used in the computation -*/
    options.add_str("REFERENCE", "RHF", "RHF ROHF UHF CUHF RKS UKS");
    /*- The maximum number of iterations -*/
    options.add_int("MAXITER", 100);

    /*- Whether to perturb the Hamiltonian or not -*/
    options.add_bool("PERTURB_H", false);
    /*- How big is the perturbation? -*/
    options.add_double("LAMBDA", 0.0);
    /*- The operator used to perturb the Hamiltonian, if requested -*/
    options.add_str("PERTURB_WITH", "DIPOLE_X", "DIPOLE_X DIPOLE_Y DIPOLE_Z");
    /*- Look for an ExternalPotential object in Python or not? -*/
    options.add_bool("EXTERN", false);

    /*- The storage scheme for the three index tensors in density fitting -*/
    options.add_str("RI_SCF_STORAGE", "DEFAULT", "DEFAULT CORE DISK");
    /*- Should we make sure to save restart information for RI SCF? -*/
    options.add_bool("RI_SCF_SAVE",false);
    /*- Should we try to restart -*/
    options.add_bool("RI_SCF_RESTART",false);
    /*- Should we find the raw fitting metric condition number? -*/
    options.add_bool("FIND_RAW_J_COND",false);
    /*- Max J basis condition number to be allowed -*/
    options.add_double("RI_MAX_COND",1E8);
    /*- SCF Fitting Type -*/
    options.add_str("RI_FITTING_TYPE", "FINISHED", "FINISHED RAW CHOLESKY");
    /*- Max Number of threads for integrals (may be turned down if memory is an issue). 0 is blank -*/
    options.add_int("RI_INTS_NUM_THREADS",0);
    /*- IO caching for CP corrections, etc -*/
    options.add_str("RI_INTS_IO", "NONE", "NONE SAVE LOAD");

    /*- SO orthogonalization: symmetric or canonical? -*/
    options.add_str("S_ORTHOGONALIZATION","SYMMETRIC","SYMMETRIC CANONICAL");
    /*- Minimum S matrix eigenvalue to be used before compensating for linear dependencies -*/
    options.add_double("S_MIN_EIGENVALUE",1E-7);

    /*- The iteration to start MOM on (or 0 for no MOM) -*/
    options.add_int("MOM_START", 0);
    /*- The absolute indices of orbitals to excite from in MOM (+/- for alpha/beta) -*/
    options.add("MOM_OCC", new ArrayType());
    /*- The absolute indices of orbitals to excite to in MOM (+/- for alpha/beta) -*/
    options.add("MOM_VIR", new ArrayType());

    /*- The minimum iteration to start storing DIIS vectors -*/
    options.add_int("START_DIIS_ITER", 1);
    /*- The minimum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("MIN_DIIS_VECTORS", 2);
    /*- The maximum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("MAX_DIIS_VECTORS", 10);
    /*- Whether DIIS extrapolation is used to accelerate convergence -*/
    options.add_bool("DIIS", true);
    /*- The molecular charge -*/
    options.add_int("CHARGE", 0);
    /*- (2$\times M_s+1$), e.g. 1 for a singlet state, 2 for a doublet, 3 for a triplet, etc. -*/
    options.add_int("MULTP", 0);
    /*- Whether to print the molecular orbitals -*/
    options.add_bool("PRINT_MOS", false);
    /*- The amount of debugging information to print -*/
    options.add_int("DEBUG", false);
    /*- -Log10 of the energy convergence criterion -*/
    options.add_int("E_CONVERGE", 8);
    /*- -Log10 of the density convergence criterion -*/
    options.add_int("D_CONVERGE", 8);
    /*- Minimum absolute TEI value for seive -*/
    options.add_double("SCHWARZ_CUTOFF", 0.0);
    /*- Minimum absolute S matrix value for DF-SCF exchange -*/
    options.add_double("OVERLAP_CUTOFF", 0.0);
    /*- Minimum absolute three-index value for DF-SCF seive -*/
    options.add_double("THREE_INDEX_CUTOFF", 0.0);
    /*- Maximum number of rows to read/write in each DF-SCF operation -*/
    options.add_int("ROWS_PER_READ", 0);

    /*- The amount of SAD information to print to the output -*/
    options.add_int("SAD_PRINT", 0);
    /*- SAD Occupation Matrix Method -*/
    options.add_str("SAD_C", "CHOLESKY", "CHOLESKY ID");
    /*- SAD Guess Convergence in E -*/
    options.add_double("SAD_E_CONVERGE", 1E-5);
    /*- SAD Guess Convergence in D -*/
    options.add_double("SAD_D_CONVERGE", 1E-5);
    /*- SAD Guess Maxiter -*/
    options.add_int("SAD_MAXITER", 50);
    /*- SAD Guess F-mix Iteration Start -*/
    options.add_int("SAD_F_MIX_START", 50);
    /*- SAD Guess Schwarz Sieve (for rough molecular F) -*/
    options.add_double("SAD_SCHWARZ_CUTOFF", 1E-7);
    /*- SAD Guess Cholesky Cutoff (for eliminating redundancies) -*/
    options.add_double("SAD_CHOL_CUTOFF", 1E-7);
  }
  if (name == "MP2"|| options.read_globals()) {
    /*- The wavefunction type -*/
    options.add_str("WFN", "MP2", "MP2");
    /*- The reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF", "RHF UHF ROHF");
    /*- The type of job being performed -*/
    options.add_str("JOBTYPE", "SP");
    /*- The order of energy derivatives required -*/
    options.add_str("DERTYPE", "NONE");
    /*- Whether to compute the one particle density matrix, for properties -*/
    options.add_bool("COMPUTE_OPDM", false);
    /*- Whether to add relaxation terms to the one particle density matrix, for properties -*/
    options.add_bool("RELAX_OPDM", false);
    /*- The amount of cacheing of data to perform -*/
    options.add_int("CACHELEV", 2);
    /*- The criterion used to retain/release cached data -*/
    options.add_str("CACHETYPE", "LRU", "LRU LOW");
    /*- Whether to perform a spin component scaled MP2 computation -*/
    options.add_bool("SCS","false");
    /*- Whether to perform a spin component scaled (N) MP2 computation -*/
    options.add_bool("SCS_N", "false");
    /*- The scale factor used for opposite-spin pairs in SCS computations -*/
    options.add_double("SCALE_OS", 6.0/5.0);
    /*- The scale factor used for same-spin pairs in SCS computations-*/
    options.add_double("SCALE_SS", 1.0/3.0);
  }
  if(name == "TRANSQT2"|| options.read_globals()) {
    /*- -*/
    options.add_str("WFN", "");
    /*- -*/
    options.add_str("REFERENCE","RHF");
    /*- -*/
    options.add_bool("PRINT_TEI", false);
    /*- -*/
    options.add_int("TOLERANCE", 14);
    /*- -*/
    options.add_int("CACHELEV", 2);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- -*/
    options.add_bool("DELETE_TEI", true);
  }
  if(name == "TRANSQT"|| options.read_globals()) {
    /*- -*/
    options.add_int("PRINT_LVL", 1);
    /*- -*/
    options.add_str("REFERENCE","RHF");
    /*- -*/
    options.add_str("MODE", "TO_MO", "TO_MO TO_AO");
    /*- -*/
    options.add_str("WFN", "CCSD");
    /*- -*/
    options.add_bool("PSIMRCC", false);
    /*- -*/
    options.add_str("MP2R12A", "MP2R12AERI", "MP2R12AERI MP2R12AR12 MP2R12AR12T1");
    /*- -*/
    options.add_int("TOLERANCE", 14);
    /*- -*/
    options.add_int("OEI_FILE", PSIF_OEI);
    /*- -*/
    options.add_int("OEI_A_FILE", PSIF_OEI);
    /*- -*/
    options.add_int("OEI_B_FILE", PSIF_OEI);
    /*- -*/
    options.add_int("FZC_FILE", PSIF_OEI);
    /*- -*/
    options.add_int("FZC_A_FILE", PSIF_OEI);
    /*- -*/
    options.add_int("FZC_B_FILE", PSIF_OEI);
    /*- -*/
    options.add_int("SORTED_TEI_FILE", PSIF_MO_TEI);
    /*- -*/
    options.add_int("TPDM_FILE", PSIF_MO_TPDM);
    /*- -*/
    options.add_int("SO_S_FILE", PSIF_OEI);
    /*- -*/
    options.add_int("SO_T_FILE", PSIF_OEI);
    /*- -*/
    options.add_int("SO_V_FILE", PSIF_OEI);
    /*- -*/
    options.add_int("SO_TEI_FILE", PSIF_SO_TEI); // ?
    /*- -*/
    options.add_int("FIRST_TMP_FILE", 150);
    /*- -*/
    options.add_int("OPDM_IN_FILE", PSIF_MO_OPDM);
    /*- -*/
    options.add_int("OPDM_OUT_FILE", PSIF_AO_OPDM);
    /*- -*/
    options.add_int("LAG_IN_FILE", PSIF_MO_LAG);
    /*- -*/
    options.add_int("PRESORT_FILE", PSIF_SO_PRESORT);
    /*- -*/
    options.add_bool("KEEP_PRESORT", false);
    /*- -*/
    options.add_int("J_FILE", 91);
    /*- -*/
    options.add_bool("KEEP_J", false); // keep half-transformed integrals
    /*- -*/
    options.add_int("M_FILE", 0); // output integrals file; depends on direction
    /*- -*/
    options.add_int("AA_M_FILE", PSIF_MO_AA_TEI);
    /*- -*/
    options.add_int("BB_M_FILE", PSIF_MO_BB_TEI);
    /*- -*/
    options.add_int("AB_M_FILE", PSIF_MO_AB_TEI);
    /*- -*/
    options.add_int("MAX_BUCKETS", 499);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- -*/
    options.add_bool("DELETE_AO", true);
    /*- -*/
    options.add_bool("DELETE_TPDM", true);

    /*- -*/
    options.add_bool("PRINT_TE_INTEGRALS", false);
    /*- -*/
    options.add_bool("PRINT_OE_INTEGRALS", false);
    /*- -*/
    options.add_bool("PRINT_SORTED_OE_INTS", false);
    /*- -*/
    options.add_bool("PRINT_SORTED_TE_INTS", false);
    /*- -*/
    options.add_bool("PRINT_MOS", false);

    /*- -*/
    options.add_bool("LAGRAN_DOUBLE", false);
    /*- -*/
    options.add_bool("LAGRAN_HALVE", false);
    /*- -*/
    options.add_bool("DO_ALL_TEI", false);
    /*- -*/
    options.add_bool("TPDM_ADD_REF", false);
    /*- -*/
    options.add_bool("DELETE_RESTR_DOCC", true);
    /*- -*/
    options.add_bool("PRINT_REORDER", false);
    /*- -*/
    options.add_bool("PITZER", false);
    /*- -*/
    options.add_bool("REORDER", false);
    /*- -*/
    options.add_bool("CHECK_C_ORTHONORM", false);
    /*- -*/
    options.add_bool("QRHF", false);
    /*- -*/
    options.add_bool("IVO", false);
    /*- -*/
    options.add("MOORDER", new ArrayType());

    /*- An array giving the number of orbitals per irrep for RAS1 !expert -*/
    options.add("RAS1", new ArrayType());

    /*- An array giving the number of orbitals per irrep for RAS2 !expert -*/
    options.add("RAS2", new ArrayType());

    /*- An array giving the number of orbitals per irrep for RAS3 !expert -*/
    options.add("RAS3", new ArrayType());

    /*- An array giving the number of orbitals per irrep for RAS4 !expert -*/
    options.add("RAS4", new ArrayType());

    /*- An array giving the number of restricted doubly-occupied orbitals per 
    irrep (not excited in CI wavefunctions, but orbitals can be optimized
    in MCSCF) -*/
    options.add("RESTRICTED_DOCC", new ArrayType());

    /*- An array giving the number of restricted unoccupied orbitals per 
    irrep (not occupied in CI wavefunctions, but orbitals can be optimized
    in MCSCF) -*/
    options.add("RESTRICTED_UOCC", new ArrayType());

  }
  if(name == "CUSP"|| options.read_globals()){
  }
  if(name == "CCSORT"|| options.read_globals()) {
    /*- -*/
    options.add_str("WFN", "");
    /*- -*/
    options.add_str("REFERENCE", "RHF");
    /*- -*/
    options.add_str("PROPERTY", "POLARIZABILITY");
    /*- -*/
    options.add_bool("LOCAL", false);
    /*- -*/
    options.add_double("LOCAL_CUTOFF", 0.02);
    /*- -*/
    options.add_double("LOCAL_CPHF_CUTOFF", 0.10);
    /*- -*/
    options.add_double("LOCAL_CORE_CUTOFF",0.05);
    /*- -*/
    options.add_str("LOCAL_METHOD","WERNER");
    /*- -*/
    options.add_str("LOCAL_WEAKP","NONE");
    /*- -*/
    options.add_str("LOCAL_PAIRDEF","BP");
    /*- -*/
    options.add_bool("LOCAL_DOMAIN_POLAR", false);
    /*- -*/
    options.add_bool("LOCAL_DOMAIN_MAG", false);
    /*- -*/
    options.add_bool("LOCAL_DOMAIN_SEP", false);
    /*- -*/
    options.add_bool("LOCAL_FILTER_SINGLES", false);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- -*/
    options.add_str("EOM_REFERENCE","RHF");
    /*- -*/
    options.add_bool("KEEP_TEIFILE", false);
    /*- -*/
    options.add_bool("KEEP_OEIFILE", false);
    /*- -*/
    options.add_int("TOLERANCE", 14);
    /*- -*/
    options.add_int("CACHELEV", 2);
    /*- -*/
    options.add_bool("LOCAL", false);
    /*- -*/
    options.add("OMEGA", new ArrayType());
    /*- -*/
    options.add_str("OMEGA_UNITS", "AU", "AU HZ EV NM");
  }
  if(name == "CCTRIPLES"|| options.read_globals()) {
    /*- The wavefunction type -*/
    options.add_str("WFN", "");
    /*- The number of threads to use on multi-core machines -*/
    options.add_int("NTHREADS",1);
    /*- The reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
  }
  if(name == "CCDENSITY"|| options.read_globals()) {
    /*- The wavefunction type -*/
    options.add_str("WFN", "SCF");
    /*- The reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
    /*- Integrals are neglected if their magnitude is below $10^{-tolerance}$ -*/
    options.add_int("TOLERANCE",14);
    /*- The amount of cacheing of data to perform -*/
    options.add_int("CACHELEV",2);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- The approximate excitation level, {\it c.f.} Stanton and Bartlett, JCP, 98, 1993, 7034. !expert -*/
    options.add_bool("AEL",false);
    /*- The type of gauge to use for properties -*/
    options.add_str("GAUGE","LENGTH");
    /*- Whether to relax the one-particle density matrix -*/
    options.add_bool("RELAX_OPDM",false);
    /*- Whether $\bar{H}$ and R must be connected !expert -*/
    options.add_bool("CONNECT_XI",false);
    /*- The number of electronic states to computed, per irreducible representation -*/
    options.add("STATES_PER_IRREP", new ArrayType());
    /*- Whether to compute all relaxed excited states -*/
    options.add_bool("PROP_ALL",false);
    /*- The symmetry of states -*/
    options.add_int("PROP_SYM", 0);
    /*- -*/
    options.add_int("PROP_ROOT", 0);
  }
  if(name == "CCLAMBDA"|| options.read_globals()) {
    /*- -*/
    options.add_str("WFN","SCF");
    /*- -*/
    options.add_int("CONVERGENCE",7);
    /*- -*/
    options.add_bool("RESTART",false);
    /*- -*/
    options.add_int("CACHELEV",2);
    /*- -*/
    options.add_bool("SEKINO",false);
    /*- -*/
    options.add_bool("DIIS",true);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- -*/
    options.add_str("ABCD","NEW");
    /*- -*/
    options.add_int("NUM_AMPS",10);
    /*- -*/
    options.add_str("JOBTYPE","");
    /*- -*/
    options.add_bool("LOCAL",false);
    /*- -*/
    options.add_str("LOCAL_WEAKP","NONE");
    /*- -*/
    options.add_double("LOCAL_CUTOFF",0.02);
    /*- -*/
    options.add_str("LOCAL_METHOD","WERNER");
    /*- -*/
    options.add_bool("LOCAL_FILTER_SINGLES",true);
    /*- -*/
    options.add_double("LOCAL_CPHF_CUTOFF",0.10);
    /*- -*/
    options.add_str("LOCAL_PAIRDEF","");
    /*- -*/
    options.add("STATES_PER_IRREP", new ArrayType());
    /*- -*/
    options.add_bool("PROP_ALL",false);
    /*- -*/
    options.add_int("PROP_SYM",1);
    /*- -*/
    options.add_int("PROP_ROOT",1);
    /*- -*/
    options.add_int("MAXITER",50);
  }
  if(name == "CLAG"|| options.read_globals()) {
    /*- -*/
    options.add_bool("WRITE_CAS_FILES",0);
    /*- -*/
    options.add_str("WFN","NONE");
    /*- -*/
    options.add_int("ROOT",1);
  }
  if(name == "STABLE"|| options.read_globals()) {
    /*- -*/
    options.add_int("CACHELEV",2);
    /*- -*/
    options.add_str("REFERENCE","RHF");
    /*- -*/
    options.add_bool("FOLLOW",false);
    /*- -*/
    options.add_int("NUM_EVECS_PRINT",0);
    /*- -*/
    options.add_int("ROTATION_METHOD",0);
    /*- -*/
    options.add_double("SCALE",0.5);
  }
  if(name == "OEPROP"|| options.read_globals()) {
    /*- -*/
    options.add_int("NUM_ROOTS",1);
    /*- -*/
    options.add_int("ROOT",1);
    /*- -*/
    options.add_int("GRID",0);
    /*- -*/
    options.add_str("MO_TO_PLOT","");
    /*- -*/
    options.add_int("GRID_ORIGIN",0);
    /*- -*/
    options.add_int("GRID_UNIT_X",0);
    /*- -*/
    options.add("GRID_XY0", new ArrayType());
    /*- -*/
    options.add("GRID_XY1", new ArrayType());
    /*- -*/
    options.add("GRID_XYZ0", new ArrayType());
    /*- -*/
    options.add("GRID_XYZ1", new ArrayType());
    /*- -*/
    options.add_int("NIX",0);
    /*- -*/
    options.add_int("NIY",0);
    /*- -*/
    options.add_int("NIZ",0);
    /*- -*/
    options.add_str("GRID_FORMAT","");
    /*- -*/
    options.add_double("GRID_ZMIN",0.0);
    /*- -*/
    options.add_double("GRID_ZMAX",3.0);
    /*- -*/
    options.add_int("EDGRAD_LOGSCALE",5);
    /*- -*/
    options.add_str("WFN","");
    /*- -*/
    options.add_bool("TRANSITION_DENSITY", false);
    /*- -*/
    options.add_str("REFERENCE", "RHF");
    /*- -*/
    options.add_bool("READ_OPDM", true);
    /*- -*/
    options.add_int("OPDM_FILE", 0);
    /*- -*/
    options.add_str("OPDM_BASIS", "MO", "AO MO");
    /*- -*/
    options.add_str("OPDM_FORMAT", "SQUARE");
    /*- -*/
    options.add_bool("WRTNOS", false);
    /*- -*/
    options.add_bool("ASYMM_OPDM", false);
    /*- -*/
    options.add_bool("SPIN_PROP", false);
    /*- -*/
    options.add_bool("PRINT_NOS", false);
    /*- -*/
    options.add_int("CORREL_CORR", 0);
    /*- -*/
    options.add_double("ZVEC_FILE", 0);
    /*- -*/
    options.add_int("DELETE_ZVEC", 0);
    /*- -*/
    options.add_int("MPMAX", 1);
    /*- -*/
    options.add("MP_REF_XYZ", new ArrayType());
    /*- -*/
    options.add_int("MP_REF", 0);
    /*- -*/
    options.add("LM_REF_XYZ", new ArrayType());
    /*- -*/
    options.add_bool("NUC_ESP", true);
    /*- -*/
    options.add_double("FINE_STRUCTURE_ALPHA", 1.0);
    /*- -*/
    options.add_bool("QED_DARWIN", false);
  }
  if(name == "CCHBAR"|| options.read_globals()) {
    /*- -*/
    options.add_bool("TAMPLITUDE",false);
    /*- -*/
    options.add_int("CACHELEV",2);
    /*- -*/
    options.add_str("WFN", "SCF");
    /*- -*/
    options.add_bool("WABEI_LOWDISK", false);
    /*- -*/
    options.add_str("EOM_REFERENCE","RHF");
  }
  if(name == "CCRESPONSE"|| options.read_globals()) {
    /*- -*/
    options.add_str("WFN", "SCF");
    /*- -*/
    options.add_int("CACHELEV",2);
    /*- -*/
    options.add_str("REFERENCE","RHF");
    /*- -*/
    options.add_str("GAUGE","LENGTH");
    /*- -*/
    options.add_int("MAXITER",50);
    /*- -*/
    options.add_int("CONVERGENCE",7);
    /*- -*/
    options.add_bool("DIIS",1);
    /*- -*/
    options.add_str("PROPERTY","POLARIZABILITY");
    /*- -*/
    options.add_str("ABCD","NEW");
    /*- -*/
    options.add_bool("RESTART",1);
    /*- -*/
    options.add_bool("LOCAL",0);
    /*- -*/
    options.add_double("LOCAL_CUTOFF",0.01);
    /*- -*/
    options.add_str("LOCAL_METHOD","WERNER");
    /*- -*/
    options.add_str("LOCAL_WEAKP","NONE");
    /*- -*/
    options.add_bool("LOCAL_FILER_SINGLES", false);
    /*- -*/
    options.add_double("LOCAL_CPHF_CUTOFF",0.10);
    /*- -*/
    options.add_str("LOCAL_PAIRDEF","NONE");
    /*- -*/
    options.add_bool("ANALYZE",0);
    /*- -*/
    options.add_int("NUM_AMPS",5);
    /*- -*/
    options.add_bool("SEKINO",0);
    /*- -*/
    options.add_bool("LINEAR",0);
    /*- -*/
    options.add("OMEGA",new ArrayType());
    /*- -*/
    options.add("MU_IRREPS",new ArrayType());
  }
  if(name == "MVO"|| options.read_globals()) {
    /*- -*/
   options.add_str("WFN","CCSD");
    /*- -*/
   options.add_int("FZC_FILE", PSIF_OEI);
    /*- -*/
   options.add_bool("PRINT_MOS",false);
    /*- -*/
   options.add_bool("OEI_ERASE",false);
    /*- -*/
   options.add_bool("FZC",true);
    /*- -*/
   options.add_bool("DELETE_RESTR_DOCC",true);
    /*- -*/
   options.add_bool("MP2NOS",false);
    /*- -*/
   options.add_bool("UNOS",false);
    /*- -*/
   options.add_double("FZC_FOCK_COEFF",1.0);
    /*- -*/
   options.add_double("FOCK_COEFF",0.0);
    /*- -*/
   options.add_bool("IVO",false);
    /*- -*/
   options.add_bool("CANONICAL",false);
    /*- -*/
   options.add("RESTRICTED_DOCC", new ArrayType());
    /*- -*/
   options.add("RESTRICTED_UOCC", new ArrayType());
    /*- -*/
   options.add("DOCC_VIRT", new ArrayType());
  }
  if(name == "RESPONSE"|| options.read_globals()){
    /*- -*/
    options.add_str("REFERENCE", "RHF");
    /*- -*/
    options.add("OMEGA", new ArrayType());
    /*- -*/
    options.add_str("PROPERTY","POLARIZABILITY");
  }
  if(name == "MCSCF"|| options.read_globals()) {
    /*- The molecular charge -*/
    options.add_int("CHARGE", 0);
    /*- (2$\times M_s+1$), e.g. 1 for a singlet state, 2 for a doublet, 3 for a triplet, etc. -*/
    options.add_int("MULTP", 1);
    /*- -*/
    options.add_int("CONVERGENCE",9);
    /*- Level shift to aid convergence -*/
    options.add_int("LEVELSHIFT",0);
    /*- The amount of debugging information to print -*/
    options.add_int("DEBUG", false);
    /*- -Log10 of the energy convergence criterion -*/
    options.add_int("E_CONVERGE", 12);
    /*- -Log10 of the density convergence criterion -*/
    options.add_int("D_CONVERGE", 12);
    /*- -*/
    options.add_int("MAXITER",100);
    /*- -*/
    options.add_int("NDIIS",7);
    /*- -*/
    options.add_int("ROOT",1);
    /*- -*/
    options.add_int("START_FAVG",5);
    /*- -*/
    options.add_int("TURN_ON_ACTV",0);
    /*- -*/
    options.add_int("ROTATE_MO_ANGLE",0);
    /*- -*/
    options.add_int("ROTATE_MO_IRREP",1);  // IRREP is one-based
    /*- -*/
    options.add_int("ROTATE_MO_P",1);      // P and Q are one-based
    /*- -*/
    options.add_int("ROTATE_MO_Q",2);
    /*- -*/
    options.add_bool("CI_DIIS",false);
    /*- -*/
    options.add_bool("USE_DIIS",true);
    /*- -*/
    options.add_bool("READ_MOS",true);
    /*- -*/
    options.add_bool("USE_FAVG",false);
    /*- -*/
    options.add_bool("CANONICALIZE_ACTIVE_FAVG",false);
    /*- -*/
    options.add_bool("CANONICALIZE_INACTIVE_FAVG",false);
    /*- -*/
    options.add_bool("INTERNAL_ROTATIONS",true);
    /*- -*/
    options.add_bool("FORCE_TWOCON",false);
    /*- The number of active orbitals, per irrep -*/
    options.add("ACTIVE", new ArrayType());
    /*- The number of active orbitals, per irrep (alternative name for ACTIVE) -*/
    options.add("ACTV", new ArrayType());


    /*- -*/
    options.add_str("REFERENCE","RHF","RHF ROHF UHF TWOCON MCSCF GENERAL");
    /*- -*/
    options.add_str("WFN_SYM","1","A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");
  }
  if(name == "CCENERGY"|| options.read_globals()) {
    /*- -*/
    options.add_bool("NEWTRIPS", 1);
    /*- -*/
    options.add_str("WFN", "NONE", "CCSD CCSD_T EOM_CCSD LEOM_CCSD BCCD BCCD_T CC2 CC3 EOM_CC2 EOM_CC3 CCSD_MVD");
    /*- -*/
    options.add_str("REFERENCE", "RHF");
    /*- -*/
    options.add_bool("ANALYZE", 0);
    /*- -*/
    options.add_int("MAXITER", 50);
    /*- -*/
    options.add_int("CONVERGENCE", 7);
    /*- -*/
    options.add_bool("RESTART",1);
    /*- -*/
    options.add_bool("FORCE_RESTART", 0);
//#warning CCEnergy ao_basis keyword type was changed.
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- -*/
    options.add_int("CACHELEV", 2);
    /*- -*/
    options.add_str("CACHETYPE", "LOW", "LOW LRU");
    /*- -*/
    options.add_int("NTHREADS",1);
    /*- -*/
    options.add_bool("DIIS", true);
    /*- -*/
    options.add_bool("T2_COUPLED", false);
    /*- -*/
    options.add_str("PROPERTY", "POLARIZABILITY", "POLARIZABILITY ROTATION MAGNETIZABILITY ROA ALL");
    /*- -*/
    options.add_str("ABCD", "NEW", "NEW OLD");
    /*- -*/
    options.add_bool("LOCAL", 0);
    /*- -*/
    options.add_double("LOCAL_CUTOFF", 0.02);
    /*- -*/
    options.add_double("LOCAL_MOS", 0);
    /*- -*/
    options.add_str("LOCAL_METHOD", "WERNER", "WERNER AOBASIS");
    /*- -*/
    options.add_str("LOCAL_WEAKP", "NONE", "NONE NEGLECT MP2");
    //options.add_int("LOCAL_FILTER_SINGLES", 1);
    options.add_double("LOCAL_CPHF_CUTOFF", 0.10);
    /*- -*/
    options.add_str("LOCAL_PAIRDEF", "BP", "BP RESPONSE");
    /*- -*/
    options.add_int("NUM_AMPS", 10);
    /*- -*/
    options.add_int("BRUECKNER_CONV", 5);
    /*- -*/
    options.add_bool("PRINT_MP2_AMPS", 0);
    /*- -*/
    options.add_bool("PRINT_PAIR_ENERGIES", 0);
    /*- -*/
    options.add_bool("SPINADAPT_ENERGIES", false);
    /*- -*/
    options.add_bool("T3_WS_INCORE", 0);
    /*- -*/
    options.add_bool("SCSN_MP2", 0);
    /*- -*/
    options.add_bool("SCS_MP2", 0);
    /*- -*/
    options.add_bool("SCS_CCSD", 0);
    /*- -*/
    options.add_double("MP2_SCALE_OS",1.20);
    /*- -*/
    options.add_double("MP2_SCALE_SS",1.0/3.0);
    /*- -*/
    options.add_double("CC_SCALE_OS", 1.27);
    /*- -*/
    options.add_double("CC_SCALE_SS",1.13);
  }
  if(name == "CIS"|| options.read_globals()) {
    /*- -*/
    options.add_str("WFN", "CIS", "CCSD CCSD_T EOM_CCSD CIS");
    /*- -*/
    options.add_str("REFERENCE", "RHF", "RHF ROHF UHF");
    /*- -*/
    options.add_double("LOCAL_AMP_PRINT_CUTOFF", 0.60);
    /*- -*/
    options.add_int("MAXITER", 500);
    /*- -*/
    options.add_int("CONVERGENCE", 7);
    /*- -*/
    options.add("STATES_PER_IRREP", new ArrayType());
    /*- -*/
    options.add_str("DIAG_METHOD", "DAVIDSON", "DAVIDSON FULL");
    /*- -*/
    options.add_bool("LOCAL", false);
    /*- -*/
    options.add_double("LOCAL_CUTOFF", 0.02);
    /*- -*/
    options.add_str("LOCAL_METHOD", "WERNER", "AOBASIS WERNER");
    /*- -*/
    options.add_str("LOCAL_WEAKP", "MP2", "MP2 NEGLECT NONE");
    /*- -*/
    options.add_int("LOCAL_GHOST", -1);
    /*- -*/
    options.add("DOMAINS", new ArrayType());
    /*- -*/
    options.add_bool("DOMAIN_PRINT", 0);
  }
  if(name == "LMP2"|| options.read_globals()) {
    /*- The wavefunction desired -*/
    options.add_str("RI_BASIS_MP2", "");
//    options.read_ipv1();
    /*- -*/
    if(options.get_str("RI_BASIS_MP2") != "")
    /*- -*/
      options.add_bool("RI_LMP2", true);
    else
    /*- -*/
      options.add_bool("RI_LMP2", false);
    /*- -*/
    options.add_str("WFN", "LMP2");
    /*- -*/
    options.add_str("REFERENCE", "RHF", "RHF");
    /*- -*/
    options.add_int("MAXITER", 50);
    /*- -*/
    options.add_int("ENERGY_CONV", 7);
    /*- -*/
    options.add_int("RMS_CONV", 5);
    /*- -*/
    options.add_int("FSKIP", 2);
    /*- -*/
    options.add_bool("USE_DIIS", 1);
    /*- -*/
    options.add_bool("NEGLECT_DP", 1);
    /*- -*/
    options.add_double("DISTANT_PAIR", 8.0);
    /*- -*/
    options.add_int("DIISSTART", 3);
    /*- -*/
    options.add_int("NDIIS", 5);
    /*- -*/
    options.add_double("LOCAL_CUTOFF", 0.02);
    /*- -*/
    options.add_int("MEMORY", 2000);
    /*- -*/
    options.add_bool("SCS","false");
    /*- -*/
    options.add_bool("SCS_N", "false");
    /*- -*/
    options.add_double("SCALE_OS", 6.0/5.0);
    /*- -*/
    options.add_double("SCALE_SS", 1.0/3.0);
    /*- -*/
    options.add_int("SCREENING", 7);
    /*- -*/
    options.add_bool("SCREEN_INTS", false);
    /*- -*/
    options.add_int("SCHWARTZ_TOL", 12);
   }
  if(name=="DFMP2"|| options.read_globals()) {
    options.add_int("MADMP2_SLEEP", 0);
    options.add_int("MADMP2_DEBUG", 0);
    options.add_str("MP2_ALGORITHM", "DFMP2");
    //options.add_str("WFN", "RI-MP2");
    /*- RI Basis, needed by Python -*/
    options.add_str("RI_BASIS_MP2","");
    /*- Basis, needed by Python -*/
    options.add_str("BASIS","NONE");
    /*- OS Scale  !expert -*/
    options.add_double("SCALE_OS", 6.0/5.0);
    /*- SS Scale  -*/
    options.add_double("SCALE_SS", 1.0/3.0);
    /*- % of memory for DF-MP2 three-index buffers  -*/
    options.add_double("DFMP2_MEM_FACTOR", 0.9);
    /*- Schwarz cutoff -*/
    options.add_double("SCHWARZ_CUTOFF", 0.0);
    /*- DFMP2 Fitting Type -*/
    options.add_str("RI_FITTING_TYPE", "FINISHED", "FINISHED RAW CHOLESKY");
    /*- DFMP2 Algorithm (usually for debugging)  -*/
    options.add_str("DFMP2_TYPE","DEFAULT", "DEFAULT DISK CORE OLD");
    /*- DFMP2 Fitting symmetry  -*/
    options.add_str("FITTING_SYMMETRY","SYMMETRIC", "SYMMETRIC ASYMMETRIC");
    /*- DFMP2 Fitting conditioning  -*/
    options.add_str("FITTING_CONDITIONING","FINISHED", "RAW FINISHED");
    /*- DFMP2 Fitting inversion  -*/
    options.add_str("FITTING_INVERSION","EIG", "EIG CHOLESKY SOLVE");
    /*- Max condition number in auxiliary basis -*/
    options.add_double("RI_MAX_COND", 1.0E8);
    /*- -Maximum number of rows to read/write in each DF-MP2 operation -*/
    options.add_int("ROWS_PER_READ", 0);
    /*- Number of threads to compute integrals with. 0 is wild card -*/
    options.add_int("RI_INTS_NUM_THREADS", 0);
    /*- Debugging information? -*/
    options.add_int("DEBUG",0);
    /*- Parallel algoritmh? -*/
    options.add_bool("PARALLEL_DFMP2",false);
    /*- -Log10 of the energy convergence criterion -*/
    options.add_int("E_CONVERGE", 8);
    /*- -Log10 of the density convergence criterion -*/
    options.add_int("D_CONVERGE", 8);
  }
  if(name=="DFCC"|| options.read_globals()) {
    /*- Type of wavefunction -*/
    options.add_str("WAVEFUNCTION","MP2","MP2 MP3 CCD DRPA");
    /*- MO basis -*/
    options.add_str("BASIS","NONE");
    /*- Schwarz cutoff -*/
    options.add_double("SCHWARZ_CUTOFF", 0.0);
    /*- Convergence of CC energy -*/
    options.add_int("E_CONVERGE", 8);
    /*- Convergence of cluster amplitudes (RMS change) -*/
    options.add_int("T_CONVERGE", 8);
    /*- Turn on DIIS -*/
    options.add_bool("DIIS",true);
    /*- Minimum DIIS vectors -*/
    options.add_int("MIN_DIIS_VECS", 2);
    /*- Maximum DIIS vectors -*/
    options.add_int("MAX_DIIS_VECS", 6);
    /*- The maximum number iterations allowed -*/
    options.add_int("MAXITER", 40);
    /*- Debugging information? -*/
    options.add_int("DEBUG",0);

    // => DF <= //

    /*- DF basis for MO integrals -*/
    options.add_str("RI_BASIS_CC","NONE");
    /*- Fitting metric algorithm -*/
    options.add_str("FITTING_TYPE", "EIG", "EIG CHOLESKY QR");
    /*- Desired Fitting condition (inverse of max condition number) -*/
    options.add_double("FITTING_CONDITION", 1.0E-10);

    // => PS <= //

    /*- Dealias basis for PS integrals -*/
    options.add_str("DEALIAS_BASIS_CC","");
    /*- Dealias basis beta parameter -*/
    options.add_double("DEALIAS_BETA", 3.5);
    /*- Dealias basis delta parameter -*/
    options.add_double("DEALIAS_DELTA", 2.0);
    /*- Dealias basis N core parameter -*/
    options.add_int("DEALIAS_N_CORE", 1);
    /*- Dealias basis N intercalater parameter -*/
    options.add_int("DEALIAS_N_INTERCALATER", 1);
    /*- Dealias basis N diffuse parameter -*/
    options.add_int("DEALIAS_N_DIFFUSE", 1);
    /*- Dealias basis N cap parameter -*/
    options.add_int("DEALIAS_N_CAP", 1);
    /*- Dealias basis highest delta l parameter -*/
    options.add_int("DEALIAS_N_L", 1);
    /*- Filename to read grid from -*/
    options.add_str_i("PS_GRID_FILE","");
    /*- File path to read grids from -*/
    options.add_str_i("PS_GRID_PATH","");
    /*- Maximum order of spherical grids -*/
    options.add_int("PS_ORDER_SPHERICAL", 7);
    /*- Number of radial points -*/
    options.add_int("PS_N_RADIAL", 5);
    /*- Spherical Scheme -*/
    options.add_str("PS_SPHERICAL_SCHEME", "LEBEDEV", "LEBEDEV");
    /*- Radial Scheme -*/
    options.add_str("PS_RADIAL_SCHEME", "TREUTLER", "TREUTLER BECKE MULTIEXP EM MURA");
    /*- Nuclear Scheme -*/
    options.add_str("PS_NUCLEAR_SCHEME", "TREUTLER", "TREUTLER BECKE NAIVE STRATMANN");
    /*- Pruning Scheme -*/
    options.add_str("PS_PRUNING_SCHEME", "FLAT", "FLAT P_GAUSSIAN D_GAUSSIAN P_SLATER D_SLATER LOG_GAUSSIAN LOG_SLATER");
    /*- Factor for effective BS radius in radial grid -*/
    options.add_double("PS_BS_RADIUS_ALPHA",1.0);
    /*- Spread alpha for logarithmic pruning -*/
    options.add_double("PS_PRUNING_ALPHA",1.0);
    /*- The number of grid points per evaluation block -*/
    options.add_int("PS_MAX_POINTS",5000);
    /*- The number of grid points per evaluation block -*/
    options.add_int("PS_MIN_POINTS",0);
    /*- The DFT basis cutoff -*/
    options.add_double("PS_BASIS_CUTOFF", 0.0);
    /*- Minumum eigenvalue for primary basis -*/
    options.add_double("PS_MIN_S_PRIMARY",1.0E-7);
    /*- Minumum eigenvalue for dealias basis -*/
    options.add_double("PS_MIN_S_DEALIAS",1.0E-7);
    /*- Fitting algorithm to use for pseudospectral -*/
    options.add_str("PS_FITTING_ALGORITHM", "CONDITIONED", "DEALIASED RENORMALIZED QUADRATURE");
    /*- Pseudospectral range-separation parameter -*/
    options.add_double("PS_OMEGA", 1.0);
    /*- Pseudospectral partition alpha -*/
    options.add_double("PS_ALPHA", 1.0);
    /*- Use range-separation procedure in PS? -*/
    options.add_bool("PS_USE_OMEGA", true);

    // => DENOMINATOR <= //

    /*- Denominator algorithm for PT methods -*/
    options.add_str("DENOMINATOR_ALGORITHM", "LAPLACE", "LAPLACE CHOLESKY");
    /*- Maximum denominator error allowed (Max error norm in Delta tensor) -*/
    options.add_double("DENOMINATOR_DELTA", 1.0E-6);

    // => MP2 <= //

    /*- MP2 Algorithm:
            \begin{tabular}{ccc}
            Algorithm Keyword  &  MP2J        &  MP2K  \\
             \hline
                MP2            &   MP2        &  MP2   \\
                DF             &   DF         &  DF    \\
                PS             &   PS         &  PS    \\
                PS1            &   DF         &  PS/DF \\
                PS2            &   DF         &  PS/PS \\
                PS3            &   PS         &  PS/DF \\
                PS4            &   PS         &  PS/PS \\
                TEST_DENOM     &   Test       &  Test  \\
                TEST_PS        &   Test       &  Test  \\
                TEST_PS_OMEGA  &   Test       &  Test  \\
                TEST_DPS_OMEGA &   Test       &  Test  \\
                TEST_DF        &   Test       &  Test  \\
            \end{tabular}
    -*/
    options.add_str("MP2_ALGORITHM", "DF", "MP2 DF PS PS1 PS2 PS3 PS4 TEST_DENOM TEST_PS TEST_PS_OMEGA TEST_DPS_OMEGA TEST_DF");
    /*- OS Scale  -*/
    options.add_double("SCALE_OS", 6.0/5.0);
    /*- SS Scale  -*/
    options.add_double("SCALE_SS", 1.0/3.0);

    // => RPA <= //

    /*- RPA algorithm:
        \begin{tabular}{cc}
        DF & $\mathcal{O}(N^5)$ \\
        CD & $\mathcal{O}(N^4)$ \\
        \end{tabular}
    -*/
    options.add_str("RPA_ALGORITHM", "CD", "CD DF");
    /*- RPA Cholesky delta -*/
    options.add_double("RPA_DELTA", 1.0E-6);
    /*- Continue RPA even if T's are not numerically SPD? -*/
    options.add_bool("RPA_RISKY",false);
    /*- Continue RPA numerical SPD Tolerance -*/
    options.add_double("RPA_PLUS_EPSILON",1.0E-12);
    /*- RPA alpha parameter -*/
    options.add_double("RPA_ALPHA", 1.0);

  }
  if(name == "PSIMRCC"|| options.read_globals()) {
    /*- The molecular charge of the target state -*/
    options.add_int("CORR_CHARGE",0);
    /*- Amount of debugging output to produce !expert -*/
    options.add_int("DEBUG",0);
    /*- The amount of damping to apply to the amplitude updates.  If this is set to 0, the full update is performed.
        A value of 1000 retains the amplitudes from the previous iteration always being used.  A value in between these extremes
        can help in cases where oscillatory convergence is observed -*/
    options.add_int("DAMPING_FACTOR",0);
    /*- The number of DIIS vectors to use in extrapolations -*/
    options.add_int("MAXDIIS",7);
    /*- The number of threads to use on multi-core machines -*/
    options.add_int("NUM_THREADS",1);
    /*- The number of electrons -*/
    options.add_int("NEL",0);
    /*- Which root of the effective hamiltonian is the target state? -*/
    options.add_int("ROOT",1);
    /*- The number of digits after the decimal to converge the energy to -*/
    options.add_int("CONVERGENCE",9);
    /*- The maximum number of iterations allowed to determine the amplitudes -*/
    options.add_int("MAXITER",100);
    /*- The number of DIIS vectors needed before extrapolation is performed -*/
    options.add_int("START_DIIS",2);
    /*- The shift to apply to the denominators ($\times$ 1000), {\it c.f.} Taube and Bartlett, JCP, 130, 144112 (2009) -*/
    options.add_int("TIKHONOW_OMEGA",0);  // Omega = TIKHONOW_OMEGA / 1000
    /*- The cycle after which Tikhonow regularization is stopped.  Set to zero to allow regularization in all iterations -*/
    options.add_int("TIKHONOW_MAX",5);
    /*- Whether to perform DIIS extrapolation for iterative triple excitation computations -*/
    options.add_bool("DIIS_TRIPLES",false);
    /*- Whether to lock onto a singlet root -*/
    options.add_bool("LOCK_SINGLET",false);
    /*- Whether to start from an MP2 guess -*/
    options.add_bool("MP2_GUESS",true);
    /*- Whether to use the averaged Fock matrix over all references in (T) computations -*/
    options.add_bool("FAVG_CCSD_T",false);
    /*- Whether to include the fourth-order contributions to the effective hamiltonian -*/
    options.add_bool("HEFF4",true);
    /*- Whether to include the off-diagonal corrections in (T) computations -*/
    options.add_bool("OFFDIAGONAL_CCSD_T",true);
    /*- Whether to include the diagonal corrections in (T) computations -*/
    options.add_bool("DIAGONAL_CCSD_T",true);
    /*- Whether to diagonalize the effective Hamiltonian -*/
    options.add_bool("DIAGONALIZE_HEFF",false);
    /*- Whether to perform DIIS extrapolation -*/
    options.add_bool("USE_DIIS",true);
    /*- Whether to use symmetry to map equivalent determinants onto each other, for efficiency -*/
    options.add_bool("USE_SPIN_SYMMETRY",true);
    /*- Whether to zero the internal amplitudes, i.e., those that map reference determinants onto each other -*/
    options.add_bool("ZERO_INTERNAL_AMPS",true);
    /*- Whether to include the terms that couple the reference determinants -*/
    options.add_bool("COUPLING_TERMS",true);
    /*- Whether to print the effective Hamiltonian -*/
    options.add_bool("PRINT_HEFF",false);
    /*- Whether to compute the perturbative corrections for basis set incompleteness !expert -*/
    options.add_bool("PERT_CBS",false);
    /*- Whether to include the terms that couple different reference determinants in
        perturbative CBS correction computations !expert -*/
    options.add_bool("PERT_CBS_COUPLING",true);
    /*- Whether to use Tikhonow regularization in (T) computations !expert -*/
    options.add_bool("TIKHONOW_TRIPLES",false);
    /*- The type of perturbation theory computation to perform -*/
    options.add_str("PT_ENERGY","SECOND_ORDER","SECOND_ORDER SCS_SECOND_ORDER PSEUDO_SECOND_ORDER SCS_PSEUDO_SECOND_ORDER");
    /*- The type of correlated wavefunction -*/
    options.add_str("CORR_WFN","CCSD","PT2 CCSD MP2-CCSD CCSD_T");
    /*- The type of CCSD(T) computation to perform -*/
    options.add_str("CORR_CCSD_T","STANDARD","STANDARD PITTNER");
    /*- The type of reference function used in MRCC computations -*/
    options.add_str("CORR_REFERENCE","GENERAL","RHF ROHF TCSCF MCSCF GENERAL");
    /*- The ansatz to use for MRCC computations -*/
    options.add_str("CORR_ANSATZ","MK","SR MK BW APBW");
    /*- The order of coupling terms to include in MRCCSDT computations -*/
    options.add_str("COUPLING","CUBIC","NONE LINEAR QUADRATIC CUBIC");
    /*- The symmetry of the target wavefunction, specified either by Sch\"onflies symbol,
        or irrep number (in Cotton ordering) -*/
    options.add_str("WFN_SYM","1","A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");
    /*- The type of algorithm to use for (T) computations -*/
    options.add_str("TRIPLES_ALGORITHM","RESTRICTED","SPIN_ADAPTED RESTRICTED UNRESTRICTED");
    /*- How to perform MP2_CCSD computations -*/
    options.add_str("MP2_CCSD_METHOD","II","I IA II");
    /*- The number of frozen occupied orbitals per irrep (same as FROZEN_DOCC)-*/
    options.add("CORR_FOCC", new ArrayType());
    /*- The number of doubly occupied orbitals per irrep (same as RESTRICTED_DOCC)-*/
    options.add("CORR_DOCC", new ArrayType());
    /*- The number of doubly occupied orbitals per irrep (same as CORR_DOCC) -*/
    options.add("RESTRICTED_DOCC", new ArrayType());
    /*- The number of active orbitals per irrep (same as ACTV) -*/
    options.add("CORR_ACTV", new ArrayType());
    /*- The number of active orbitals per irrep (same as CORR_ACTV, ACTIVE) -*/
    options.add("ACTV", new ArrayType());
    /*- The number of active orbitals per irrep (same as CORR_ACTV, ACTV) -*/
    options.add("ACTIVE", new ArrayType());
    /*- The number of frozen virtual orbitals (same as FROZEN_UOCC) -*/
    options.add("CORR_FVIR", new ArrayType());
//    /*- The number of -*/
//    options.add("ACTIVE_DOCC", new ArrayType());
  }
  if(name == "OPTKING"|| options.read_globals()) {
      /*- Specifies minimum search, transition-state search, or IRC following; allowed values = {MIN, TS, IRC} -*/
      options.add_str("OPT_TYPE", "MIN", "MIN TS IRC");
      /*- Whether to do a Newton-Raphson step, or an RFO step; allowed values = {NR, RFO} -*/
      options.add_str("STEP_TYPE", "RFO", "RFO NR");
      /*- Maximum step size in bohr or radian along an internal coordinate {double} -*/
      options.add_double("INTRAFRAGMENT_STEP_LIMIT", 0.4);
      /*- Whether to 'follow' the initial RFO vector after the first step {true, false} -*/
      options.add_bool("RFO_FOLLOW_ROOT", false);
      /*- Which RFO root to follow; 0 indicates lowest (to a minimum); {integer} -*/
      options.add_int("RFO_ROOT", 0);
      /*- When determining connectivity, a bond is assigned if interatomic distance
          is less than (this number) * sum of covalent radii {double} -*/
      options.add_double("SCALE_CONNECTIVITY", 1.3);
      /*- Whether to treat multiple molecule fragments as a single bonded molecule;
          or via interfragment coordinates ; a primary difference is that in MULTI mode,
          the interfragment coordinates are not redundant. {SINGLE, MULTI} -*/
      options.add_str("FRAGMENT_MODE", "SINGLE", "SINGLE MULTI");
      /*- whether to use fixed linear combinations of atoms as reference points for
          interfragment coordinates or whether to use principal axes {FIXED, PRINCIPAL_AXES} -*/
      options.add_str("INTERFRAGMENT_MODE", "FIXED", "FIXED INTERFRAGMENT");
      /*- Whether to only generate the internal coordinates and then stop {true, false} -*/
      options.add_bool("GENERATE_INTCOS_ONLY", false);
      /*- What model Hessian to use to guess intrafragment force constants {SCHLEGEL, FISCHER} -*/
      options.add_str("INTRAFRAGMENT_H", "FISCHER", "FISCHER SCHLEGEL");
      /*- Whether to use the default of FISCHER_LIKE force constants for the initial guess {DEFAULT, FISCHER_LIKE} -*/
      options.add_str("INTERFRAGMENT_H", "DEFAULT", "DEFAULT FISCHER_LIKE");
      /*- Whether to freeze all fragments rigid -*/
      options.add_bool("FREEZE_INTRAFRAGMENT", false);
      /*- By default, optking prints and saves the last (previous) geometry at the end of an
          optimization, i.e., the one at which a gradient was computed.  If this keyword is
          set to true, then the structure obtained from the last projected step is printed out and saved instead. -*/
      options.add_bool("WRITE_FINAL_STEP_GEOMETRY", false);
      /*- Choose from supported Hessian updates {NONE, BFGS, MS, POWELL, BOFILL} -*/
      options.add_str("H_UPDATE", "BFGS", "NONE BFGS MS POWELL BOFILL");
      /*-  How many previous steps' data to use in Hessian update; 0=use them all ; {integer} -*/
      options.add_int("H_UPDATE_USE_LAST", 6);
      /*- Whether to limit the magnitutde of changes caused by the Hessian update {true, false} -*/
      options.add_bool("H_UPDATE_LIMIT", true);
      /*- If the above is true, changes to the Hessian from the update are limited to the larger of
          (H_update_limit_scale)*(the previous value) and H_update_limit_max (in au). -*/
      options.add_double("H_UPDATE_LIMIT_MAX", 1.00);
      /*- If the above is true, changes to the Hessian from the update are limited to the larger of
          (H_update_limit_scale)*(the previous value) and H_update_limit_max (in au). -*/
      options.add_double("H_UPDATE_LIMIT_SCALE", 0.50);
      /*- Whether to use 1/R(AB) for stretching coordinate between fragments (or just R(AB)) -*/
      options.add_bool("INTERFRAGMENT_DISTANCE_INVERSE", false);
      /*- For now, this is a general maximum distance for the definition of H-bonds -*/
      options.add_double("MAXIMUM_H_BOND_DISTANCE", 4.3);
      /*- QCHEM optimization criteria: maximum force -*/
      options.add_double("CONV_MAX_FORCE", 3.0e-4);
      /*- QCHEM optimization criteria: maximum energy change -*/
      options.add_double("CONV_MAX_DE", 1.0e-6);
      /*- QCHEM optimization criteria: maximum displacement -*/
      options.add_double("CONV_MAX_DISP", 1.2e-3);
      /*- Whether to test B matrix -*/
      options.add_bool("TEST_B", false);
      /*- Whether to test derivative B matrix -*/
      options.add_bool("TEST_DERIVATIVE_B", false);
      /*- Read Cartesian Hessian -*/
      options.add_bool("READ_CARTESIAN_H", false);
      /*- Define IRC step size in bohr(amu)$^1/2$ -*/
      options.add_double("IRC_STEP_SIZE", 0.2);
      /*- Define IRC mapping direction {FORWARD, BACKWARD} -*/
      options.add_str("IRC_DIRECTION", "FORWARD", "FORWARD BACKWARD");
  }
  if(name == "FINDIF"|| options.read_globals()) {
      /*- Number of points for finite-differences (3 or 5) -*/
      options.add_int("POINTS", 3); // Can we error check integers?
      /*- Displacement size in au for finite-differences. -*/
      options.add_double("DISP_SIZE", 0.005);
  }
  return true;
}

} //end ::psi

