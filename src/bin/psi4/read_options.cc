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

  /*- An array containing the number of doubly-occupied orbitals per irrep
  (in Cotton order) -*/
  options.add("DOCC", new ArrayType());
  /*- An array containing the number of singly-occupied orbitals per irrep
  (in Cotton order) -*/
  options.add("SOCC", new ArrayType());
  /*- An array containing the number of frozen doubly-occupied orbitals per
  irrep (these are not excited in a CI, nor can they be optimized in MCSCF -*/
  options.add("FROZEN_DOCC", new ArrayType());
  /*- An array containing the number of frozen unoccupied orbitals per
  irrep (these are not populated in a CI, nor can they be optimized in MCSCF -*/
  options.add("FROZEN_UOCC", new ArrayType());
  /*- The scope of core orbitals to freeze in later correlated computations
  FROZEN_DOCC trumps this option -*/
  options.add_int("NUM_FROZEN_DOCC", 0);
  /*- The scope of virtual orbitals to freeze in later correlated computations
  FROZEN_UOCC trumps this option -*/
  options.add_int("NUM_FROZEN_UOCC", 0);
  /*- The scope of core orbitals to freeze in later correlated computations
  FROZEN_DOCC or NUM_FROZEN_DOCC trumps this option  -*/
  options.add_str("FREEZE_CORE","FALSE", "FALSE TRUE SMALL LARGE");

  /*- Do use pure angular momentum basis functions?
  If not explicitly set, the default comes from the basis set. -*/
  options.add_bool("PUREAM", true);
  /*- The amount of information to print to the output file -*/
  options.add_int("PRINT", 1);
  /*- The amount of information to print to the output file !expert -*/
  options.add_int("DEBUG", 0);
  /*- Some codes (DFT) can dump benchmarking data to separate output files -*/
  options.add_int("BENCH", 0);
  /*- Maximum number of geometry optimization steps -*/
  options.add_int("GEOM_MAXITER", 20);
  /*- Wavefunction type !expert -*/
  options.add_str("WFN", "SCF");
  /*- Derivative level !expert -*/
  options.add_str("DERTYPE", "NONE", "NONE FIRST SECOND RESPONSE");
  /*- Number of columns to print in calls to Matrix::print_mat !expert -*/
  options.add_int("MAT_NUM_COLUMN_PRINT", 5);

  // CDS-TODO: We should go through and check that the user hasn't done
  // something silly like specify frozen_docc in DETCI but not in TRANSQT.
  // That would create problems.  (This was formerly checked in DETCI
  // itself, but I don't think DETCI will have the info available to check
  // this anymore).  This problem has affected users in the past.
  // Same goes for restricted_docc, restricted_uocc, ras1, ras2, ras3,
  // frozen_uocc.

  if (name == "DETCI" || options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "DETCI", "DETCI CI ZAPTN DETCAS CASSCF RASSCF");

    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF", "RHF ROHF");

    /*- Convergence criterion for CI residual vector in the Davidson algorithm (RMS error). 
    The default is 1e-4 for energies and 1e-7 for gradients. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("R_CONVERGENCE", 1e-4);

    /*- Convergence criterion for energy. See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("E_CONVERGENCE", 1e-6);

    /*- Do a full CI (FCI)? If TRUE, overrides the value of EX_LEVEL -*/
    options.add_bool("FCI",false);

    /*- The CI excitation level -*/
    options.add_int("EX_LEVEL", 2);

    /*- The CC excitation level -*/
    options.add_int("CC_EX_LEVEL", 2);

    /*- In a RAS CI, this is the additional excitation level for allowing
    electrons out of RAS I into RAS II.  The maximum number of holes in
    RAS I is therefore EX_LEVEL + VAL_EX_LEVEL. -*/
    options.add_int("VAL_EX_LEVEL", 0);

    /*- The CC valence excitation level -*/
    options.add_int("CC_VAL_EX_LEVEL", 0);

    /*- number of CI roots to find -*/
    options.add_int("NUM_ROOTS", 1);

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

    /*- Do use the $M@@s = 0$ component of the state? Defaults to TRUE
    if closed-shell and FALSE otherwise. Related to the S option. -*/
    options.add_bool("MS0",false);

    /*- Do stop DETCI after string information is formed
    and before integrals are read? -*/
    options.add_bool("ISTOP",false);

    /*- Do print a summary of the CI blocks? -*/
    options.add_bool("CIBLKS_PRINT",false);

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
    options.add_int("REFERENCE_SYM", -1);

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

    /*- Do use coupling block in preconditioner? !expert -*/
    options.add_bool("H0_BLOCK_COUPLING",false);

    /*- Parameters which specifies the size of the coupling block
     within the generalized davidson preconditioner. !expert -*/
    options.add_int("H0_BLOCK_COUPLING_SIZE",0);

    /*- Number of important determinants to print -*/
    options.add_int("NUM_PRINT",20);

    /*- Number of important CC amplitudes per excitation level to print -*/
    options.add_int("CC_NUM_PRINT",10);

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

    /*- Do compute the diagonal elements of the Hamiltonian matrix 
      on-the-fly? Otherwise, a diagonal element vector is written
      to a separate file on disk. !expert -*/
    options.add_bool("HD_OTF",true);

    /*- Do use the last vector space in the BVEC file to write
      scratch DVEC rather than using a separate DVEC file? (Only
      possible if NUM_ROOTS = 1.) !expert -*/
    options.add_bool("NO_DFILE",false);

    /*- Do freeze core orbitals? -*/
    // CDS-TODO: Need to make DETCI compatible with normal FREEZE_CORE
    options.add_bool("DETCI_FREEZE_CORE",true);

    /*- Do store strings specifically for FCI? (Defaults to TRUE for FCI.)
        !expert -*/
    options.add_bool("FCI_STRINGS",false);

    /*- Do allow `mixed' RAS II/RAS III excitations into the CI space?
      If FALSE, then if there are any electrons
      in RAS III, then the number of holes in RAS I cannot exceed the given
      excitation level EX_LEVEL. !expert -*/
    options.add_bool("MIXED",true);

    /*- Do allow `mixed' excitations involving RAS IV into the CI space.
      Useful to specify a split-virtual
      CISD[TQ] computation.  If FALSE, then if there are any electrons
      in RAS IV, then the number of holes in RAS I cannot exceed the given
      excitation level EX_LEVEL.  !expert -*/
    options.add_bool("MIXED4",true);

    /*- Do restrict strings with $e-$ in RAS IV?  Useful to reduce the number of
      strings required if MIXED4=true, as in a split-virutal CISD[TQ]
      computation.  If more than one electron is in RAS IV, then the
      holes in RAS I cannot exceed the number of particles in
      RAS III + RAS IV (i.e., EX_LEVEL), or else the string is discarded.
      !expert -*/
    options.add_bool("R4S",false);

    /*- Do string replacements on the fly in DETCI? Can
      save a gigantic amount of memory (especially for truncated CI's) but
      is somewhat flaky and hasn't been tested for a while.  It may work
      only works for certain classes of RAS calculations.  The current
      code is very slow with this option turned on. !expert -*/
    options.add_bool("REPL_OTF",false);

    /*- Do calculate the value of $<S^2>$ for each root? -*/
    options.add_bool("CALC_SSQ",false);

    /*- Do compute the MPn series out to
    kth order where k is determined by MAX_NUM_VECS?  For open-shell systems
    (REF=ROHF, WFN = ZAPTN), DETCI will compute the ZAPTn series.
    GUESS_VECTOR must be set to UNIT, HD_OTF must be set to TRUE, and
    HD_AVE must be set to orb_ener; these should happen by default for
    MPN=TRUE. -*/
    options.add_bool("MPN",false);

    /*- If 0, save the MPn energy; if 1, save the MP(2n-1) energy (if
    available from WIGNER=true); if 2, save the MP(2n-2) energy (if
    available from WIGNER=true). !expert -*/
    options.add_int("SAVE_MPN2",0);

    /*- Do employ an orthonormal vector space rather than
      storing the kth order wavefunction? !expert -*/
    options.add_bool("MPN_SCHMIDT",false);

    /*- Do use Wigner formulas in the Empn series? !expert -*/
    options.add_bool("WIGNER",true);

    /*- $z$ in $H = H@@0 + z H@@1$ !expert -*/
    options.add_double("PERTURB_MAGNITUDE",1.0);

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
    the value of MAX_NUM_VECS, the Davidson subspace will be
    collapsed to COLLAPSE_SIZE vectors for each root.  This is very
    helpful for saving disk space.  Defaults to MAXITER * NUM_ROOTS
    + NUM_INIT_VECS. -*/
    options.add_int("MAX_NUM_VECS", 0);

    /*- Gives the number of vectors to retain when the Davidson subspace is
    collapsed (see MAX_NUM_VECS below).  If greater than one, the
    collapsed subspace retains the best estimate of the CI vector for
    the previous n iterations.   Defaults to 1. -*/
    options.add_int("COLLAPSE_SIZE", 1);

    /*- Do use least-squares extrapolation in iterative solution of CI
    vector? -*/
    options.add_bool("LSE",false);

    /*- Number of iterations between least-squares extrapolations -*/
    options.add_int("LSE_COLLAPSE", 3);

    /*- Minimum converged energy for least-squares
    extrapolation to be performed -*/
    options.add_double("LSE_TOLERANCE", 3);

    /*- Do result a DETCI iteration that
    terminated prematurely? It assumes that the CI and sigma vectors are on
    disk; the number of vectors specified by RESTART_VECS is collapsed
    down to one vector per root. -*/
    options.add_bool("RESTART",false);

    /*- Do use some routines based on the papers of Bendazzoli
    et al. to calculate sigma?  Seems to be slower and not worthwhile; may disappear
    eventually.  Works only for full CI and I don't remember if I could see
    how their clever scheme might be extended to RAS in general. !expert -*/
    options.add_bool("BENDAZZOLI", false);

    /*- Do coupled-cluster computation? -*/
    options.add_bool("CC", false);

    /*- Do compute one-particle density matrix if not otherwise required? -*/
    options.add_bool("OPDM", false);

    /*- Do compute two-particle density matrix if not otherwise required? -*/
    options.add_bool("TPDM", false);

    /*- Maximum number of iterations to diagonalize the Hamiltonian -*/
    options.add_int("MAXITER", 12);

    /*- Do print the one-particle density matrix for each root? -*/
    options.add_bool("OPDM_PRINT", false);

    /*- Do write the natural orbitals? -*/
    options.add_bool("NAT_ORBS_WRITE", false);

    /*- Do average the OPDM over several roots in
    order to obtain a state-average one-particle density matrix?  This
    density matrix can be diagonalized to obtain the CI natural orbitals. -*/
    options.add_bool("OPDM_AVE", false);

    /*- Sets the root number for which CI natural orbitals are written
    to PSIF_CHKPT.  The default value is 1 (lowest root). -*/
    options.add_int("NAT_ORBS_WRITE_ROOT", 1);

    /*- Do compute the kinetic energy contribution from the correlated part of
    the one-particle density matrix? !expert -*/
    options.add_bool("OPDM_KE", false);

    /*- Do print the two-particle density matrix? (Warning: large tensor) -*/
    options.add_bool("TPDM_PRINT", false);

    /*- The root to write out the two-particle density matrix for
    (the one-particle density matrices are written for all roots).
    Useful for a state-specific CASSCF or CI optimization on an
    excited state. -*/
    options.add_int("FOLLOW_ROOT", 1);

    /*- Do compute the transition density?  Note: only transition densities
    between roots of the same symmetry will be evaluated.  DETCI
    does not compute states of different irreps within the same
    computation; to do this, lower the symmetry of the computation.-*/
    options.add_bool("TRANSITION_DENSITY", false);

    /*- Do write the transition density? -*/
    options.add_bool("TDM_WRITE", false);

    /*- Do print the transition density? -*/
    options.add_bool("TDM_PRINT", false);

    /*- Do compute the dipole moment? -*/
    options.add_bool("DIPMOM", false);

    /*- Number of threads -*/
    options.add_int("NUM_THREADS", 1);

    /*- Do store converged vector(s) at the end of
    the run?  The vector(s) is(are) stored in a transparent format such that
    other programs can use it easily. The format is specified in
    src/lib/libqt/slaterdset.h. -*/
    options.add_bool("EXPORT_VECTOR", false);

    /*- Number of vectors to export -*/
    options.add_int("NUM_EXPORT", 1);

    /*- Do eliminate determinants not valid for spin-complete spin-flip CI's?
    [see J. S. Sears et al, J. Chem. Phys. 118, 9084-9094 (2003)] !expert -*/
    options.add_bool("SF_RESTRICT", false);

    /*- Do print the sigma overlap matrix?  Not generally useful.  !expert -*/
    options.add_bool("SIGMA_OVERLAP", false);

    /*- The value of the spin quantum number S is given by this option.
    The default is determined by the value of the multiplicity.  This is used 
    for two things: (1) determining the phase of the redundant half of the CI 
    vector when the Ms=0 component is used (i.e., Ms0 = TRUE), and (2) making 
    sure the guess vector has the desired value of $<S^2>$ (if CALC_SSQ is TRUE
    and ICORE=1). -*/
    options.add_double("S", 0.0);

    /*- An array of length EX_LEVEL specifying whether each excitation type
    (S,D,T, etc.) is allowed (1 is allowed, 0 is disallowed).  Used to
    specify non-standard CI spaces such as CIST.  !expert -*/
    options.add("EX_ALLOW", new ArrayType());

    /*- Do invoke the FILTER_GUESS options that are used to filter out some trial
    vectors which may not have the appropriate phase convention
    between two determinants?  This is useful to remove, e.g.,
    delta states when a sigma state is desired.  The user
    inputs two determinants (by giving the absolute alpha string
    number and beta string number for each), and also the
    desired phase between these two determinants for guesses
    which are to be kept.  FILTER_GUESS = TRUE turns on the filtering
    routine.  Requires additional keywords FILTER_GUESS_DET1,
    FILTER_GUESS_DET2, and FILTER_GUESS_SIGN. !expert -*/
    options.add_bool("FILTER_GUESS", false);

    /*- The required phase (1 or -1) between the two determinants specified
    by FILTER_GUESS_DET1 and FILTER_GUESS_DET2 !expert -*/
    options.add_int("FILTER_GUESS_SIGN", 1);

    /*- Array specifying the absolute alpha string number and beta string
    number for the first determinant in the filter procedure.
    (See FILTER_GUESS).  !expert -*/
    options.add("FILTER_GUESS_DET1", new ArrayType());

    /*- Array specifying the absolute alpha string number and beta string
    number for the second determinant in the filter procedure.
    (See FILTER_GUESS).  !expert -*/
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

    /*- Do export a CC vector to disk? -*/
    options.add_bool("CC_EXPORT", false);

    /*- Do import a CC vector from disk? -*/
    options.add_bool("CC_IMPORT", false);

    /*- Do fix amplitudes involving RAS I or RAS IV?  Useful in mixed
    MP2-CC methods. !expert -*/
    options.add_bool("CC_FIX_EXTERNAL", false);

    /*- Number of external indices before amplitude gets fixed by
    CC_FIX_EXTERNAL.  Experimental. !expert -*/
    options.add_int("CC_FIX_EXTERNAL_MIN", 1);

    /*- Do use variational energy expression in CC computation?
    Experimental.  !expert -*/
    options.add_bool("CC_VARIATIONAL", false);

    /*- Do ignore block if num holes in RAS I and II is > cc_ex_lvl and if
    any indices correspond to RAS I or IV (i.e., include only all-active
    higher excitations)? !expert -*/
    options.add_bool("CC_MIXED", true);

    /*- Do update T amplitudes with orbital eigenvalues? (Usually would
    do this).  Not doing this is experimental.  !expert -*/
    options.add_bool("CC_UPDATE_EPS", true);

    /*- Do use DIIS extrapolation to accelerate CC convergence? -*/
    options.add_bool("DIIS", true);

    /*- Iteration at which to start using DIIS -*/
    options.add_int("DIIS_START_ITER", 1);

    /*- How often to do a DIIS extrapolation. 1 means do DIIS every
    iteration, 2 is every other iteration, etc. -*/
    options.add_int("DIIS_FREQ", 1);

    /*- Minimum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("DIIS_MIN_VECS", 2);

    /*- Maximum number of error vectors stored for DIIS extrapolation -*/
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
    /*- MODULEDESCRIPTION Performs symmetry adapted perturbation theory (SAPT) analysis
        to quantitatively analyze noncovalent interactions. -*/
    /*- The level of theory for SAPT -*/
    options.add_str("SAPT_LEVEL","SAPT0","SAPT0 SAPT2 SAPT2+ SAPT2+3 MP2C");
    /*- Convergence criterion for energy (change) in the SAPT Ind20 term. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("E_CONVERGENCE",1e-10);
    /*- Convergence criterion for density in the SAPT Ind20 term. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("D_CONVERGENCE",1e-8);
    /*- Don't solve the CPHF equations? -*/
    options.add_bool("NO_RESPONSE",false);
    /*- Do use asynchronous I/O in the CPHF solver? -*/
    options.add_bool("AIO_CPHF",false);
    /*- Do use asynchronous I/O in the DF integral formation? -*/
    options.add_bool("AIO_DF_INTS",false);
    /*- Maxmum number of CPHF iterations -*/
    options.add_int("MAXITER",50);
    /*- Do compute third-order corrections? -*/
    options.add_bool("DO_THIRD_ORDER",false);
    /*- Do compute natural orbitals? -*/
    options.add_bool("NAT_ORBS",false);
    /*- Do use natural orbitals for T2's? -*/
    options.add_bool("NAT_ORBS_T2",false);
    /*- Minimum occupation below which natural orbitals are neglected. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("OCC_TOLERANCE",1.0E-6);
    /*- Minimum absolute value below which integrals are neglected. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("INTS_TOLERANCE",1.0E-12);
    /*- Memory safety -*/
    options.add_double("SAPT_MEM_SAFETY",0.9);
    /*- Do force SAPT2 and higher to die if it thinks there isn't enough memory? -*/
    options.add_bool("SAPT_MEM_CHECK",true);
    /*- Auxiliary basis set for SAPT density fitting computations. Defaults to BASIS-RI. -*/
    options.add_str("DF_BASIS_SAPT", "");
    /*- Auxiliary basis set for SAPT Elst10 and Exch10 density fitting computations. Defaults to BASIS-RI. -*/
    options.add_str("DF_BASIS_ELST", "");
    /*- Maximum denominator error allowed (Max error norm in Delta tensor) -*/
    options.add_double("DENOMINATOR_DELTA", 1.0E-6);
    /*- Denominator algorithm for PT methods -*/
    options.add_str("DENOMINATOR_ALGORITHM", "LAPLACE", "LAPLACE CHOLESKY");
    /*- Number of Omega points for Casimir-Polder -*/
    options.add_int("OMEGA_POINTS",8);
    /*- Omega (atomic wavenumbers) to center Casimir-Polder on -*/
    options.add_double("OMEGA_CENTER", 0.4);
    /*- The scale factor used for opposite-spin pairs in SCS computations -*/
    options.add_double("SCALE_OS", 6.0/5.0);
    /*- The scale factor used for same-spin pairs in SCS computations-*/
    options.add_double("SCALE_SS", 1.0/3.0);
  }
  if(name == "DCFT"|| options.read_globals()) {
//      ip_cwk_add(":DCFT");
      /*- How to cache quantities within the DPD library -*/
      options.add_int("CACHELEVEL", 2);
      /*- The shift applied to the denominator -*/
      options.add_double("REGULARIZER", 0.0);
      /*- Maximum number of lambda iterations per macro-iteration -*/
      options.add_int("LAMBDA_MAXITER", 50);
      /*- Maximum number of SCF iterations per cycle -*/
      options.add_int("SCF_MAXITER", 50);
      /*- Maximum number of iterations -*/
      options.add_int("MAXITER", 40);
      /*- Do compute the full two particle density matrix at the end of the computation, for properties? -*/
      options.add_bool("COMPUTE_TPDM", 0);
      /*- Convergence criterion for the SCF density (RMS error). 
      See the note at the beginning of Section \ref{keywords}. -*/
      options.add_double("SCF_D_CONVERGENCE", 1e-8);
      /*- Convergence criterion for residuals (RMS error) in density cummulant equations.
      See the note at the beginning of Section \ref{keywords}. -*/
      options.add_double("R_CONVERGENCE", 1e-10);
      /*- Do relax the orbitals? -*/
      options.add_bool("MO_RELAX", true);
      /*- The amount (percentage) of damping to apply to the initial SCF procedures
      0 will result in a full update, 100 will completely stall the update. A
      value around 20 (which corresponds to 20\% of the previous iteration's
      density being mixed into the current iteration)
      can help in cases where oscillatory convergence is observed. -*/
      options.add_double("DAMPING_PERCENTAGE",0.0);
      /*- Don't include the tau terms? -*/
      options.add_bool("IGNORE_TAU", false);
      /*- Do compute the DCFT energy with the $\tau^{2}$ correction to $\tau$? -*/
      options.add_bool("TAU_SQUARED", false);
      /*- Minimum absolute value below which integrals are neglected. 
      See the note at the beginning of Section \ref{keywords}. -*/
      options.add_double("INTS_TOLERANCE", 1e-14);
      /*- DIIS starts when the  RMS lambda and SCF errors are less than $10^{-diis_start}$ -*/
      options.add_double("DIIS_START_CONVERGENCE", 1e-3);
      /*- Maximum number of error vectors stored for DIIS extrapolation -*/
      options.add_int("DIIS_MAX_VECS", 6);
      /*- Minimum number of error vectors stored for DIIS extrapolation -*/
      options.add_int("DIIS_MIN_VECS", 3);
      /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
      options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
      /*- The algorithm to use for lambda and orbital updates -*/
      options.add_str("ALGORITHM", "SIMULTANEOUS", "TWOSTEP SIMULTANEOUS");
      /*- Do force the occupation to be that of the SCF starting point? -*/
      options.add_bool("LOCK_OCCUPATION", true);
  }
  if (name == "MINTS"|| options.read_globals()) {
      /*- Primary basis set -*/
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
    /*- The boxing scheme for DFT -*/
    options.add_str("DFT_BOXING_SCHEME","NAIVE","NAIVE OCTREE");
    /*- DFT basis cutoff. See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("DFT_BASIS_TOLERANCE", 0.0);
    /*- The DFT combined functional name (for now) -*/
    options.add_str("DFT_FUNCTIONAL", "");
    /*- The DFT Range-separation parameter (only used if changed by the user) -*/
    options.add_double("DFT_OMEGA", 0.0);
  }
  if (name == "SCF"|| options.read_globals()) {

    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "SCF", "SCF");

    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF", "RHF ROHF UHF CUHF RKS UKS");

    /*- The dimension sizes of the processor grid -*/
    options.add("PROCESS_GRID", new ArrayType());

    /*- The tile size for the distributed matrices -*/
    options.add_int("TILE_SZ", 512);

    /*- The dimension sizes of the distributed matrix -*/
    options.add("DISTRIBUTED_MATRIX", new ArrayType());

    /*- Are going to do SAPT? If so, what part?  -*/
    options.add_str("SAPT","FALSE","FALSE 2-DIMER 2-MONOMER_A 2-MONOMER_B 3-TRIMER 3-DIMER_AB 3-DIMER_BC 3-DIMER_AC 3-MONOMER_A 3-MONOMER_B 3-MONOMER_C");

    /*- Primary basis set -*/
    options.add_str("BASIS", "");
    /*- Auxiliary basis set for SCF density fitting computations. Defaults to BASIS-JKFIT. -*/
    options.add_str("DF_BASIS_SCF", "");

    /*- Atomic Charge cutoff (for primary domain) -*/
    options.add_double("CHARGE_CUTOFF",0.05);
    /*- Extended domain radius, Angstrom -*/
    options.add_double("R_EXT",3.0);
    /*- Iterations per full Pipek-Mizey Localization -*/
    options.add_int("LOCAL_ITER",1);
    /*- What algorithm to use for the SCF computation -*/
    options.add_str("SCF_TYPE","PK","PK OUT_OF_CORE DIRECT DF PSEUDOSPECTRAL POISSON L_DF CD 1C_CD");
    /*- Do run in parallel? -*/
    options.add_bool("PARALLEL", false);

    /*- Primary basis set -*/
    options.add_str("BASIS","");

    /*- The type of guess orbitals -*/
    options.add_str("GUESS", "CORE", "CORE GWH SAD READ");
    /*- Maximum number of iterations -*/
    options.add_int("MAXITER", 100);

    /*- The amount (percentage) of damping to apply to the early density updates.
        0 will result in a full update, 100 will completely stall the update.  A
        value around 20 (which corresponds to 20\% of the previous iteration's
        density being mixed into the current density)
        could help to solve problems with oscillatory convergence. -*/
    options.add_double("DAMPING_PERCENTAGE", 100.0);
    /*- The density convergence threshold after which damping is no longer performed, if it is enabled.
        It is recommended to leave damping on until convergence, which is the default.
        See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("DAMPING_CONVERGENCE", 1.0E-18);

    /*- Do perturb the Hamiltonian? -*/
    options.add_bool("PERTURB_H", false);
    /*- Size of the perturbation -*/
    options.add_double("PERTURB_MAGNITUDE", 0.0);
    /*- The operator used to perturb the Hamiltonian, if requested -*/
    options.add_str("PERTURB_WITH", "DIPOLE_X", "DIPOLE_X DIPOLE_Y DIPOLE_Z");
    /*- Do look for an ExternalPotential object in Python? -*/
    options.add_bool("EXTERN", false);

    /*- The storage scheme for the three index tensors in density fitting -*/
    options.add_str("DF_SCF_STORAGE", "DEFAULT", "DEFAULT CORE DISK");
    /*- Do save restart information for RI SCF? -*/
    options.add_bool("DF_SCF_SAVE",false);
    /*- Do try to restart? -*/
    options.add_bool("DF_SCF_RESTART",false);
    /*- Do find the raw fitting metric condition number? -*/
    options.add_bool("FIND_RAW_J_COND",false);
    /*- Max J basis condition number to be allowed -*/
    options.add_double("DF_MAX_COND",1E8);
    /*- SCF Fitting Type -*/
    options.add_str("DF_FITTING_TYPE", "FINISHED", "FINISHED RAW CHOLESKY");
    /*- Number of threads for integrals (may be turned down if memory is an issue). 0 is blank -*/
    options.add_int("DF_INTS_NUM_THREADS",0);
    /*- IO caching for CP corrections, etc -*/
    options.add_str("DF_INTS_IO", "NONE", "NONE SAVE LOAD");

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
    options.add_int("DIIS_START", 1);
    /*- Minimum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("DIIS_MIN_VECS", 2);
    /*- Maximum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("DIIS_MAX_VECS", 10);
    /*- Do use DIIS extrapolation to accelerate convergence? -*/
    options.add_bool("DIIS", true);
    /*- Convergence criterion for SCF energy. See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("E_CONVERGENCE", 1e-8);
    /*- Convergence criterion for SCF density. See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("D_CONVERGENCE", 1e-8);
    /*- Minimum absolute value below which TEI are neglected. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("INTS_TOLERANCE", 0.0);
    /*- Maximum number of rows to read/write in each DF-SCF operation -*/
    options.add_int("ROWS_PER_READ", 0);

    /*- The amount of SAD information to print to the output -*/
    options.add_int("SAD_PRINT", 0);
    /*- Convergence criterion for SCF energy in SAD Guess. See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("SAD_E_CONVERGENCE", 1E-5);
    /*- Convergence criterion for SCF density in SAD Guess. See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("SAD_D_CONVERGENCE", 1E-5);
    /*- Maximum number of SAD guess iterations -*/
    options.add_int("SAD_MAXITER", 50);
    /*- SAD Guess F-mix Iteration Start -*/
    options.add_int("SAD_F_MIX_START", 50);
    /*- SAD Guess Cholesky Cutoff (for eliminating redundancies).
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("SAD_CHOL_TOLERANCE", 1E-7);
  }
  if (name == "MP2"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "MP2", "MP2");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF", "RHF UHF ROHF");
    /*- Type of job being performed !expert -*/
    options.add_str("JOBTYPE", "SP");
    /*- Do compute the one particle density matrix, for properties? -*/
    options.add_bool("COMPUTE_OPDM", false);
    /*- Do add relaxation terms to the one particle density matrix, for properties? -*/
    options.add_bool("RELAX_OPDM", false);
    /*- The amount of cacheing of data to perform -*/
    options.add_int("CACHELEVEL", 2);
    /*- The criterion used to retain/release cached data -*/
    options.add_str("CACHETYPE", "LRU", "LRU LOW");
    /*- Do perform a spin component scaled MP2 computation? -*/
    options.add_bool("SCS", false);
    /*- Do perform a spin component scaled (N) MP2 computation? -*/
    options.add_bool("SCS_N", false);
    /*- The scale factor used for opposite-spin pairs in SCS computations -*/
    options.add_double("SCALE_OS", 6.0/5.0);
    /*- The scale factor used for same-spin pairs in SCS computations-*/
    options.add_double("SCALE_SS", 1.0/3.0);
  }
  // Options of this module not standardized since it's bound for deletion
  if(name == "TRANSQT2"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
    /*- Do ? -*/
    options.add_bool("PRINT_TEI", false);
    /*- Minimum absolute value below which integrals are neglected. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("INTS_TOLERANCE", 1e-14);
    /*- -*/
    options.add_int("CACHELEVEL", 2);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- Boolean to delete the SO-basis two-electron integral file after the transformation -*/
    options.add_bool("DELETE_TEI", true);
    /*- Convert ROHF MOs to semicanonical MOs -*/
    options.add_bool("SEMICANONICAL", true);
  }
  // Options of this module not standardized since it's bound for deletion
  if(name == "TRANSQT"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "CCSD");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
    /*- -*/
    options.add_int("PRINT_LVL", 1);
    /*- -*/
    options.add_str("MODE", "TO_MO", "TO_MO TO_AO");
    /*- Do ? -*/
    options.add_bool("PSIMRCC", false);
    /*- -*/
    options.add_str("MP2R12A", "MP2R12AERI", "MP2R12AERI MP2R12AR12 MP2R12AR12T1");
    /*- Minimum absolute value below which integrals are neglected. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("INTS_TOLERANCE", 1e-14);
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
    /*- Do ? -*/
    options.add_bool("KEEP_PRESORT", false);
    /*- -*/
    options.add_int("J_FILE", 91);
    /*- Do keep half-transformed integrals? -*/
    options.add_bool("KEEP_J", false);
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
    /*- Don't ? -*/
    options.add_bool("DELETE_AO", true);
    /*- Don't ? -*/
    options.add_bool("DELETE_TPDM", true);
    /*- Do ? -*/
    options.add_bool("PRINT_TE_INTEGRALS", false);
    /*- Do ? -*/
    options.add_bool("PRINT_OE_INTEGRALS", false);
    /*- Do ? -*/
    options.add_bool("PRINT_SORTED_OE_INTS", false);
    /*- Do ? -*/
    options.add_bool("PRINT_SORTED_TE_INTS", false);
    /*- Do ? -*/
    options.add_bool("PRINT_MOS", false);

    /*- Do ? -*/
    options.add_bool("LAGRAN_DOUBLE", false);
    /*- Do ? -*/
    options.add_bool("LAGRAN_HALVE", false);
    /*- Do ? -*/
    options.add_bool("DO_ALL_TEI", false);
    /*- Do ? -*/
    options.add_bool("TPDM_ADD_REF", false);
    /*- Don't ? -*/
    options.add_bool("DELETE_RESTR_DOCC", true);
    /*- Do ? -*/
    options.add_bool("PRINT_REORDER", false);
    /*- Do ? -*/
    options.add_bool("PITZER", false);
    /*- Do ? -*/
    options.add_bool("REORDER", false);
    /*- Do ? -*/
    options.add_bool("CHECK_C_ORTHONORM", false);
    /*- Do ? -*/
    options.add_bool("QRHF", false);
    /*- Do ? -*/
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
  if(name == "CCSORT"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF");
    /*- -*/
    options.add_str("EOM_REFERENCE","RHF");
    /*- -*/
    options.add_str("PROPERTY", "POLARIZABILITY");
    /*- Do ? -*/
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
    /*- Do ? -*/
    options.add_bool("LOCAL_DOMAIN_POLAR", false);
    /*- Do ? -*/
    options.add_bool("LOCAL_DOMAIN_MAG", false);
    /*- Do ? -*/
    options.add_bool("LOCAL_DOMAIN_SEP", false);
    /*- Do ? -*/
    options.add_bool("LOCAL_FILTER_SINGLES", false);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- Do ? -*/
    options.add_bool("KEEP_TEIFILE", false);
    /*- Do ? -*/
    options.add_bool("KEEP_OEIFILE", false);
    /*- Minimum absolute value below which integrals are neglected. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("INTS_TOLERANCE", 1e-14);
    /*- -*/
    options.add_int("CACHELEVEL", 2);
    /*- Do ? -*/
    options.add_bool("LOCAL", false);
    /*- -*/
    options.add("OMEGA", new ArrayType());
    /*- -*/
    options.add_str("OMEGA_UNITS", "AU", "AU HZ EV NM");
    /*- Convert ROHF MOs to semicanonical MOs -*/
    options.add_bool("SEMICANONICAL", true);
  }
  if(name == "CCTRIPLES"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "SCF");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
    /*- Number of threads -*/
    options.add_int("NUM_THREADS",1);
    /*- Convert ROHF MOs to semicanonical MOs -*/
    options.add_bool("SEMICANONICAL", true);
  }
  if(name == "CCDENSITY"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "SCF");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
    /*- Minimum absolute value below which integrals are neglected. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("INTS_TOLERANCE",1e-14);
    /*- The amount of cacheing of data to perform -*/
    options.add_int("CACHELEVEL",2);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- Do compute the approximate excitation level? See Stanton and Bartlett, JCP, 98, 1993, 7034. !expert -*/
    options.add_bool("AEL",false);
    /*- The type of gauge to use for properties -*/
    options.add_str("GAUGE","LENGTH");
    /*- Do relax the one-particle density matrix? -*/
    options.add_bool("RELAX_OPDM",false);
    /*- Do require $\bar{H}$ and $R$ to be connected? !expert -*/
    options.add_bool("CONNECT_XI",false);
    /*- The number of electronic states to computed, per irreducible representation -*/
    options.add("STATES_PER_IRREP", new ArrayType());
    /*- Do compute all relaxed excited states? -*/
    options.add_bool("PROP_ALL",false);
    /*- The symmetry of states -*/
    options.add_int("PROP_SYM", 1);
    /*- -*/
    options.add_int("PROP_ROOT", 1);
  }
  if(name == "CCLAMBDA"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN","SCF");
    /*- Convergence criterion for wavefunction (change) in CC lambda-amplitude equations.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("R_CONVERGENCE",1e-7);
    /*- Do ? -*/
    options.add_bool("RESTART",false);
    /*- -*/
    options.add_int("CACHELEVEL",2);
    /*- Do ? -*/
    options.add_bool("SEKINO",false);
    /*- Do use DIIS extrapolation to accelerate convergence? -*/
    options.add_bool("DIIS",true);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- -*/
    options.add_str("ABCD","NEW");
    /*- -*/
    options.add_int("NUM_AMPS_PRINT",10);
    /*- Type of job being performed !expert -*/
    options.add_str("JOBTYPE","");
    /*- Do ? -*/
    options.add_bool("LOCAL",false);
    /*- -*/
    options.add_str("LOCAL_WEAKP","NONE");
    /*- -*/
    options.add_double("LOCAL_CUTOFF",0.02);
    /*- -*/
    options.add_str("LOCAL_METHOD","WERNER");
    /*- Do ? -*/
    options.add_bool("LOCAL_FILTER_SINGLES",true);
    /*- -*/
    options.add_double("LOCAL_CPHF_CUTOFF",0.10);
    /*- -*/
    options.add_str("LOCAL_PAIRDEF","");
    /*- -*/
    options.add("STATES_PER_IRREP", new ArrayType());
    /*- Do ? -*/
    options.add_bool("PROP_ALL",false);
    /*- -*/
    options.add_int("PROP_SYM",1);
    /*- -*/
    options.add_int("PROP_ROOT",1);
    /*- Maximum number of iterations -*/
    options.add_int("MAXITER",50);
  }
  if(name == "CLAG"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN","NONE");
    /*- Do write the OEI, TEI, OPDM, TPDM, and Lagrangian files in canonical form, Pitzer order? -*/
    options.add_bool("CAS_FILES_WRITE",0);
    /*- Root to get OPDM -*/
    options.add_int("FOLLOW_ROOT",1);
  }
  if(name == "STABLE"|| options.read_globals()) {
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
    /*- -*/
    options.add_int("CACHELEVEL",2);
    /*- Do ? -*/
    options.add_bool("FOLLOW",false);
    /*- -*/
    options.add_int("NUM_EVECS_PRINT",0);
    /*- -*/
    options.add_int("ROTATION_METHOD",0);
    /*- -*/
    options.add_double("SCALE",0.5);
  }
  if(name == "ADC" || options.read_globals()) {
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "");
    /*- How to cache quantities within the DPD library -*/
    options.add_int("CACHELEVEL", 2);
    /*- The amount of memory available (in Mb) -*/
    options.add_int("MEMORY", 1000);
    /*- The convergence criterion for pole searching step. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("NEWTON_CONVERGENCE", 1e-7);
    /*- The maximum numbers of the pole searching iteration  -*/
    options.add_int("POLE_MAX", 20);
    /*- Maximu iteration number in simultaneous expansion method -*/
    options.add_int("SEM_MAX", 30);
    /*- The cutoff norm of residual vector in SEM step.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("NORM_TOLERANCE", 1e-6);
    /*- The poles per irrep vector -*/
    options.add("STATES_PER_IRREP", new ArrayType());
    /*- Do use the partial renormalization scheme for the ground state wavefunction? -*/
    options.add_bool("PR", false);
  }
  if(name == "CCHBAR"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "SCF");
    /*- -*/
    options.add_str("EOM_REFERENCE","RHF");
    /*- Do compute the Tamplitude equation matrix elements? -*/
    options.add_bool("T_AMPS",false);
    /*- -*/
    options.add_int("CACHELEVEL",2);
    /*- Do use the minimal-disk algorithm for Wabei? It's VERY slow! -*/
    options.add_bool("WABEI_LOWDISK", false);
  }
  if(name == "CCEOM"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "EOM_CCSD", "EOM_CCSD EOM_CC2 EOM_CC3");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF", "RHF ROHF UHF");
    /*- -*/
    options.add_str("EOM_REFERENCE","RHF", "RHF ROHF UHF");
    /*- Do ? -*/
    options.add_bool("FULL_MATRIX",false);
    /*- -*/
    options.add_int("CACHELEVEL",2);
    /*- -*/
    options.add_str("CACHETYPE", "LRU", "LOW LRU");
    /*- Number of threads -*/
    options.add_int("NUM_THREADS", 1);
    /*- -*/
    options.add_str("ABCD", "NEW", "NEW OLD");
    /*- Do ? -*/
    options.add_bool("T3_WS_INCORE", false);
    /*- Do ? -*/
    options.add_bool("LOCAL", false);
    /*- -*/
    options.add_double("LOCAL_CUTOFF", 0.02);
    /*- -*/
    options.add_str("LOCAL_METHOD", "WERNER", "WERNER AOBASIS");
    /*- -*/
    options.add_str("LOCAL_WEAKP", "NONE", "NONE MP2 NEGLECT");
    /*- -*/
    options.add_str("LOCAL_PRECONDITIONER", "HBAR", "HBAR FOCK");
    /*- -*/
    options.add_int("LOCAL_GHOST", -1);
    /*- Do ? -*/
    options.add_bool("LOCAL_DO_SINGLES", true);
    /*- Do ? -*/
    options.add_bool("LOCAL_FILTER_SINGLES", true);
    /*- Do ? -*/
    options.add_bool("NEW_TRIPLES", true);
    /*- Number of excited states per irreducible representation for EOM-CC
    and CC-LR calculations. Irreps denote the final state symmetry, not the
    symmetry of the transtion. -*/
    options.add("STATES_PER_IRREP", new ArrayType());
    /*- Maximum number of iterations -*/
    options.add_int("MAXITER", 80);
    /*- -*/
    options.add_int("PROP_SYM", 1);
    /*- -*/
    options.add_int("PROP_ROOT", 0);
    /*- Do ? -*/
    options.add_bool("CC3_FOLLOW_ROOT", false);
    /*- Do ? -*/
    options.add_bool("RHF_TRIPLETS", false);
    /*- The depth into the occupied and valence spaces from which one-electron 
    excitations are seeded into the Davidson guess to the CIS (the default of 2 
    includes all single excitations between HOMO-1, HOMO, LUMO, and LUMO+1). This 
    CIS is in turn the Davidson guess to the EOM-CC. Expand to capture more exotic 
    excited states in the EOM-CC calculation !expert -*/
    options.add_int("EXCITATION_RANGE", 2);
    /*- Do print information on the iterative solution to the single-excitation
    EOM-CC problem used as a guess to full EOM-CC? -*/
    options.add_bool("SINGLES_PRINT", false);
    /*- -*/
    options.add_int("SS_VECS_PER_ROOT", 5);
    /*- -*/
    options.add_int("VECS_PER_ROOT", 12);
    /*- -*/
    options.add_int("VECS_CC3", 10);
    /*- Do ? -*/
    options.add_bool("COLLAPSE_WITH_LAST", true);
    /*- See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("COMPLEX_TOLERANCE", 1E-12);
    /*- Convergence criterion for norm of the residual vector in the Davidson algorithm for CC-EOM.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("R_CONVERGENCE", 1E-6);
    /*- Convergence criterion for norm of the residual vector in the Davidson algorithm for the CIS guess to CC-EOM.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("SS_R_CONVERGENCE", 1E-6);
    /*- Convergence criterion for excitation energy (change) in the Davidson algorithm for CC-EOM.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("E_CONVERGENCE", 1E-8);
    /*- Convergence criterion for excitation energy (change) in the Davidson algorithm for the CIS guess to CC-EOM.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("SS_E_CONVERGENCE", 1E-6);
    /*- Number of important CC amplitudes to print -*/
    options.add_int("NUM_AMPS_PRINT", 5);
    /*- Minimum absolute value above which a guess vector to a root is added to the Davidson algorithm.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("SCHMIDT_ADD_RESIDUAL_TOLERANCE", 1E-3);
    /*- Do ? -*/
    options.add_bool("SKIP_DIAGSS", false);
    /*- Do ? -*/
    options.add_bool("RESTART_EOM_CC3", false);
    /*- -*/
    options.add_str("EOM_GUESS", "SINGLES", "SINGLES DISK INPUT");
    /*- Convert ROHF MOs to semicanonical MOs -*/
    options.add_bool("SEMICANONICAL", true);
  }
  if(name == "CCRESPONSE"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "SCF");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
    /*- Cacheing level for libdpd -*/
    options.add_int("CACHELEVEL",2);
    /*- Gauge for optical rotation -*/
    options.add_str("GAUGE","LENGTH");
    /*- Maximum number of iterations to converge perturbed amplitude equations -*/
    options.add_int("MAXITER",50);
    /*- Convergence criterion for wavefunction (change) in perturbed CC equations.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("R_CONVERGENCE",1e-7);
    /*- Do use DIIS extrapolation to accelerate convergence? -*/
    options.add_bool("DIIS",1);
    /*- -*/
    options.add_str("PROPERTY","POLARIZABILITY");
    /*- -*/
    options.add_str("ABCD","NEW");
    /*- Do restart from on-disk amplitudes? -*/
    options.add_bool("RESTART",1);
    /*- Do simulate local correlation? -*/
    options.add_bool("LOCAL",0);
    /*- -*/
    options.add_double("LOCAL_CUTOFF",0.01);
    /*- -*/
    options.add_str("LOCAL_METHOD","WERNER");
    /*- -*/
    options.add_str("LOCAL_WEAKP","NONE");
    /*- Do ? -*/
    options.add_bool("LOCAL_FILTER_SINGLES", false);
    /*- -*/
    options.add_double("LOCAL_CPHF_CUTOFF",0.10);
    /*- -*/
    options.add_str("LOCAL_PAIRDEF","NONE");
    /*- Do ? -*/
    options.add_bool("ANALYZE",0);
    /*- -*/
    options.add_int("NUM_AMPS_PRINT",5);
    /*- Do Sekino-Bartlett size-extensive model-III? -*/
    options.add_bool("SEKINO",0);
    /*- Do Bartlett size-extensive linear model? -*/
    options.add_bool("LINEAR",0);
    /*- Energy of applied field for dynamic polarizabilities -*/
    options.add("OMEGA",new ArrayType());
    /*- -*/
    options.add("MU_IRREPS",new ArrayType());
  }
  if(name == "RESPONSE"|| options.read_globals()){
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF");
    /*- -*/
    options.add("OMEGA", new ArrayType());
    /*- -*/
    options.add_str("PROPERTY","POLARIZABILITY");
  }
  if(name == "MCSCF"|| options.read_globals()) {
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF","RHF ROHF UHF TWOCON MCSCF GENERAL");
    /*- Level shift to aid convergence -*/
    options.add_int("LEVELSHIFT",0);
    /*- Convergence criterion for energy. See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("E_CONVERGENCE", 1e-12);
    /*- Convergence criterion for density. See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("D_CONVERGENCE", 1e-12);
    /*- Maximum number of iterations -*/
    options.add_int("MAXITER",100);
    /*- Maximum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("DIIS_MAX_VECS",7);
    /*- Which solution of the SCF equations to find, where 1 is the SCF ground state-*/
    options.add_int("FOLLOW_ROOT",1);
    /*- Iteration at which to begin using the averaged Fock matrix-*/
    options.add_int("START_FAVG",5);
    /*- -*/
    options.add_int("TURN_ON_ACTV",0);
    /*- For orbital rotations after convergence, the angle (in degrees) by which to rotate. !expert -*/
    options.add_int("ROTATE_MO_ANGLE",0);
    /*- For orbital rotations after convergence, irrep (1-based, Cotton order) of the orbitals to rotate. !expert -*/
    options.add_int("ROTATE_MO_IRREP",1);
    /*- For orbital rotations after convergence, number of the first orbital (1-based) to rotate. !expert -*/
    options.add_int("ROTATE_MO_P",1);
    /*- For orbital rotations after convergence, number of the second orbital (1-based) to rotate. !expert -*/
    options.add_int("ROTATE_MO_Q",2);
    /*- Do use DIIS extrapolation to accelerate convergence of the CI coefficients? -*/
    options.add_bool("CI_DIIS",false);
    /*- Do use DIIS extrapolation to accelerate convergence of the SCF energy (MO coefficients only)? -*/
    options.add_bool("DIIS",true);
    /*- Do read in from file the MOs from a previous computation? -*/
    options.add_bool("MO_READ",true);
    /*- Do use the average Fock matrix during the SCF optimization? -*/
    options.add_bool("USE_FAVG",false);
    /*- Do canonicalize the active orbitals such that the average Fock matrix is diagonal? -*/
    options.add_bool("CANONICALIZE_ACTIVE_FAVG",false);
    /*- Do canonicalize the inactive (DOCC and Virtual) orbitals such that the average Fock matrix is diagonal? -*/
    options.add_bool("CANONICALIZE_INACTIVE_FAVG",false);
    /*- Do ? -*/
    options.add_bool("INTERNAL_ROTATIONS",true);
    /*- Do attempt to force a two configruation solution by starting with CI coefficents of $\pm \sqrt{\frac{1}{2}}$ -*/
    options.add_bool("FORCE_TWOCON",false);
    /*- The number of singly occupied orbitals, per irrep -*/
    options.add("SOCC", new ArrayType());
    /*- The number of doubly occupied orbitals, per irrep -*/
    options.add("DOCC", new ArrayType());
//    /*- The number of active orbitals, per irrep (alternative name for ACTIVE) -*/
//    options.add("ACTV", new ArrayType());
    /*- The symmetry of the SCF wavefunction.-*/
    options.add_str("WFN_SYM","1","A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");
  }
  if(name == "CCENERGY"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "NONE", "CCSD CCSD_T EOM_CCSD LEOM_CCSD BCCD BCCD_T CC2 CC3 EOM_CC2 EOM_CC3 CCSD_MVD");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF");
    /*- Do ? -*/
    options.add_bool("NEW_TRIPLES", 1);
    /*- Do ? -*/
    options.add_bool("ANALYZE", 0);
    /*- Maximum number of iterations to solve the CC equations -*/
    options.add_int("MAXITER", 50);
    /*- Convergence criterion for wavefunction (change) in CC amplitude equations.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("R_CONVERGENCE", 1e-7);
    /*- Do restart the coupled-cluster iterations from old $t@@1$ and $t@@2$ amplitudes? -*/
    options.add_bool("RESTART",1);
    /*- Do restart the coupled-cluster iterations even if MO phases are screwed up? !expert -*/
    options.add_bool("FORCE_RESTART", 0);
//#warning CCEnergy ao_basis keyword type was changed.
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- Cacheing level for libdpd governing the storage of amplitudes, 
    integrals, and intermediates in the CC procedure. A value of 0 retains 
    no quantities in cache, while a level of 6 attempts to store all 
    quantities in cache.  For particularly large calculations, a value of 
    0 may help with certain types of memory problems.  The default is 2, 
    which means that all four-index quantites with up to two virtual-orbital 
    indices (e.g., <ij|ab> integrals) may be held in the cache. -*/
    options.add_int("CACHELEVEL", 2);
    /*- Selects the priority type for maintaining the automatic memory 
    cache used by the libdpd codes. A value of LOW selects a "low priority" 
    scheme in which the deletion of items from the cache is based on 
    pre-programmed priorities. A value of LRU selects a "least recently used" 
    scheme in which the oldest item in the cache will be the first one deleted. -*/
    options.add_str("CACHETYPE", "LOW", "LOW LRU");
    /*- Number of threads -*/
    options.add_int("NUM_THREADS",1);
    /*- Do use DIIS extrapolation to accelerate convergence? -*/
    options.add_bool("DIIS", true);
    /*- Do ? -*/
    options.add_bool("T2_COUPLED", false);
    /*- -*/
    options.add_str("PROPERTY", "POLARIZABILITY", "POLARIZABILITY ROTATION MAGNETIZABILITY ROA ALL");
    /*- -*/
    options.add_str("ABCD", "NEW", "NEW OLD");
    /*- Do simulate the effects of local correlation techniques? -*/
    options.add_bool("LOCAL", 0);
    /*- Value (always between one and zero) for the Broughton-Pulay completeness 
    check used to contruct orbital domains for local-CC calculations. See 
    J. Broughton and P. Pulay, J. Comp. Chem. 14, 736-740 (1993) and C. Hampel 
    and H.-J. Werner, J. Chem. Phys. 104, 6286-6297 (1996). -*/
    options.add_double("LOCAL_CUTOFF", 0.02);
    /*- Type of local-CCSD scheme to be simulated. WERNER selects the method 
    developed by H.-J. Werner and co-workers, and AOBASIS selects the method 
    developed by G.E. Scuseria and co-workers (currently inoperative). -*/
    options.add_str("LOCAL_METHOD", "WERNER", "WERNER AOBASIS");
    /*- Desired treatment of "weak pairs" in the local-CCSD method. A value of 
    NEGLECT ignores weak pairs entirely. A value of NONE treats weak pairs in 
    the same manner as strong pairs. A value of MP2 uses second-order perturbation 
    theory to correct the local-CCSD energy computed with weak pairs ignored. -*/
    options.add_str("LOCAL_WEAKP", "NONE", "NONE NEGLECT MP2");
    //options.add_int("LOCAL_FILTER_SINGLES", 1);
    options.add_double("LOCAL_CPHF_CUTOFF", 0.10);
    /*- -*/
    options.add_str("LOCAL_PAIRDEF", "BP", "BP RESPONSE");
    /*- Number of important $t@@1$ and $t@@2$ amplitudes to print -*/
    options.add_int("NUM_AMPS_PRINT", 10);
    /*- Convergence criterion for Breuckner orbitals. See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("BRUECKNER_ORBS_R_CONVERGENCE", 1e-5);
    /*- Do print the MP2 amplitudes which are the starting guesses for RHF and UHF reference functions? -*/
    options.add_bool("MP2_AMPS_PRINT", 0);
    /*- Do print MP2 and CCSD pair energies for RHF references? -*/
    options.add_bool("PAIR_ENERGIES_PRINT", 0);
    /*- Do print spin-adapted pair energies? -*/
    options.add_bool("SPINADAPT_ENERGIES", false);
    /*- Do ? -*/
    options.add_bool("T3_WS_INCORE", 0);
    /*- Do ? -*/
    options.add_bool("SCSN_MP2", 0);
    /*- Do ? -*/
    options.add_bool("SCS_MP2", 0);
    /*- Do ? -*/
    options.add_bool("SCS_CCSD", 0);
    /*- -*/
    options.add_double("MP2_SCALE_OS",1.20);
    /*- -*/
    options.add_double("MP2_SCALE_SS",1.0/3.0);
    /*- -*/
    options.add_double("CC_SCALE_OS", 1.27);
    /*- -*/
    options.add_double("CC_SCALE_SS",1.13);
    /*- Convert ROHF MOs to semicanonical MOs -*/
    options.add_bool("SEMICANONICAL", true);
  }
  if(name == "CIS"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "CIS", "CCSD CCSD_T EOM_CCSD CIS");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF", "RHF ROHF UHF");
    /*- -*/
    options.add_double("LOCAL_AMPS_PRINT_CUTOFF", 0.60);
    /*- Maximum number of iterations -*/
    options.add_int("MAXITER", 500);
    /*- Convergence criterion for CIS wavefunction.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("R_CONVERGENCE", 1e-7);
    /*- -*/
    options.add("STATES_PER_IRREP", new ArrayType());
    /*- -*/
    options.add_str("DIAG_METHOD", "DAVIDSON", "DAVIDSON FULL");
    /*- Do ? -*/
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
    /*- Do ? -*/
    options.add_bool("DOMAIN_PRINT", 0);
  }
  if(name == "LMP2"|| options.read_globals()) {
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "LMP2");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF", "RHF");
    /*- Auxiliary basis set for MP2 density fitting calculations -*/
    options.add_str("DF_BASIS_MP2", "");
    /*- Do use density fitting? Turned on with specification of fitting basis. -*/
    if(options.get_str("DF_BASIS_MP2") != "")
      options.add_bool("DF_LMP2", true);
    else
      options.add_bool("DF_LMP2", false);
    /*- Maximum number of iterations -*/
    options.add_int("MAXITER", 50);
    /*- Convergence criterion for energy (change).
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("E_CONVERGENCE", 1e-7);
    /*- Convergence criterion for T2 amplitudes (RMS change). 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("R_CONVERGENCE", 1e-5);
    /*- Minimum absolute value below which parts of the Fock matrix are skipped.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("FOCK_TOLERANCE", 1e-2);
    /*- Do use DIIS extrapolation to accelerate convergence? -*/
    options.add_bool("DIIS", 1);
    /*- Do neglect distant pairs? -*/
    options.add_bool("NEGLECT_DISTANT_PAIR", 1);
    /*-  Distant pair cutoff -*/
    options.add_double("DISTANT_PAIR_CUTOFF", 8.0);
    /*- Iteration at which to start DIIS extrapolation -*/
    options.add_int("DIIS_START_ITER", 3);
    /*- Maximum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("DIIS_MAX_VECS", 5);
    /*- Localization cutoff -*/
    options.add_double("LOCAL_CUTOFF", 0.02);
    /*- -*/
    options.add_int("MEMORY", 2000);
    /*- Do ? -*/
    options.add_bool("SCS", false);
    /*- Do ? -*/
    options.add_bool("SCS_N", false);
    /*- -*/
    options.add_double("SCALE_OS", 6.0/5.0);
    /*- -*/
    options.add_double("SCALE_SS", 1.0/3.0);
    /*- Do screen integrals? -*/
    options.add_bool("SCREEN_INTS", false);
    /*- Minimum absolute value below which integrals are neglected. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("INTS_TOLERANCE", 1e-7);
    /*- Do exit after printing the domains? -*/
    options.add_bool("DOMAIN_PRINT_EXIT", 0);
   }
  if(name=="DFMP2"|| options.read_globals()) {
    /*- -*/
    options.add_int("MADMP2_SLEEP", 0);
    /*- Primary basis set -*/
    options.add_str("BASIS","NONE");
    /*- Auxiliary basis set for MP2 density fitting computations. Defaults to BASIS-RI. -*/
    options.add_str("DF_BASIS_MP2","");
    /*- OS Scale -*/
    options.add_double("SCALE_OS", 6.0/5.0);
    /*- SS Scale  -*/
    options.add_double("SCALE_SS", 1.0/3.0);
    /*- \% of memory for DF-MP2 three-index buffers  -*/
    options.add_double("DFMP2_MEM_FACTOR", 0.9);
    /*- Minimum absolute value below which integrals are neglected. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("INTS_TOLERANCE", 0.0);
    /*- Number of threads to compute integrals with. 0 is wild card -*/
    options.add_int("DF_INTS_NUM_THREADS", 0);
  }
  if(name=="DFCC"|| options.read_globals()) {
    /*- Type of wavefunction -*/
    options.add_str("WAVEFUNCTION","MP2","MP2 MP3 CCD DRPA");
    /*- Primary basis set -*/
    options.add_str("BASIS","NONE");
    /*- Minimum absolute value below which integrals are neglected. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("INTS_TOLERANCE", 0.0);
    /*- Convergence criterion for CC energy. 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("E_CONVERGENCE", 1e-8);
    /*- Convergence criterion for cluster amplitudes (RMS change). 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("R_CONVERGENCE", 1e-8);
    /*- Do use DIIS extrapolation to accelerate convergence? -*/
    options.add_bool("DIIS",true);
    /*- Minimum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("DIIS_MIN_VECS", 2);
    /*- Maximum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("DIIS_MAX_VECS", 6);
    /*- Maximum number iterations -*/
    options.add_int("MAXITER", 40);

    // => DF <= //

    /*- Auxiliary basis set for density fitting MO integrals. Defaults to BASIS-RI. -*/
    options.add_str("DF_BASIS_CC","NONE");
    /*- Fitting metric algorithm -*/
    options.add_str("FITTING_TYPE", "EIG", "EIG CHOLESKY QR");
    /*- Desired Fitting condition (inverse of max condition number) -*/
    options.add_double("FITTING_COND", 1.0E-10);

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
    /*- The DFT basis cutoff.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("PS_BASIS_TOLERANCE", 0.0);
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
    /*- Do use range-separation procedure in PS? -*/
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
    /*- Do continue RPA even if T's are not numerically SPD? -*/
    options.add_bool("RPA_RISKY",false);
    /*- Continue RPA numerical SPD Tolerance -*/
    options.add_double("RPA_PLUS_EPSILON",1.0E-12);
    /*- RPA alpha parameter -*/
    options.add_double("RPA_ALPHA", 1.0);

  }
  if(name == "PSIMRCC"|| options.read_globals()) {
    /*- The multiplicity, $M@@S(M@@S+1)$, of the target state.  Must be specified if different from the reference $M@@s$. -*/
      options.add_int("CORR_MULTP",1);
    /*- The molecular charge of the target state -*/
      options.add_int("CORR_CHARGE",0);
    /*- The amount (percentage) of damping to apply to the amplitude updates. 
        0 will result in a full update, 100 will completely stall the update. A
        value around 20 (which corresponds to 20\% of the amplitudes from the
        previous iteration being mixed into the current iteration)
        can help in cases where oscillatory convergence is observed. -*/
    options.add_double("DAMPING_PERCENTAGE",0.0);
    /*- Maximum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("DIIS_MAX_VECS",7);
    /*- Number of threads -*/
    options.add_int("NUM_THREADS",1);
    /*- Which root of the effective hamiltonian is the target state? -*/
    options.add_int("FOLLOW_ROOT",1);
    /*- Convergence criterion for energy.
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("E_CONVERGENCE",1e-9);
    /*- Convergence criterion for amplitudes (residuals).
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("R_CONVERGENCE",1e-9);
    /*- Maximum number of iterations to determine the amplitudes -*/
    options.add_int("MAXITER",100);
    /*- The number of DIIS vectors needed before extrapolation is performed -*/
    options.add_int("DIIS_START",2);
    /*- The shift to apply to the denominators ($\times$ 1000), {\it c.f.} Taube and Bartlett, JCP, 130, 144112 (2009) -*/
    options.add_int("TIKHONOW_OMEGA",0);  // Omega = TIKHONOW_OMEGA / 1000
    /*- The cycle after which Tikhonow regularization is stopped.  Set to zero to allow regularization in all iterations -*/
    options.add_int("TIKHONOW_MAX",5);
    /*- Do use DIIS extrapolation to accelerate convergence for iterative triples excitations? -*/
    options.add_bool("TRIPLES_DIIS",false);
    /*- Do lock onto a singlet root? -*/
    options.add_bool("LOCK_SINGLET",false);
    /*- Do start from a MP2 guess? -*/
    options.add_bool("MP2_GUESS",true);
    /*- Do use the averaged Fock matrix over all references in (T) computations? -*/
    options.add_bool("FAVG_CCSD_T",false);
    /*- Do include the fourth-order contributions to the effective Hamiltonian? -*/
    options.add_bool("HEFF4",true);
    /*- Do include the off-diagonal corrections in (T) computations? -*/
    options.add_bool("OFFDIAGONAL_CCSD_T",true);
    /*- Do include the diagonal corrections in (T) computations? -*/
    options.add_bool("DIAGONAL_CCSD_T",true);
    /*- Do diagonalize the effective Hamiltonian? -*/
    options.add_bool("DIAGONALIZE_HEFF",false);
    /*- Do use symmetry to map equivalent determinants onto each other, for efficiency? -*/
    options.add_bool("USE_SPIN_SYM",true);
    /*- Do zero the internal amplitudes, i.e., those that map reference determinants onto each other? -*/
    options.add_bool("ZERO_INTERNAL_AMPS",true);
    /*- Do include the terms that couple the reference determinants? -*/
    options.add_bool("COUPLING_TERMS",true);
    /*- Do print the effective Hamiltonian? -*/
    options.add_bool("HEFF_PRINT",false);
    /*- Do compute the perturbative corrections for basis set incompleteness? !expert -*/
    options.add_bool("PERTURB_CBS",false);
    /*- Do include the terms that couple different reference determinants in
        perturbative CBS correction computations? !expert -*/
    options.add_bool("PERTURB_CBS_COUPLING",true);
    /*- Do use Tikhonow regularization in (T) computations? !expert -*/
    options.add_bool("TIKHONOW_TRIPLES",false);
    /*- The type of perturbation theory computation to perform -*/
    options.add_str("PT_ENERGY","SECOND_ORDER","SECOND_ORDER SCS_SECOND_ORDER PSEUDO_SECOND_ORDER SCS_PSEUDO_SECOND_ORDER");
    /*- The type of correlated wavefunction -*/
    options.add_str("CORR_WFN","CCSD","PT2 CCSD MP2-CCSD CCSD_T");
    /*- The type of CCSD(T) computation to perform -*/
    options.add_str("CORR_CCSD_T","STANDARD","STANDARD PITTNER");
    /*- Reference wavefunction type used in MRCC computations -*/
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
    /*- The number of frozen occupied orbitals per irrep -*/
    options.add("FROZEN_DOCC", new ArrayType());
    /*- The number of doubly occupied orbitals per irrep -*/
    options.add("RESTRICTED_DOCC", new ArrayType());
    /*- The number of active orbitals per irrep -*/
    options.add("ACTIVE", new ArrayType());
    /*- The number of frozen virtual orbitals per irrep -*/
    options.add("FROZEN_UOCC", new ArrayType());
    /*- -*/
    options.add_int("SMALL_CUTOFF", 0);
    /*- Do ? -*/
    options.add_bool("NO_SINGLES", false);
  }
  if(name == "OPTKING"|| options.read_globals()) {
      /*- Specifies minimum search, transition-state search, or IRC following; allowed values = {MIN, TS, IRC} -*/
      options.add_str("OPT_TYPE", "MIN", "MIN TS IRC");
      /*- Whether to do a Newton-Raphson step, or an RFO step; allowed values = {NR, RFO} -*/
      options.add_str("STEP_TYPE", "RFO", "RFO NR SD");
      /*- Initial maximum step size in bohr or radian along an internal coordinate {double} -*/
      options.add_double("INTRAFRAGMENT_STEP_LIMIT", 0.4);
      /*- Lower bound for dynamic trust radius in au {double} -*/
      options.add_double("INTRAFRAGMENT_STEP_LIMIT_MIN", 0.001);
      /*- Upper bound for dynamic trust radius in au {double} -*/
      options.add_double("INTRAFRAGMENT_STEP_LIMIT_MAX", 1.0);
      /*- Do 'follow' the initial RFO vector after the first step? -*/
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
      /*- Do only generate the internal coordinates and then stop? -*/
      options.add_bool("GENERATE_INTCOS_ONLY", false);
      /*- What model Hessian to use to guess intrafragment force constants {SCHLEGEL, FISCHER} -*/
      options.add_str("INTRAFRAGMENT_H", "FISCHER", "FISCHER SCHLEGEL LINDH SIMPLE");
      /*- Whether to use the default of FISCHER_LIKE force constants for the initial guess {DEFAULT, FISCHER_LIKE} -*/
      options.add_str("INTERFRAGMENT_H", "DEFAULT", "DEFAULT FISCHER_LIKE");
      /*- Do freeze all fragments rigid? -*/
      options.add_bool("FREEZE_INTRAFRAGMENT", false);
      /*- Do add bond coordinates at nearby atoms for non-bonded systems? -*/
      options.add_bool("ADD_AUXILIARY_BONDS", false);
      /*- Do save and print the geometry from the last projected step at the end
      of a geometry optimization? Otherwise (and by default), save and print
      the previous geometry at which was computed the gradient that satisfied 
      the convergence criteria. -*/
      options.add_bool("FINAL_GEOM_WRITE", false);
      /*- Choose from supported Hessian updates {NONE, BFGS, MS, POWELL, BOFILL} -*/
      options.add_str("H_UPDATE", "BFGS", "NONE BFGS MS POWELL BOFILL");
      /*-  How many previous steps' data to use in Hessian update; 0=use them all ; {integer} -*/
      options.add_int("H_UPDATE_USE_LAST", 6);
      /*- Do limit the magnitude of changes caused by the Hessian update? -*/
      options.add_bool("H_UPDATE_LIMIT", true);
      /*- If the above is true, changes to the Hessian from the update are limited to the larger of
          (H_update_limit_scale)*(the previous value) and H_update_limit_max (in au). -*/
      options.add_double("H_UPDATE_LIMIT_MAX", 1.00);
      /*- If the above is true, changes to the Hessian from the update are limited to the larger of
          (H_update_limit_scale)*(the previous value) and H_update_limit_max (in au). -*/
      options.add_double("H_UPDATE_LIMIT_SCALE", 0.50);
      /*- Do use $\frac{1}{R@@{AB}}$ for the stretching coordinate between fragments? Otherwise, use $R@@{AB}$. -*/
      options.add_bool("INTERFRAGMENT_DISTANCE_INVERSE", false);
      /*- For now, this is a general maximum distance for the definition of H-bonds -*/
      options.add_double("MAXIMUM_H_BOND_DISTANCE", 4.3);
      /*- QCHEM optimization criteria: maximum force. See the note at the beginning of Section \ref{keywords}. -*/
      options.add_double("MAX_FORCE_G_CONVERGENCE", 3.0e-4);
      /*- QCHEM optimization criteria: rms force. See the note at the beginning of Section \ref{keywords}. -*/
      options.add_double("RMS_FORCE_G_CONVERGENCE", 3.0e-4);
      /*- QCHEM optimization criteria: maximum energy change. See the note at the beginning of Section \ref{keywords}. -*/
      options.add_double("MAX_ENERGY_G_CONVERGENCE", 1.0e-6);
      /*- QCHEM optimization criteria: maximum displacement. See the note at the beginning of Section \ref{keywords}. -*/
      options.add_double("MAX_DISP_G_CONVERGENCE", 1.2e-3);
      /*- QCHEM optimization criteria: rms displacement. See the note at the beginning of Section \ref{keywords}. -*/
      options.add_double("RMS_DISP_G_CONVERGENCE", 1.2e-3);
      /*- Do test B matrix? -*/
      options.add_bool("TEST_B", false);
      /*- Do test derivative B matrix? -*/
      options.add_bool("TEST_DERIVATIVE_B", false);
      /*- Do read Cartesian Hessian? -*/
      options.add_bool("READ_CARTESIAN_H", false);
      /*- Define IRC step size in bohr(amu)$^{1/2}$ -*/
      options.add_double("IRC_STEP_SIZE", 0.2);
      /*- Define IRC mapping direction {FORWARD, BACKWARD} -*/
      options.add_str("IRC_DIRECTION", "FORWARD", "FORWARD BACKWARD");
      /*- Set number of consecutive backward steps allowed in optimization -*/
      options.add_int("CONSECUTIVE_BACKSTEPS", 1);
  }
  if(name == "FINDIF"|| options.read_globals()) {
      /*- Number of points for finite-differences (3 or 5) -*/
      options.add_int("POINTS", 3); // Can we error check integers?
      /*- Displacement size in au for finite-differences. -*/
      options.add_double("DISP_SIZE", 0.005);
  }
  if (name == "OMP2"|| options.read_globals()) {
    //options.add_int("MEMORY", 256);
    //options.add_str("REFERENCE", "UHF", "UHF");
    
    /*- 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("E_CONVERGENCE",1e-8);
    /*-
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("R_CONVERGENCE",1e-5);
    /*-
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("RMS_MOGRAD_CONVERGENCE",1e-5);
    /*- 
    See the note at the beginning of Section \ref{keywords}. -*/
    options.add_double("MAX_MOGRAD_CONVERGENCE",1e-4);
    /*- -*/
    options.add_int("CC_MAXITER",50);
    /*- -*/
    options.add_int("MO_MAXITER",50);
    /*- -*/
    options.add_int("CACHELEVEL",2);
    /*- Number of vectors used in DIIS -*/
    options.add_int("DIIS_MAX_VECS",4);
    /*- -*/
    options.add_int("CUTOFF",14);
    
    /*- -*/
    options.add_double("MO_STEP_MAX",0.5);
    /*- -*/
    options.add_double("LSHIFT_PARAMETER",0.02);
    /*- -*/
    options.add_double("OS_SCALE",1.2);
    /*- -*/
    options.add_double("SS_SCALE",0.3333333333333333);
    /*- -*/
    options.add_double("SOS_SCALE",1.3); // It is used for MP2 (for SOS-MP2 recommended value is 1.3, but for SOS-OO-MP2 (O2) it is 1.2)
    /*- -*/
    options.add_double("SOS_SCALE2",1.2); // It is used for OMP2 (for SOS-MP2 recommended value is 1.3, but for SOS-OO-MP2 (O2) it is 1.2)
    /*- Do ? -*/ 
    options.add_bool("LEVEL_SHIFT",true);
    //options.add_str("LINEQ","CDGESV","CDGESV FLIN POPLE");
    options.add_str("ORTH_TYPE","MGS","GS MGS");
    //options.add_str("STABILITY","FALSE","TRUE FALSE");
    /*- Do ? -*/ 
    options.add_bool("NAT_ORBS",false);
    /*- -*/
    options.add_str("OPT_METHOD","DIIS","SD DIIS");    
    /*- -*/
    options.add_str("HESS_TYPE","NONE","NONE");    
    /*- Do ? -*/ 
    options.add_bool("OMP2_ORBS_PRINT",false);  
    /*- Do ? -*/ 
    options.add_bool("DO_SCS",false); 
    /*- Do ? -*/ 
    options.add_bool("DO_SOS",false); 
    /*- Do write coefficient matrices to psi files? -*/
    options.add_bool("MO_WRITE",false);   
    /*- Do read coefficient matrices from psi files? -*/ 
    options.add_bool("MO_READ",false);
  }
  return true;
}

} //end ::psi

