/*! \file read_calculation_options
    \defgroup PSI4
*/

#include <liboptions/liboptions.h>
#include <liboptions/python.h>
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

  // name == "GLOBALS" fake line to make document_options_and_tests.pl generate a GLOBALS doc section

  /*- Units used in geometry specification -*/
  options.add_str("UNITS", "ANGSTROMS", "BOHR AU A.U. ANGSTROMS ANG ANGSTROM");

  /*- An array containing the number of doubly-occupied orbitals per irrep
  (in Cotton order) -*/
  options.add("DOCC", new ArrayType());
  /*- An array containing the number of singly-occupied orbitals per irrep
  (in Cotton order) -*/
  options.add("SOCC", new ArrayType());
  /*- An array containing the number of frozen doubly-occupied orbitals per
  irrep (these are not excited in a correlated wavefunction, nor can they be
  optimized in MCSCF -*/
  options.add("FROZEN_DOCC", new ArrayType());
  /*- An array containing the number of frozen unoccupied orbitals per
  irrep (these are not populated in a correlated wavefunction, nor can they be
  optimized in MCSCF -*/
  options.add("FROZEN_UOCC", new ArrayType());
  /*- The number of core orbitals to freeze in later correlated computations.
  FROZEN_DOCC trumps this option -*/
  options.add_int("NUM_FROZEN_DOCC", 0);
  /*- The number of virtual orbitals to freeze in later correlated computations.
  FROZEN_UOCC trumps this option -*/
  options.add_int("NUM_FROZEN_UOCC", 0);
  /*- Specifies how many core orbitals to freeze in correlated computations.
  ``TRUE`` will default to freezing the standard default number of core orbitals.
  For heavier elements, there can be some ambiguity in how many core
  orbitals to freeze; in such cases, ``SMALL`` picks the most conservative
  standard setting (freezes fewer orbitals), and ``LARGE`` picks the least
  conservative standard setting (freezes more orbitals).  More precise
  control over the number of frozen orbitals can be attained by using
  the keywords |globals__num_frozen_docc| (gives the total number of orbitals to
  freeze, program picks the lowest-energy orbitals) or |globals__frozen_docc| (gives
  the number of orbitals to freeze per irreducible representation) -*/
  options.add_str("FREEZE_CORE","FALSE", "FALSE TRUE SMALL LARGE");

  /*- Do use pure angular momentum basis functions?
  If not explicitly set, the default comes from the basis set. -*/
  options.add_bool("PUREAM", true);
  /*- The amount of information to print to the output file.  1 prints
  basic information, and higher levels print more information. A value
  of 5 will print very large amounts of debugging information. -*/
  options.add_int("PRINT", 1);
  /*- The amount of information to print to the output file !expert -*/
  options.add_int("DEBUG", 0);
  /*- Some codes (DFT) can dump benchmarking data to separate output files -*/
  options.add_int("BENCH", 0);
  /*- Wavefunction type !expert -*/
  options.add_str("WFN", "SCF");
  /*- Derivative level !expert -*/
  options.add_str("DERTYPE", "NONE", "NONE FIRST SECOND RESPONSE");
  /*- Number of columns to print in calls to ``Matrix::print_mat``. !expert -*/
  options.add_int("MAT_NUM_COLUMN_PRINT", 5);
  /*- List of properties to compute -*/
  options.add("PROPERTIES", new ArrayType());

  /*- PSI4 dies if energy does not converge. !expert -*/
  options.add_bool("DIE_IF_NOT_CONVERGED", true);

  // CDS-TODO: We should go through and check that the user hasn't done
  // something silly like specify frozen_docc in DETCI but not in TRANSQT.
  // That would create problems.  (This was formerly checked in DETCI
  // itself, but I don't think DETCI will have the info available to check
  // this anymore).  This problem has affected users in the past.
  // Same goes for restricted_docc, restricted_uocc, ras1, ras2, ras3,
  // frozen_uocc.

  if (name == "DETCI" || options.read_globals()) {
    /*- MODULEDESCRIPTION Performs configuration interaction (CI)
    computations of various types, including restricted-active-space
    (RAS) CI, full CI, the CI component of multi-configuration
    self-consistent-field (MCSCF) and complete-active-space 
    self-consistent-field (CASSCF) computations, and arbitrary-order
    perturbation theory and arbitrary-order coupled-cluster
    computations for small molecules. -*/

    /*- SUBSECTION General Options -*/

    /*- Wavefunction type.  This should be set automatically from
    the calling Psithon function.  !expert -*/
    options.add_str("WFN", "DETCI", "DETCI CI ZAPTN DETCAS CASSCF RASSCF");

    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF", "RHF ROHF");

    /*- Convergence criterion for CI residual vector in the Davidson 
    algorithm (RMS error).
    The default is 1e-4 for energies and 1e-7 for gradients. -*/
    options.add_double("R_CONVERGENCE", 1e-4);

    /*- Convergence criterion for energy. -*/
    options.add_double("E_CONVERGENCE", 1e-6);

    /*- Maximum number of iterations to diagonalize the Hamiltonian -*/
    options.add_int("MAXITER", 12);

    /*- Do a full CI (FCI)? If TRUE, overrides the value of |detci__ex_level|. -*/
    options.add_bool("FCI",false);

    /*- The CI excitation level -*/
    options.add_int("EX_LEVEL", 2);

    /*- In a RAS CI, this is the additional excitation level for allowing
    electrons out of RAS I into RAS II.  The maximum number of holes in
    RAS I is therefore |detci__ex_level| + VAL_EX_LEVEL. -*/
    options.add_int("VAL_EX_LEVEL", 0);

    /*- number of CI roots to find -*/
    options.add_int("NUM_ROOTS", 1);

    /*- Do stop DETCI after string information is formed
    and before integrals are read? -*/
    options.add_bool("ISTOP",false);

    /*- Do print a summary of the CI blocks? -*/
    options.add_bool("CIBLKS_PRINT",false);

    /*- Number of important determinants to print -*/
    options.add_int("NUM_DETS_PRINT",20);

    /*- Do freeze core orbitals? -*/
    // CDS-TODO: Need to make DETCI compatible with normal FREEZE_CORE
    options.add_bool("DETCI_FREEZE_CORE",true);

    /*- Do calculate the value of $\langle S^2\rangle$ for each root? -*/
    options.add_bool("S_SQUARED",false);

    /*- Specifies how to handle buffering of CI vectors.  A value of 0
    makes the program perform I/O one RAS subblock at a time; 1
    uses entire CI vectors at a time; and 2 uses one irrep block
    at a time.  Values of 0 or 2 cause some inefficiency in the I/O
    (requiring multiple reads of the C vector when constructing
    H in the iterative subspace if |detci__diag_method| = SEM), but require
    less core memory. -*/
    options.add_int("ICORE", 1);

    /*- Number of threads for DETCI. -*/
    options.add_int("CI_NUM_THREADS", 1);

    /*- Do print the sigma overlap matrix?  Not generally useful.  !expert -*/
    options.add_bool("SIGMA_OVERLAP", false);

    /*- Array giving the root numbers of the states to average in a
    state-averaged procedure such as SA-CASSCF. Root numbering starts
    from 1. -*/
    options.add("AVG_STATES", new ArrayType());

    /*- Array giving the weights for each state in a state-averaged
    procedure -*/
    // CDS:TODO - Does this work for doubles??
    options.add("AVG_WEIGHTS", new ArrayType());


    /*- SUBSECTION Specifying the CI Space -*/

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

    /*- The value of the spin quantum number $S$ is given by this option.
    The default is determined by the value of the multiplicity.  This is used
    for two things: (1) determining the phase of the redundant half of the CI
    vector when the $M@@s = 0$ component is used (i.e., |detci__ms0| = ``TRUE``), and (2) making
    sure the guess vector has the desired value of $\langle S^2\rangle$ 
    (if |detci__s_squared| is ``TRUE`` and |detci__icore| = ``1``). -*/
    options.add_double("S", 0.0);

    /*- Do use the $M@@s = 0$ component of the state? Defaults to TRUE
    if closed-shell and FALSE otherwise. Related to the |detci__s| option. -*/
    options.add_bool("MS0",false);

    /*- An array of length |detci__ex_level| specifying whether each excitation type
    (S,D,T, etc.) is allowed (1 is allowed, 0 is disallowed).  Used to
    specify non-standard CI spaces such as CIST.  !expert -*/
    options.add("EX_ALLOW", new ArrayType());

    /*- Do eliminate determinants not valid for spin-complete spin-flip CI's?
    [see J. S. Sears et al, J. Chem. Phys. 118, 9084-9094 (2003)] !expert -*/
    options.add_bool("SF_RESTRICT", false);

    /*- maximum number of alpha electrons in RAS III -*/
    options.add_int("A_RAS3_MAX",-1);

    /*- maximum number of beta electrons in RAS III -*/
    options.add_int("B_RAS3_MAX",-1);

    /*- maximum number of electrons in RAS III -*/
    options.add_int("RAS3_MAX",-1);

    /*- maximum number of electrons in RAS IV -*/
    options.add_int("RAS4_MAX",-1);

    /*- maximum number of electrons in RAS III + IV -*/
    options.add_int("RAS34_MAX",-1);

    /*- Do allow "mixed" RAS II/RAS III excitations into the CI space?
    If FALSE, then if there are any electrons
    in RAS III, then the number of holes in RAS I cannot exceed the given
    excitation level |detci__ex_level|. !expert -*/
    options.add_bool("MIXED",true);

    /*- Do allow "mixed" excitations involving RAS IV into the CI space.
    Useful to specify a split-virtual
    CISD[TQ] computation.  If FALSE, then if there are any electrons
    in RAS IV, then the number of holes in RAS I cannot exceed the given
    excitation level |detci__ex_level|.  !expert -*/
    options.add_bool("MIXED4",true);

    /*- Do restrict strings with $e-$ in RAS IV?  Useful to reduce the number 
    of strings required if MIXED4=true, as in a split-virutal CISD[TQ]
    computation.  If more than one electron is in RAS IV, then the
    holes in RAS I cannot exceed the number of particles in
    RAS III + RAS IV (i.e., |detci__ex_level|), or else the string is discarded.
    !expert -*/
    options.add_bool("R4S",false);

    /*- SUBSECTION Diagonalization Methods -*/

    /*- This specifies which method is to be used in diagonalizing the
    Hamiltonian.  The valid options are: ``RSP``, to form the entire H
    matrix and diagonalize using libciomr to obtain all eigenvalues
    (n.b. requires HUGE memory); ``OLSEN``, to use Olsen's preconditioned
    inverse subspace method (1990); ``MITRUSHENKOV``, to use a 2x2
    Olsen/Davidson method; and ``DAVIDSON`` (or ``SEM``) to use Liu's
    Simultaneous Expansion Method, which is identical to the Davidson method
    if only one root is to be found.  There also exists a SEM debugging mode,
    ``SEMTEST``.  The ``SEM`` method is the most robust, but it also
    requires $2NM+1$ CI vectors on disk, where $N$ is the maximum number of
    iterations and $M$ is the number of roots. -*/
    options.add_str("DIAG_METHOD", "SEM", "RSP OLSEN MITRUSHENKOV DAVIDSON SEM SEMTEST");

    /*- This specifies the type of preconditioner to use in the selected
    diagonalization method.  The valid options are: ``DAVIDSON`` which
    approximates the Hamiltonian matrix by the diagonal elements;
    ``H0BLOCK_INV`` which uses an exact Hamiltonian of |detci__h0_blocksize| and
    explicitly inverts it; ``GEN_DAVIDSON`` which does a spectral
    decomposition of H0BLOCK; ``ITER_INV`` using an iterative approach
    to obtain the correction vector of H0BLOCK.  The ``H0BLOCK_INV``, ``GEN_DAVIDSON``,
    and ``ITER_INV`` approaches are all formally equivalent but the ``ITER_INV`` is
    less computationally expensive.  Default is ``DAVIDSON``. -*/
    options.add_str("PRECONDITIONER", "DAVIDSON", "LANCZOS DAVIDSON GEN_DAVIDSON H0BLOCK H0BLOCK_INV ITER_INV H0BLOCK_COUPLING EVANGELISTI");

    /*- The update or correction vector formula, either ``DAVIDSON`` (default)
    or ``OLSEN``. -*/
    options.add_str("UPDATE", "DAVIDSON", "DAVIDSON OLSEN");

    /*- How to average H diag energies over spin coupling sets.
    ``HD_EXACT`` uses the exact diagonal energies which results in expansion
    vectors which break spin symmetry. ``HD_KAVE`` averages the diagonal
    energies over a spin-coupling set yielding spin pure expansion vectors.
    ``ORB_ENER`` employs the sum of orbital energy approximation giving
    spin pure expansion vectors but usually doubles the number of Davidson
    iterations. ``EVANGELISTI`` uses the sums and differences of orbital
    energies with the SCF reference energy to produce spin pure expansion
    vectors. ``LEININGER`` approximation which subtracts the one-electron
    contribution from the orbital energies, multiplies by 0.5, and adds
    the one-electron contribution back in, producing spin pure expansion
    vectors and developed by Matt Leininger and works as well as
    ``EVANGELISTI``. !expert -*/
    options.add_str("HD_AVG", "EVANGELISTI", "EVANGELISTI HD_EXACT HD_KAVE ORB_ENER LEININGER Z_KAVE");

    /*- This parameter specifies the size of the H0 block of the Hamiltonian
    which is solved exactly.  The n determinants with the lowest SCF
    energy are selected, and a submatrix of the Hamiltonian is formed
    using these determinants.  This submatrix is used to accelerate
    convergence of the CI iterations in the OLSEN and MITRUSHENKOV
    iteration schemes, and also to find a good starting guess for the
    SEM method if |detci__guess_vector| is ``H0_BLOCK``.  Defaults to 400.
    Note that the program may change the given size for Ms=0 cases
    (|detci__ms0| is TRUE) if it determines that the H0 block includes only
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

    /*- Do use least-squares extrapolation in iterative solution of CI
    vector? -*/
    options.add_bool("LSE",false);

    /*- Number of iterations between least-squares extrapolations -*/
    options.add_int("LSE_COLLAPSE", 3);

    /*- Minimum converged energy for least-squares
    extrapolation to be performed -*/
    options.add_double("LSE_TOLERANCE", 3);


    /*- SUBSECTION Density Matrices -*/

    /*- Do compute one-particle density matrix if not otherwise required? -*/
    options.add_bool("OPDM", false);

    /*- Do compute two-particle density matrix if not otherwise required? -*/
    options.add_bool("TPDM", false);

    /*- Do print the one-particle density matrix for each root? -*/
    options.add_bool("OPDM_PRINT", false);

    /*- Do write the natural orbitals? -*/
    options.add_bool("NAT_ORBS_WRITE", false);

    /*- Do average the OPDM over several roots in
    order to obtain a state-average one-particle density matrix?  This
    density matrix can be diagonalized to obtain the CI natural orbitals. -*/
    options.add_bool("OPDM_AVG", false);

    /*- Sets the root number for which CI natural orbitals are written
    to PSIF_CHKPT.  The default value is 1 (lowest root). -*/
    options.add_int("NAT_ORBS_WRITE_ROOT", 1);

    /*- Do compute the kinetic energy contribution from the correlated part of
    the one-particle density matrix? !expert -*/
    options.add_bool("OPDM_KE", false);

    /*- Do print the two-particle density matrix? (Warning: large tensor) -*/
    options.add_bool("TPDM_PRINT", false);

    /*- Do compute the transition density?  Note: only transition densities
    between roots of the same symmetry will be evaluated.  DETCI
    does not compute states of different irreps within the same
    computation; to do this, lower the symmetry of the computation.-*/
    options.add_bool("TDM", false);

    /*- Do write the transition density? -*/
    options.add_bool("TDM_WRITE", false);

    /*- Do print the transition density? -*/
    options.add_bool("TDM_PRINT", false);

    /*- Do compute the dipole moment? -*/
    options.add_bool("DIPMOM", false);


    /*- SUBSECTION Root Following -*/

    /*- The root to write out the two-particle density matrix for
    (the one-particle density matrices are written for all roots).
    Useful for a state-specific CASSCF or CI optimization on an
    excited state. -*/
    options.add_int("FOLLOW_ROOT", 1);

    /*- In following a particular root (see |detci__follow_root|), sometimes the
    root number changes.  To follow a root of a particular character,
    one can specify a list of determinants and their coefficients,
    and the code will follow the root with the closest overlap.  The
    user specifies arrays containing the absolute alpha string indices
    (A_i below), absolute beta indices (B_i below), and CI coefficients
    (C_i below) to form the desired vector.
    The format is FOLLOW_VECTOR = [ [[A_1, B_1], C_1], [[A_2, B_2], C_2], ...].
    !expert -*/
    options.add("FOLLOW_VECTOR", new ArrayType());

    /*- SUBSECTION Guess Vectors -*/

    /*- Guess vector type.  Accepted values are ``UNIT`` for a unit vector
    guess (|detci__num_roots| and |detci__num_init_vecs| must both be 1); ``H0_BLOCK`` to use
    eigenvectors from the H0 BLOCK submatrix (default); ``DFILE`` to use
    NUM_ROOTS previously converged vectors in the D file; ``IMPORT`` to
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

    /*- Do restart a DETCI iteration that
    terminated prematurely? It assumes that the CI and sigma vectors are on
    disk; the number of vectors specified by RESTART_VECS (obsolete) is collapsed
    down to one vector per root. -*/
    options.add_bool("RESTART",false);

    /*- Do invoke the FILTER_GUESS options that are used to filter out some 
    trial vectors which may not have the appropriate phase convention
    between two determinants?  This is useful to remove, e.g.,
    delta states when a sigma state is desired.  The user
    inputs two determinants (by giving the absolute alpha string
    number and beta string number for each), and also the
    desired phase between these two determinants for guesses
    which are to be kept.  FILTER_GUESS = TRUE turns on the filtering
    routine.  Requires additional keywords |detci__filter_guess_det1|,
    |detci__filter_guess_det2|, and |detci__filter_guess_sign|. !expert -*/
    options.add_bool("FILTER_GUESS", false);

    /*- The required phase (1 or -1) between the two determinants specified
    by |detci__filter_guess_det1| and |detci__filter_guess_det2|. !expert -*/
    options.add_int("FILTER_GUESS_SIGN", 1);

    /*- Array specifying the absolute alpha string number and beta string
    number for the first determinant in the filter procedure.
    (See |detci__filter_guess|).  !expert -*/
    options.add("FILTER_GUESS_DET1", new ArrayType());

    /*- Array specifying the absolute alpha string number and beta string
    number for the second determinant in the filter procedure.
    (See |detci__filter_guess|).  !expert -*/
    options.add("FILTER_GUESS_DET2", new ArrayType());

    /*- If present, the code will try to filter out a particular determinant
    by setting its CI coefficient to zero.  FILTER_ZERO_DET = [alphastr,
    betastr] specifies the absolute alpha and beta string numbers of the
    target determinant. This could be useful for trying to exclude states
    that have a nonzero CI coefficient for the given determinant.  However,
    this option was experimental and may not be effective.  !expert -*/
    options.add("FILTER_ZERO_DET", new ArrayType());


    /*- SUBSECTION File Handling -*/

    /*- Maximum number of Davidson subspace vectors which can
    be held on disk for the CI coefficient and sigma vectors.  (There
    is one H(diag) vector and the number of D vectors is equal to the
    number of roots).  When the number of vectors on disk reaches
    the value of MAX_NUM_VECS, the Davidson subspace will be
    collapsed to |detci__collapse_size| vectors for each root.  This is very
    helpful for saving disk space.  Defaults to |detci__maxiter| * |detci__num_roots|
    + |detci__num_init_vecs|. -*/
    options.add_int("MAX_NUM_VECS", 0);

    /*- Gives the number of vectors to retain when the Davidson subspace is
    collapsed (see |detci__max_num_vecs|).  If greater than one, the
    collapsed subspace retains the best estimate of the CI vector for
    the previous n iterations.   Defaults to 1. -*/
    options.add_int("COLLAPSE_SIZE", 1);

    /*- Do store converged vector(s) at the end of
    the run?  The vector(s) is(are) stored in a transparent format such that
    other programs can use it easily. The format is specified in
    :source:`src/lib/libqt/slaterdset.h` . -*/
    options.add_bool("VECS_WRITE", false);

    /*- Number of vectors to export -*/
    options.add_int("NUM_VECS_WRITE", 1);

    /*- Do compute the diagonal elements of the Hamiltonian matrix
    on-the-fly? Otherwise, a diagonal element vector is written
    to a separate file on disk. !expert -*/
    options.add_bool("HD_OTF",true);

    /*- Do use the last vector space in the BVEC file to write
    scratch DVEC rather than using a separate DVEC file? (Only
    possible if |detci__num_roots| = 1.) !expert -*/
    options.add_bool("NO_DFILE",false);

    /*- SUBSECTION General-Order Perturbation Theory -*/

    /*- Do compute the MPn series out to
    kth order where k is determined by |detci__max_num_vecs| ?  For open-shell systems
    (|detci__reference| is ROHF, |detci__wfn| is ZAPTN), DETCI will compute the ZAPTn series.
    |detci__guess_vector| must be set to UNIT, |detci__hd_otf| must be set to TRUE, and
    |detci__hd_avg| must be set to orb_ener; these should happen by default for
    MPN = TRUE. -*/
    options.add_bool("MPN",false);

    /*- If 0, save the MPn energy; if 1, save the MP(2n-1) energy (if
    available from |detci__mpn_wigner| = true); if 2, save the MP(2n-2) energy (if
    available from |detci__mpn_wigner| = true). !expert -*/
    options.add_int("MPN_ORDER_SAVE",0);

    /*- Do employ an orthonormal vector space rather than
      storing the kth order wavefunction? !expert -*/
    options.add_bool("MPN_SCHMIDT",false);

    /*- Do use Wigner formulas in the $E_{text{mp}n}$ series? !expert -*/
    options.add_bool("MPN_WIGNER",true);

    /*- The magnitude of perturbation $z$ in $H = H@@0 + z H@@1$ !expert -*/
    options.add_double("PERTURB_MAGNITUDE",1.0);


    /*- SUBSECTION General-Order Coupled-Cluster -*/

    /*- Do coupled-cluster computation? -*/
    options.add_bool("CC", false);

    /*- The CC excitation level -*/
    options.add_int("CC_EX_LEVEL", 2);

    /*- The CC valence excitation level -*/
    options.add_int("CC_VAL_EX_LEVEL", 0);

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

    /*- Number of important CC amplitudes per excitation level to print.
    CC analog to |detci__num_dets_print|. -*/
    options.add_int("NUM_AMPS_PRINT",10);

    /*- maximum number of alpha electrons in RAS III, for CC -*/
    options.add_int("CC_A_RAS3_MAX",-1);

    /*- maximum number of beta electrons in RAS III, for CC -*/
    options.add_int("CC_B_RAS3_MAX",-1);

    /*- maximum number of electrons in RAS III, for CC -*/
    options.add_int("CC_RAS3_MAX",-1);

    /*- maximum number of electrons in RAS IV, for CC -*/
    options.add_int("CC_RAS4_MAX",-1);

    /*- maximum number of electrons in RAS III + IV, for CC -*/
    options.add_int("CC_RAS34_MAX",-1);

    /*- Do export a CC vector to disk? -*/
    options.add_bool("CC_VECS_WRITE", false);

    /*- Do import a CC vector from disk? -*/
    options.add_bool("CC_VECS_READ", false);

    /*- Do fix amplitudes involving RAS I or RAS IV?  Useful in mixed
    MP2-CC methods. !expert -*/
    options.add_bool("CC_FIX_EXTERNAL", false);

    /*- Number of external indices before amplitude gets fixed by
    |detci__cc_fix_external|.  Experimental. !expert -*/
    options.add_int("CC_FIX_EXTERNAL_MIN", 1);

    /*- Do use variational energy expression in CC computation?
    Experimental.  !expert -*/
    options.add_bool("CC_VARIATIONAL", false);

    /*- Do ignore block if num holes in RAS I and II is $>$ cc_ex_lvl and if
    any indices correspond to RAS I or IV (i.e., include only all-active
    higher excitations)? !expert -*/
    options.add_bool("CC_MIXED", true);

    /*- Do update T amplitudes with orbital eigenvalues? (Usually would
    do this).  Not doing this is experimental.  !expert -*/
    options.add_bool("CC_UPDATE_EPS", true);

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


    /*- SUBSECTION Alternative Algorithms -*/

    /*- Do store strings specifically for FCI? (Defaults to TRUE for FCI.)
        !expert -*/
    options.add_bool("FCI_STRINGS",false);

    /*- Do string replacements on the fly in DETCI? Can
    save a gigantic amount of memory (especially for truncated CI's) but
    is somewhat flaky and hasn't been tested for a while.  It may work
    only works for certain classes of RAS calculations.  The current
    code is very slow with this option turned on. !expert -*/
    options.add_bool("REPL_OTF",false);

    /*- Do use some routines based on the papers of Bendazzoli et al. 
    to calculate sigma?  Seems to be slower and not worthwhile; may disappear
    eventually.  Works only for full CI and I don't remember if I could see
    how their clever scheme might be extended to RAS in general. !expert -*/
    options.add_bool("BENDAZZOLI", false);
  }

  if (name == "SAPT"|| options.read_globals()) {
    /*- MODULEDESCRIPTION Performs symmetry adapted perturbation theory (SAPT) 
    analysis to quantitatively analyze noncovalent interactions. -*/

    /*- The level of theory for SAPT -*/
    options.add_str("SAPT_LEVEL","SAPT0","SAPT0 SAPT2 SAPT2+ SAPT2+3");

    /*- Convergence criterion for energy (change) in the SAPT 
    $E@@{ind,resp}^{(20)}$ term during solution of the CPHF equations. -*/

    options.add_double("E_CONVERGENCE",1e-10);

    /*- Convergence criterion for residual of the CPHF coefficients in the SAPT
    $E@@{ind,resp}^{(20)}$ term. -*/
    options.add_double("D_CONVERGENCE",1e-8);

    /*- Don't solve the CPHF equations? Evaluate $E@@{ind}^{(20)}$ and
    $E@@{exch-ind}^{(20)}$ instead of their response-including coupterparts.
    Only turn on this option if the induction energy is not going to be 
    used. -*/
    options.add_bool("NO_RESPONSE",false);

    /*- Do use asynchronous disk I/O in the solution of the CPHF equations?
    Use may speed up the computation slightly at the cost of spawning an 
    additional thread. -*/
    options.add_bool("AIO_CPHF",false);

    /*- Do use asynchronous disk I/O in the formation of the DF integrals?
    Use may speed up the computation slightly at the cost of spawning an 
    additional thread. -*/
    options.add_bool("AIO_DF_INTS",false);

    /*- Maxmum number of CPHF iterations -*/
    options.add_int("MAXITER",50);
    /*- Do compute third-order corrections? !expert -*/
    options.add_bool("DO_THIRD_ORDER",false);
    /*- Do natural orbitals to speed up evaluation of the triples
    contribution to dispersion by truncating the virtual orbital space?
    Recommended true for all SAPT computations. -*/
    options.add_bool("NAT_ORBS",false);
    /*- Do use MP2 natural orbital approximations for the $v^4$ block of
    two-electron integrals in the evaluation of second-order T2 amplitudes?
    This approximation is promising for accuracy and computational savings,
    but it has not been rigorously tested. -*/
    options.add_bool("NAT_ORBS_T2",false);
    /*- Minimum occupation (eigenvalues of the MP2 OPDM) below which virtual
    natural orbitals are discarded for evaluating the triples contribution
    to dispersion. -*/
    options.add_double("OCC_TOLERANCE",1.0E-6);
    /*- Minimum absolute value below which all three-index DF integrals
    and those contributing to four-index integrals are neglected. The
    default is conservative, but there isn't much to be gained from
    loosening it, especially for higher-order SAPT. -*/
    options.add_double("INTS_TOLERANCE",1.0E-12);
    /*- Memory safety -*/
    options.add_double("SAPT_MEM_SAFETY",0.9);
    /*- Do force SAPT2 and higher to die if it thinks there isn't enough 
    memory?  Turning this off is ill-advised. -*/
    options.add_bool("SAPT_MEM_CHECK",true);
    /*- Primary basis set, describes the monomer molecular orbitals -*/
    options.add_str("BASIS", "");
    /*- Auxiliary basis set for SAPT density fitting computations.
    :ref:`Defaults <apdx:basisFamily>` to a RI basis. -*/
    options.add_str("DF_BASIS_SAPT", "");
    /*- Auxiliary basis set for SAPT Elst10 and Exch10 density fitting 
    computations, may be important if heavier elements are involved.
    Defaults to |sapt__df_basis_sapt|. -*/
    options.add_str("DF_BASIS_ELST", "");
    /*- Maximum error allowed (Max error norm in Delta tensor)
    in the approximate energy denominators employed for most of the
    $E@@{disp}^{(20)}$ and $E@@{exch-disp}^{(20)}$ evaluation. -*/
    options.add_double("DENOMINATOR_DELTA", 1.0E-6);
    /*- Denominator algorithm for PT methods. Laplace transformations
    are slightly more efficient. -*/
    options.add_str("DENOMINATOR_ALGORITHM", "LAPLACE", "LAPLACE CHOLESKY");
    /*- The scale factor used for opposite-spin pairs in SCS computations. 
    SS/OS decomposition performed for $E@@{disp}^{(20)}$ and 
    $E@@{exch-disp}^{(20)}$ terms. -*/
    options.add_double("SAPT_OS_SCALE", 6.0/5.0);
    /*- The scale factor used for same-spin pairs in SCS computations. SS/OS
    decomposition performed for $E@@{disp}^{(20)}$ and $E@@{exch-disp}^{(20)}$ 
    terms. -*/
    options.add_double("SAPT_SS_SCALE", 1.0/3.0);
    /*- The scope of core orbitals to freeze in evaluation of SAPT
    $E@@{disp}^{(20)}$ and $E@@{exch-disp}^{(20)}$ terms. Recommended true
    for all SAPT computations -*/
    options.add_str("FREEZE_CORE","FALSE", "FALSE TRUE SMALL LARGE");
    /*- The amount of information to print to the output file for the sapt
    module. For 0, only the header and final results are printed. For 1,
    (recommended for large calculations) some intermediate quantities are also 
    printed. -*/
    options.add_int("PRINT", 1);
  }

  if(name == "DCFT"|| options.read_globals()) {
      /*-MODULEDESCRIPTION Performs Density Cumulant Functional Theory 
      computations -*/

      /*- The algorithm to use for the density cumulant and orbital updates in the energy computation.
      Two-step algorithm (default) is generally more efficient and shows better convergence than simultaneous -*/
      options.add_str("ALGORITHM", "TWOSTEP", "TWOSTEP SIMULTANEOUS QC");
      /*- The algorithm to use for the solution of the response equations for the analytic gradients and properties.
      Two-step algorithm is generally more efficient than simultaneous and is used by default-*/
      options.add_str("RESPONSE_ALGORITHM", "TWOSTEP", "TWOSTEP SIMULTANEOUS");
      /*- Convergence criterion for the RMS of the residual vector in the density cumulant updates as well as
      the solution of the density cumulant and orbital response equations. In the orbital updates controls
      the RMS of the SCF error vector -*/
      options.add_double("R_CONVERGENCE", 1e-10);
      /*- Maximum number of density cumulant update micro-iterations per
      macro-iteration (for ALOGRITHM = TWOSTEP). Same keyword controls the
      maximum number of density cumulant response micro-iterations per
      macro-iteration for the solution of the response equations
      (for RESPONSE_ALOGRITHM = TWOSTEP) -*/
      options.add_int("LAMBDA_MAXITER", 50);
      /*- Maximum number of orbital update micro-iterations per
      macro-iteration (for ALOGRITHM = TWOSTEP). Same keyword controls the
      maximum number of orbital response micro-iterations per
      macro-iteration for the solution of the response equations
      (for RESPONSE_ALOGRITHM = TWOSTEP) -*/
      options.add_int("SCF_MAXITER", 50);
      /*- Maximum number of macro-iterations for both energy and the solution of the response equations -*/
      options.add_int("MAXITER", 40);
      /*- Value of RMS of the density cumulant residual and SCF error vector below which DIIS extrapolation starts.
      Same keyword controls the DIIS extrapolation for the solution of the response equations. -*/
      options.add_double("DIIS_START_CONVERGENCE", 1e-3);
      /*- Maximum number of error vectors stored for DIIS extrapolation -*/
      options.add_int("DIIS_MAX_VECS", 6);
      /*- Minimum number of error vectors stored for DIIS extrapolation -*/
      options.add_int("DIIS_MIN_VECS", 3);
      /*- Controls whether to avoid the AO->MO transformation of the two-electron integrals for the four-virtual case
      (<VV||VV>) by computing the corresponding terms in the AO basis. AO_BASIS = DISK algorithm reduces the memory
      requirements. It is, however, less efficient due to the extra I/O, so the default algorithm is preferred. -*/
      options.add_str("AO_BASIS", "NONE", "NONE DISK");
      /*- The amount (percentage) of damping to apply to the orbital update procedure:
      0 will result in a full update, 100 will completely stall the
      update. A value around 20 (which corresponds to 20\% of the previous 
      iteration's density being mixed into the current iteration)
      can help in cases where oscillatory convergence is observed. -*/
      options.add_double("DAMPING_PERCENTAGE",0.0);
      /*- The shift applied to the denominator in the density cumulant update iterations -*/
      options.add_double("TIKHONOW_OMEGA", 0.0);
      /*- Controls whether to compute the DCFT energy with the Tau^2 correction to Tau !expert-*/
      options.add_bool("TAU_SQUARED", false);
      /*- Controls whether to compute unrelaxed two-particle density matrix at the end of the energy computation !expert-*/
      options.add_bool("TPDM", false);
      /*- Controls whether to relax the orbitals during the energy computation or not (for debug puproses only).
      For practical applications only the default must be used !expert-*/
      options.add_bool("MO_RELAX", true);
      /*- Controls whether to ignore terms containing non-idempotent contribution to OPDM or not (for debug puproses only).
      For practical applications only the default must be used !expert-*/
      options.add_bool("IGNORE_TAU", false);
      /*- Controls how to cache quantities within the DPD library !expert-*/
      options.add_int("CACHELEVEL", 2);
      /*- Minimum absolute value below which integrals are neglected !expert-*/
      options.add_double("INTS_TOLERANCE", 1e-14);
      /*- Controls whether to force the occupation to be that of the SCF guess.
      For practical applications only the default must be used !expert-*/
      options.add_bool("LOCK_OCC", true);
      /*- Whether to read the orbitals from a previous computation, or to compute
          an MP2 guess !expert -*/
      options.add_str("DCFT_GUESS", "MP2", "CC BCC MP2");
      /*- Controls whether to relax the guess orbitals by taking the guess density cumulant
      and performing orbital update on the first macroiteration (for ALOGRITHM = TWOSTEP only) !expert-*/
      options.add_bool("RELAX_GUESS_ORBITALS", false);

  }
  if (name == "MINTS"|| options.read_globals()) {
      /*- MODULEDESCRIPTION Called at the beginning of SCF computations, 
      whenever disk-based molecular integrals are required. -*/

      /*- Primary basis set -*/
      options.add_str("BASIS","");
  }
  if (name == "SCF"|| options.read_globals()) {
    /*- MODULEDESCRIPTION Performs self consistent field (Hartree-Fock and 
    Density Functional Theory) computations.  These are the starting 
    points for most computations, so this code is called in most cases. -*/

    /*- SUBSECTION General Wavefunction Info -*/

    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "SCF", "SCF");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF", "RHF ROHF UHF CUHF RKS UKS");
    /*- Primary basis set -*/
    options.add_str("BASIS", "");
    /*- Auxiliary basis set for SCF density fitting computations. 
    :ref:`Defaults <apdx:basisFamily>` to a JKFIT basis. -*/
    options.add_str("DF_BASIS_SCF", "");
    /*- What algorithm to use for the SCF computation -*/
    options.add_str("SCF_TYPE", "PK", "DIRECT DF PK OUT_OF_CORE PS");
    /*- Keep JK object for later use? -*/
    options.add_bool("SAVE_JK", false);
    /*- Memory safety factor for allocating JK -*/
    options.add_double("SCF_MEM_SAFETY_FACTOR",0.75);
    /*- SO orthogonalization: symmetric or canonical? -*/
    options.add_str("S_ORTHOGONALIZATION","SYMMETRIC","SYMMETRIC CANONICAL");
    /*- Minimum S matrix eigenvalue to be used before compensating for linear 
    dependencies. -*/
    options.add_double("S_TOLERANCE",1E-7);
    /*- Minimum absolute value below which TEI are neglected. -*/
    options.add_double("INTS_TOLERANCE", 0.0);
    /*- The type of guess orbitals -*/
    options.add_str("GUESS", "CORE", "CORE GWH SAD READ");
    /*- The name of a molden-style output file which is only generated
    if the user specifies one -*/
    options.add_str("MOLDEN_FILE", "");
    /*- Flag to print the molecular orbitals. -*/
    options.add_bool("PRINT_MOS", false);
    /*- Flag to print the basis set. -*/
    options.add_bool("PRINT_BASIS", false);

    /*- SUBSECTION Convergence Control/Stabilization -*/

    /*- Maximum number of iterations -*/
    options.add_int("MAXITER", 100);
    /*- Convergence criterion for SCF energy. -*/
    options.add_double("E_CONVERGENCE", 1e-8);
    /*- Convergence criterion for SCF density. -*/
    options.add_double("D_CONVERGENCE", 1e-8);
    /*- The amount (percentage) of damping to apply to the early density updates.
        0 will result in a full update, 100 will completely stall the update.  A
        value around 20 (which corresponds to 20\% of the previous iteration's
        density being mixed into the current density)
        could help to solve problems with oscillatory convergence. -*/
    options.add_double("DAMPING_PERCENTAGE", 100.0);
    /*- The density convergence threshold after which damping is no longer performed, if it is enabled.
        It is recommended to leave damping on until convergence, which is the default. -*/
    options.add_double("DAMPING_CONVERGENCE", 1.0E-18);
    /*- The minimum iteration to start storing DIIS vectors -*/
    options.add_int("DIIS_START", 1);
    /*- Minimum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("DIIS_MIN_VECS", 2);
    /*- Maximum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("DIIS_MAX_VECS", 10);
    /*- Do use DIIS extrapolation to accelerate convergence? -*/
    options.add_bool("DIIS", true);
    /*- The iteration to start MOM on (or 0 for no MOM) -*/
    options.add_int("MOM_START", 0);
    /*- The absolute indices of orbitals to excite from in MOM (+/- for alpha/beta) -*/
    options.add("MOM_OCC", new ArrayType());
    /*- The absolute indices of orbitals to excite to in MOM (+/- for alpha/beta) -*/
    options.add("MOM_VIR", new ArrayType());
    /*- Whether to perform stability analysis after convergence.  NONE prevents analysis being
        performed. CHECK will print out the analysis of the wavefunction stability at the end of
        the computation.  FOLLOW will perform the analysis and, if a totally symmetric instability
        is found, will attemp to follow the eigenvector and re-run the computations to find a stable
        solution. -*/
    options.add_str("STABILITY_ANALYSIS", "NONE", "NONE CHECK FOLLOW");
    /*- When using STABILITY_ANALYSIS = FOLLOW, how much to scale the step along the eigenvector
        by. !expert -*/
    options.add_double("FOLLOW_STEP_SCALE", 0.5);

    /*- SUBSECTION Fractional Occupation UHF/UKS -*/

    /*- The iteration to start fractionally occupying orbitals (or 0 for no fractional occupation) -*/
    options.add_int("FRAC_START", 0);
    /*- The absolute indices of occupied orbitals to fractionally occupy (+/- for alpha/beta) -*/
    options.add("FRAC_OCC", new ArrayType());
    /*- The occupations of the orbital indices specified above ($0.0\ge occ \ge 1.0$) -*/
    options.add("FRAC_VAL", new ArrayType());
    /*- Do use DIIS extrapolation to accelerate convergence in frac? -*/
    options.add_bool("FRAC_DIIS", true);
    /*- Do renormalize C matrices prior to writing to checkpoint? -*/
    options.add_bool("FRAC_RENORMALIZE", true);
    /*- Do recompute guess from stored orbitals? -*/
    options.add_bool("FRAC_LOAD", false);

    /*- SUBSECTION Environmental Effects -*/

    /*- Do perturb the Hamiltonian? -*/
    options.add_bool("PERTURB_H", false);
    /*- Size of the perturbation (applies only to dipole perturbations) -*/
    options.add_double("PERTURB_MAGNITUDE", 0.0);
    /*- The operator used to perturb the Hamiltonian, if requested -*/
    options.add_str("PERTURB_WITH", "DIPOLE_X", "DIPOLE_X DIPOLE_Y DIPOLE_Z EMBPOT SPHERE");
    /*- An ExternalPotential (built by Python or NULL/None) -*/
    options.add("EXTERN", new PythonDataType());

    /*- Radius (bohr) of a hard-sphere external potential -*/
    options.add_double("RADIUS", 10.0); // bohr
    /*- Thickness (bohr) of a hard-sphere external potential -*/
    options.add_double("THICKNESS", 20.0); // bohr
    /*- Number of radial grid points for sphereical potential integration -*/
    options.add_int("R_POINTS", 100);
    /*- Number of colatitude grid points for sphereical potential integration -*/
    options.add_int("THETA_POINTS", 360);
    /*- Number of azimuthal grid points for sphereical potential integration -*/
    options.add_int("PHI_POINTS", 360);


    /*- SUBSECTION Parallel Runtime -*/

    /*- The dimension sizes of the processor grid !expert -*/
    options.add("PROCESS_GRID", new ArrayType());
    /*- The tile size for the distributed matrices !expert -*/
    options.add_int("TILE_SZ", 512);
    /*- The dimension sizes of the distributed matrix !expert -*/
    options.add("DISTRIBUTED_MATRIX", new ArrayType());
    /*- Do run in parallel? !expert -*/
    options.add_bool("PARALLEL", false);

    /*- SUBSECTION Misc. -*/

    /*- Are going to do SAPT? If so, what part? !expert -*/
    options.add_str("SAPT","FALSE","FALSE 2-DIMER 2-MONOMER_A 2-MONOMER_B 3-TRIMER 3-DIMER_AB 3-DIMER_BC 3-DIMER_AC 3-MONOMER_A 3-MONOMER_B 3-MONOMER_C");

    /*- SUBSECTION DFSCF Algorithm -*/

    /*- Number of threads for integrals (may be turned down if memory is an issue). 0 is blank -*/
    options.add_int("DF_INTS_NUM_THREADS",0);
    /*- IO caching for CP corrections, etc !expert -*/
    options.add_str("DF_INTS_IO", "NONE", "NONE SAVE LOAD");
    /*- Fitting Condition !expert -*/
    options.add_double("DF_FITTING_CONDITION", 1.0E-12);

    /*- SUBSECTION SAD Guess Algorithm -*/

    /*- The amount of SAD information to print to the output !expert -*/
    options.add_int("SAD_PRINT", 0);
    /*- Convergence criterion for SCF energy in SAD Guess. -*/
    options.add_double("SAD_E_CONVERGENCE", 1E-5);
    /*- Convergence criterion for SCF density in SAD Guess. -*/
    options.add_double("SAD_D_CONVERGENCE", 1E-5);
    /*- Maximum number of SAD guess iterations !expert -*/
    options.add_int("SAD_MAXITER", 50);
    /*- SAD Guess F-mix Iteration Start !expert -*/
    options.add_int("SAD_F_MIX_START", 50);
    /*- SAD Guess Cholesky Cutoff (for eliminating redundancies). !expert -*/
    options.add_double("SAD_CHOL_TOLERANCE", 1E-7);

    /*- SUBSECTION DFT -*/

    /*- The DFT combined functional name, e.g. B3LYP, or GEN to use a python reference to a 
        custom functional specified by DFT_CUSTOM_FUNCTIONAL. -*/
    options.add_str("DFT_FUNCTIONAL", "");
    /*- A custom DFT functional object (built by Python or NULL/None) -*/
    options.add("DFT_CUSTOM_FUNCTIONAL", new PythonDataType());
    /*- The DFT Range-separation parameter -*/
    options.add_double("DFT_OMEGA", 0.0);
    /*- The DFT Exact-exchange parameter -*/
    options.add_double("DFT_ALPHA", 0.0);
    /*- Number of spherical points (A :ref:`Lebedev Points <table:lebedevorder>` number). -*/
    options.add_int("DFT_SPHERICAL_POINTS", 302);
    /*- Number of radial points. -*/
    options.add_int("DFT_RADIAL_POINTS", 99);
    /*- Spherical Scheme. -*/
    options.add_str("DFT_SPHERICAL_SCHEME", "LEBEDEV", "LEBEDEV");
    /*- Radial Scheme. -*/
    options.add_str("DFT_RADIAL_SCHEME", "TREUTLER", "TREUTLER BECKE MULTIEXP EM MURA");
    /*- Nuclear Scheme. -*/
    options.add_str("DFT_NUCLEAR_SCHEME", "TREUTLER", "TREUTLER BECKE NAIVE STRATMANN");
    /*- Factor for effective BS radius in radial grid. -*/
    options.add_double("DFT_BS_RADIUS_ALPHA",1.0);
    /*- DFT basis cutoff. -*/
    options.add_double("DFT_BASIS_TOLERANCE", 1.0E-12);
    /*- The DFT grid specification, such as SG1.!expert -*/
    options.add_str("DFT_GRID_NAME","","SG1");
    /*- Pruning Scheme. !expert -*/
    options.add_str("DFT_PRUNING_SCHEME", "FLAT", "FLAT P_GAUSSIAN D_GAUSSIAN P_SLATER D_SLATER LOG_GAUSSIAN LOG_SLATER");
    /*- Spread alpha for logarithmic pruning. !expert -*/
    options.add_double("DFT_PRUNING_ALPHA",1.0);
    /*- The maximum number of grid points per evaluation block. !expert -*/
    options.add_int("DFT_BLOCK_MAX_POINTS",5000);
    /*- The minimum number of grid points per evaluation block. !expert -*/
    options.add_int("DFT_BLOCK_MIN_POINTS",1000);
    /*- The maximum radius to terminate subdivision of an octree block [au]. !expert -*/ 
    options.add_double("DFT_BLOCK_MAX_RADIUS",3.0);
    /*- The blocking scheme for DFT. !expert -*/
    options.add_str("DFT_BLOCK_SCHEME","OCTREE","NAIVE OCTREE");
  }
  if (name == "CPHF"|| options.read_globals()) {
    /*- The amount of information printed
        to the output file -*/
    options.add_int("PRINT", 1);
    /*- The amount of debug information printed
        to the output file -*/
    options.add_int("DEBUG", 0);
    /*- What app to test?
      -*/
    options.add_str("MODULE", "RCIS", "RCIS RCPHF RTDHF RCPKS RTDA RTDDFT");
    /*- Do singlet states? Default true
     -*/
    options.add_bool("DO_SINGLETS", true);
    /*- Do triplet states? Default true
     -*/
    options.add_bool("DO_TRIPLETS", true);
    /*- Do explicit hamiltonian only? -*/
    options.add_bool("EXPLICIT_HAMILTONIAN", false);
    /*- Minimum singles amplitude to print in 
        CIS analysis
     -*/
    options.add_double("CIS_AMPLITUDE_CUTOFF", 0.15);
    /*- Memory safety factor for allocating JK
    -*/
    options.add_double("TDHF_MEM_SAFETY_FACTOR",0.75);
    /*- Memory safety factor for allocating JK
    -*/
    options.add_double("CIS_MEM_SAFETY_FACTOR",0.75);
    /*- Memory safety factor for allocating JK
    -*/
    options.add_double("CPHF_MEM_SAFETY_FACTOR",0.75);
    /*- Which states to save AO OPDMs for?
     *   Positive - Singlets
     *   Negative - Triplets
     * -*/
    options.add("CIS_OPDM_STATES", new ArrayType());
    /*- Which states to save AO transition OPDMs for?
     *   Positive - Singlets
     *   Negative - Triplets
     * -*/
    options.add("CIS_TOPDM_STATES", new ArrayType());
    /*- Which states to save AO difference OPDMs for?
     *   Positive - Singlets
     *   Negative - Triplets
     * -*/
    options.add("CIS_DOPDM_STATES", new ArrayType());
    /*- Which states to save AO Natural Orbitals for?
     *   Positive - Singlets
     *   Negative - Triplets
     * -*/
    options.add("CIS_NO_STATES", new ArrayType());
    /*- Which states to save AD Matrices for?
     *   Positive - Singlets
     *   Negative - Triplets
     * -*/
    options.add("CIS_AD_STATES", new ArrayType());
    /*- Which tasks to run CPHF For
     *  Valid choices:
     *  -Polarizability
     * -*/
    options.add("CPHF_TASKS", new ArrayType());
    /*- The maximum number of integral threads (0 for omp_get_max_threads()) 
     -*/
    options.add_int("OMP_N_THREAD", 0);
    /*- The schwarz cutoff value 
     -*/
    options.add_double("SCHWARZ_CUTOFF", 1.0E-12);
    /*- The maximum reciprocal condition allowed in the fitting metric 
     -*/
    options.add_double("FITTING_CONDITION", 1.0E-12);
    /*- Fitting algorithm (0 for old, 1 for new)
     -*/
    options.add_int("FITTING_ALGORITHM", 0);
    /*- SCF Type 
     -*/
    options.add_str("SCF_TYPE", "DIRECT", "DIRECT DF PK OUT_OF_CORE PS");
    /*- Auxiliary basis for SCF 
     -*/
    options.add_str("DF_BASIS_SCF", ""); 
    /*- Solver maximum iterations
     -*/
    options.add_int("SOLVER_MAXITER",100);
    /*- Solver convergence threshold (max 2-norm). -*/
    options.add_double("SOLVER_CONVERGENCE",1.0E-6);
    /*- DL Solver number of roots 
     -*/
    options.add_int("SOLVER_N_ROOT",1);
    /*- DL Solver number of guesses
     -*/
    options.add_int("SOLVER_N_GUESS",1);
    /*- DL Solver number of subspace vectors to collapse to
     -*/
    options.add_int("SOLVER_MIN_SUBSPACE",2);
    /*- DL Solver maximum number of subspace vectors 
     -*/
    options.add_int("SOLVER_MAX_SUBSPACE",6);
    /*- DL Solver minimum corrector norm to add to subspace
     -*/
    options.add_double("SOLVER_NORM",1.0E-6);
    /*- Solver precondition type
     -*/
    options.add_str("SOLVER_PRECONDITION","JACOBI","SUBSPACE JACOBI NONE");
    /*- Solver type (for interchangeable solvers)
     -*/
    options.add_str("SOLVER_TYPE", "DL", "DL RAYLEIGH"); 
    /*- Solver precondtion max steps
    -*/
    options.add_int("SOLVER_PRECONDITION_MAXITER", 1);
    /*- Solver precondition step type
    -*/
    options.add_str("SOLVER_PRECONDITION_STEPS", "TRIANGULAR", "CONSTANT TRIANGULAR");   
    /*- Solver residue or eigenvector delta
    -*/
    options.add_str("SOLVER_QUANTITY", "RESIDUAL", "EIGENVECTOR RESIDUAL");
    /*- Solver exact diagonal or eigenvalue difference?
    -*/
    options.add_bool("SOLVER_EXACT_DIAGONAL", false);

  }
  if (name == "MP2"|| options.read_globals()) {
      /*- MODULEDESCRIPTION Performs second order Moller-Plesset perturbation theory (MP2) computations.  This code can
          compute RHF/ROHF/UHF energies, and RHF gradient/property computations.  However, given the small errors introduced,
          we recommend using the new density fitted MP2 codes instead, which are much more efficient. -*/
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "MP2", "MP2");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF", "RHF UHF ROHF");
    /*- Type of job being performed !expert -*/
    options.add_str("JOBTYPE", "SP");
    /*- Do compute the one particle density matrix, for properties? -*/
    options.add_bool("OPDM", false);
    /*- Do add relaxation terms to the one particle density matrix, for properties? -*/
    options.add_bool("OPDM_RELAX", false);
    /*- The amount of cacheing of data to perform -*/
    options.add_int("CACHELEVEL", 2);
    /*- The criterion used to retain/release cached data -*/
    options.add_str("CACHETYPE", "LRU", "LRU LOW");
    /*- Do perform a spin component scaled MP2 computation? -*/
    options.add_bool("SCS", false);
    /*- Do perform a spin component scaled (N) MP2 computation? -*/
    options.add_bool("SCS_N", false);
    /*- The scale factor used for opposite-spin pairs in SCS computations -*/
    options.add_double("MP2_OS_SCALE", 6.0/5.0);
    /*- The scale factor used for same-spin pairs in SCS computations-*/
    options.add_double("MP2_SS_SCALE", 1.0/3.0);
  }
  // Options of this module not standardized since it's bound for deletion
  if(name == "TRANSQT2"|| options.read_globals()) {
      /*- MODULEDESCRIPTION Performs transformations of integrals into the molecular orbital (MO) basis.  This
          module is currently used by the (non-density fitted) MP2 and coupled cluster codes, but it is being phased
          out. -*/
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
    /*- Do print two-electron integrals (TEIs)? -*/
    options.add_bool("PRINT_TEI", false);
    /*- Minimum absolute value below which integrals are neglected. -*/
    options.add_double("INTS_TOLERANCE", 1e-14);
    /*- Controls how to cache quantities within the DPD library !expert-*/
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
      /*- MODULEDESCRIPTION The predecessor to Transqt2.  Currently used by the configuration interaction codes, but
          it is being phased out. -*/
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "CCSD");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
    /*- The amount of information to print to the output file.  1 prints
     basic information, and higher levels print more information. A value
     of 5 will print very large amounts of debugging information. -*/
    options.add_int("PRINT_LVL", 1);
    /*- The way of transformation, from ao basis to mo basis or vice versa -*/
    options.add_str("MODE", "TO_MO", "TO_MO TO_AO");
    /*- Do specific arrangements for PSIMRCC? -*/
    options.add_bool("PSIMRCC", false);
    /*- Transformations for explicitly-correlated MP2 methods -*/
    options.add_str("MP2R12A", "MP2R12AERI", "MP2R12AERI MP2R12AR12 MP2R12AR12T1");
    /*- Minimum absolute value below which integrals are neglected. -*/
    options.add_double("INTS_TOLERANCE", 1e-14);
    /*- One-electron parameters file -*/
    options.add_int("OEI_FILE", PSIF_OEI);
    /*- Alpha-spin one-electron parameters file -*/
    options.add_int("OEI_A_FILE", PSIF_OEI);
    /*- Beta-spin one-electron parameters file -*/
    options.add_int("OEI_B_FILE", PSIF_OEI);
    /*- Frozen-core file -*/
    options.add_int("FZC_FILE", PSIF_OEI);
    /*- Alpha-spin frozen-core file -*/
    options.add_int("FZC_A_FILE", PSIF_OEI);
    /*- Beta-spin frozen-core file -*/
    options.add_int("FZC_B_FILE", PSIF_OEI);
    /*- MO-basis sorted two-electron integrals file -*/
    options.add_int("SORTED_TEI_FILE", PSIF_MO_TEI);
    /*- MO-basis two-particle density matrix file -*/
    options.add_int("TPDM_FILE", PSIF_MO_TPDM);
    /*- SO basis overlap matrix file -*/
    options.add_int("SO_S_FILE", PSIF_OEI);
    /*- SO basis kinetic energy matrix file -*/
    options.add_int("SO_T_FILE", PSIF_OEI);
    /*- SO basis potential energy matrix file -*/
    options.add_int("SO_V_FILE", PSIF_OEI);
    /*- SO basis two-electron integrals file -*/
    options.add_int("SO_TEI_FILE", PSIF_SO_TEI); // ?
    /*- First temporary file -*/
    options.add_int("FIRST_TMP_FILE", 150);
    /*- MO-basis one-particle density matrix file -*/
    options.add_int("OPDM_IN_FILE", PSIF_MO_OPDM);
    /*- AO-basis one-particle density matrix file -*/
    options.add_int("OPDM_OUT_FILE", PSIF_AO_OPDM);
    /*- MO-basis MO-lagrangian file -*/
    options.add_int("LAG_IN_FILE", PSIF_MO_LAG);
    /*- SO-basis presort file -*/
    options.add_int("PRESORT_FILE", PSIF_SO_PRESORT);
    /*- Do keep presort file? -*/
    options.add_bool("KEEP_PRESORT", false);
    /*- Half-transformed integrals -*/
    options.add_int("J_FILE", 91);
    /*- Do keep half-transformed integrals? -*/
    options.add_bool("KEEP_J", false);
    /*- Output integrals file -*/
    options.add_int("M_FILE", 0); // output integrals file; depends on direction
    /*- MO basis (PQ|RS) type two-electron integrals file -*/
    options.add_int("AA_M_FILE", PSIF_MO_AA_TEI);
    /*- MO basis (pq|rs) type two-electron integrals file -*/
    options.add_int("BB_M_FILE", PSIF_MO_BB_TEI);
    /*- MO basis (PQ|rs) type two-electron integrals file -*/
    options.add_int("AB_M_FILE", PSIF_MO_AB_TEI);
    /*- Maximum buckets -*/
    options.add_int("MAX_BUCKETS", 499);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- Do delete AO integral files? -*/
    options.add_bool("DELETE_AO", true);
    /*- Do delete TPDM file? -*/
    options.add_bool("DELETE_TPDM", true);
    /*- Do print two-electron integrals? -*/
    options.add_bool("PRINT_TE_INTEGRALS", false);
    /*- Do print one-electron integrals? -*/
    options.add_bool("PRINT_OE_INTEGRALS", false);
    /*- Do print sorted one-electron integrals? -*/
    options.add_bool("PRINT_SORTED_OE_INTS", false);
    /*- Do print sorted two-electron integrals (TEIs)? -*/
    options.add_bool("PRINT_SORTED_TE_INTS", false);
    /*- Do print MOs? -*/
    options.add_bool("PRINT_MOS", false);

    /*- Do multiply the MO-lagrangian by 2.0? -*/
    options.add_bool("LAGRAN_DOUBLE", false);
    /*- Do divide the MO-lagrangian by 2.0? -*/
    options.add_bool("LAGRAN_HALVE", false);
    /*- Do transform all TEIs -*/
    options.add_bool("DO_ALL_TEI", false);
    /*- Do add reference contribution to TPDM? -*/
    options.add_bool("TPDM_ADD_REF", false);
    /*- Do delete restricted doubly occupieds? -*/
    options.add_bool("DELETE_RESTR_DOCC", true);
    /*- Do print reordered MOs? -*/
    options.add_bool("PRINT_REORDER", false);
    /*- Do use Pitzer ordering? -*/
    options.add_bool("PITZER", false);
    /*- Do reorder MOs? -*/
    options.add_bool("REORDER", false);
    /*- Do check MO orthogonality condition? -*/
    options.add_bool("CHECK_C_ORTHONORM", false);
    /*- Do form quasi RHF (QRHF) orbitals? -*/
    options.add_bool("QRHF", false);
    /*- Do form improved virtual orbitals (IVO)? -*/
    options.add_bool("IVO", false);
    /*- Numbering of MOs for reordering requests?  -*/
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
      /*- MODULEDESCRIPTION Sorts integrals for efficiency. Called before (non density-fitted) MP2 and
          coupled cluster computations. -*/
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF");
    /*- Reference wavefunction type for EOM computations -*/
    options.add_str("EOM_REFERENCE","RHF");
    /*- The response property desired.  The unique acceptable values is ``POLARIZABILITY``
    for dipole-polarizabilitie. -*/
    options.add_str("PROPERTY", "POLARIZABILITY");
    /*- Do simulate the effects of local correlation techniques? -*/
    options.add_bool("LOCAL", false);
    /*- Value (always between one and zero) for the Broughton-Pulay completeness
    check used to contruct orbital domains for local-CC calculations. See
    J. Broughton and P. Pulay, J. Comp. Chem. 14, 736-740 (1993) and C. Hampel
    and H.-J. Werner, J. Chem. Phys. 104, 6286-6297 (1996). -*/
    options.add_double("LOCAL_CUTOFF", 0.02);
    /*- Cutoff value for local-coupled-perturbed-Hartree-Fock -*/
    options.add_double("LOCAL_CPHF_CUTOFF", 0.10);
    /*- Local core cutoff value -*/
    options.add_double("LOCAL_CORE_CUTOFF",0.05);
    /*- Type of local-CCSD scheme to be simulated. ``WERNER`` (unique avaliable option) selects the method
    developed by H.-J. Werner and co-workers. -*/
    options.add_str("LOCAL_METHOD","WERNER");
    /*- Desired treatment of "weak pairs" in the local-CCSD method. The value of ``NONE`` (unique avaliable option) treats weak pairs in
    the same manner as strong pairs. -*/
    options.add_str("LOCAL_WEAKP","NONE");
    /*- Definition of local pair domains, unique avaliable option is BP, Boughton-Pulay. -*/
    options.add_str("LOCAL_PAIRDEF","BP");
    /*- Do use augment domains with polarized orbitals? -*/
    options.add_bool("LOCAL_DOMAIN_POLAR", false);
    /*- Do generate magnetic-field CPHF solutions for local-CC? -*/
    options.add_bool("LOCAL_DOMAIN_MAG", false);
    /*- -*/
    options.add_bool("LOCAL_DOMAIN_SEP", false);
    /*- Do apply local filtering to single excitation amplitudes? -*/
    options.add_bool("LOCAL_FILTER_SINGLES", false);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- Do retain the input two-electron integrals? -*/
    options.add_bool("KEEP_TEIFILE", false);
    /*- Do retain the input one-electron integrals? -*/
    options.add_bool("KEEP_OEIFILE", false);
    /*- Minimum absolute value below which integrals are neglected. -*/
    options.add_double("INTS_TOLERANCE", 1e-14);
    /*- Cacheing level for libdpd governing the storage of amplitudes,
    integrals, and intermediates in the CC procedure. A value of 0 retains
    no quantities in cache, while a level of 6 attempts to store all
    quantities in cache.  For particularly large calculations, a value of
    0 may help with certain types of memory problems.  The default is 2,
    which means that all four-index quantites with up to two virtual-orbital
    indices (e.g., $\langle ij | ab \rangle>$ integrals) may be held in the cache. -*/
    options.add_int("CACHELEVEL", 2);
    /*- Energy of applied field [au] for dynamic properties -*/
    options.add("OMEGA", new ArrayType());
    /*- Convert ROHF MOs to semicanonical MOs -*/
    options.add_bool("SEMICANONICAL", true);
  }
  if(name == "CCTRIPLES"|| options.read_globals()) {
     /*- MODULEDESCRIPTION Computes the triples component of CCSD(T) energies (and gradients, if necessary). -*/
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "SCF");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
    /*- Number of threads -*/
    options.add_int("CC_NUM_THREADS",1);
    /*- Convert ROHF MOs to semicanonical MOs -*/
    options.add_bool("SEMICANONICAL", true);
  }
  if(name == "CCDENSITY"|| options.read_globals()) {
     /*- MODULEDESCRIPTION Computes the coupled cluster density matrices. Called whenever CC properties and/or
         gradients are required. -*/
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "SCF");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
    /*- Minimum absolute value below which integrals are neglected. -*/
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
    options.add_bool("OPDM_RELAX",false);
    /*- Do require $\bar{H}$ and $R$ to be connected? !expert -*/
    options.add_bool("XI_CONNECT",false);
    /*- The number of electronic states to computed, per irreducible
    representation -*/
    options.add("ROOTS_PER_IRREP", new ArrayType());
    /*- Compute non-relaxed properties for all excited states. -*/
    options.add_bool("PROP_ALL",true);
    /*- The symmetry of states -*/
    options.add_int("PROP_SYM", 1);
    /*- Root number (within its irrep) for computing properties -*/
    options.add_int("PROP_ROOT", 1);
    /*- Do compute Xi? -*/
    options.add_bool("XI", false);
    /*- Do use zeta?  -*/
    options.add_bool("ZETA",false);
    /*- Do compute one-particle density matrix? -*/
    options.add_bool("ONEPDM",false);
    /*- Write one-particle density matrix on a grid to file opdm.dx -*/
    options.add_bool("ONEPDM_GRID_DUMP",false);
    /*- Cutoff (e/A^3) for printing one-particle density matrix values on a grid -*/
    options.add_double("ONEPDM_GRID_CUTOFF", 1.0e-30);
    /*- Stepsize (Angstrom) for one-particle density matrix values on a grid -*/
    options.add_double("ONEPDM_GRID_STEPSIZE", 0.1);
  }
  if(name == "CCLAMBDA"|| options.read_globals()) {
     /*- MODULEDESCRIPTION Solves for the Lagrange multipliers, which are needed whenever coupled cluster properties
         or gradients are requested. -*/
    /*- Wavefunction type !expert -*/
    options.add_str("WFN","SCF");
    /*- Convergence criterion for wavefunction (change) in CC lambda-amplitude equations. -*/
    options.add_double("R_CONVERGENCE",1e-7);
    /*- Do restart the coupled-cluster iterations from old $\lambda@@1$ and $\lambda@@2$
    amplitudes? -*/
    options.add_bool("RESTART",false);
    /*- Cacheing level for libdpd governing the storage of amplitudes,
    integrals, and intermediates in the CC procedure. A value of 0 retains
    no quantities in cache, while a level of 6 attempts to store all
    quantities in cache.  For particularly large calculations, a value of
    0 may help with certain types of memory problems.  The default is 2,
    which means that all four-index quantites with up to two virtual-orbital
    indices (e.g., $\langle ij | ab \rangle>$ integrals) may be held in the cache. -*/
    options.add_int("CACHELEVEL",2);
    /*- Do Sekino-Bartlett size-extensive model-III? -*/
    options.add_bool("SEKINO",false);
    /*- Do use DIIS extrapolation to accelerate convergence? -*/
    options.add_bool("DIIS",true);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- Type of ABCD algorithm will be used -*/
    options.add_str("ABCD","NEW");
    /*- Number of important CC amplitudes per excitation level to print.
    CC analog to |detci__num_dets_print|. -*/
    options.add_int("NUM_AMPS_PRINT",10);
    /*- Type of job being performed !expert -*/
    options.add_str("JOBTYPE","");
    /*- Do simulate the effects of local correlation techniques? -*/
    options.add_bool("LOCAL",false);
    /*- Desired treatment of "weak pairs" in the local-CCSD method. The value of ``NONE`` (unique avaliable option) treats weak pairs in
    the same manner as strong pairs. -*/
    options.add_str("LOCAL_WEAKP","NONE");
    /*- Value (always between one and zero) for the Broughton-Pulay completeness
    check used to contruct orbital domains for local-CC calculations. See
    J. Broughton and P. Pulay, J. Comp. Chem. 14, 736-740 (1993) and C. Hampel
    and H.-J. Werner, J. Chem. Phys. 104, 6286-6297 (1996). -*/
    options.add_double("LOCAL_CUTOFF",0.02);
    /*- Type of local-CCSD scheme to be simulated. ``WERNER`` (unique avaliable option) selects the method
    developed by H.-J. Werner and co-workers. -*/
    options.add_str("LOCAL_METHOD","WERNER");
    /*- Do apply local filtering to single de-excitation ($\lambda 1$ amplitudes? -*/
    options.add_bool("LOCAL_FILTER_SINGLES",true);
    /*- Cutoff value for local-coupled-perturbed-Hartree-Fock -*/
    options.add_double("LOCAL_CPHF_CUTOFF",0.10);
    /*- Definition of local pair domains -*/
    options.add_str("LOCAL_PAIRDEF","");
    /*- The number of electronic states to computed, per irreducible
    representation -*/
    options.add("ROOTS_PER_IRREP", new ArrayType());
    /*- Compute unrelaxed properties for all excited states. -*/
    options.add_bool("PROP_ALL",true);
    /*- The symmetry of states -*/
    options.add_int("PROP_SYM",1);
    /*- Root number (within its irrep) for computing properties -*/
    options.add_int("PROP_ROOT",1);
    /*- Maximum number of iterations -*/
    options.add_int("MAXITER",50);
    /*- Do use zeta?  -*/
    options.add_bool("ZETA",false);
  }
  if(name == "CLAG"|| options.read_globals()) {
     /*- MODULEDESCRIPTION Solves for the CI Lagrangian. Called whenever CI properties or gradients are requested. -*/
    /*- Wavefunction type !expert -*/
    options.add_str("WFN","NONE");
    /*- Do write the OEI, TEI, OPDM, TPDM, and Lagrangian files in canonical form, Pitzer order? -*/
    options.add_bool("CAS_FILES_WRITE",0);
    /*- Root to get OPDM -*/
    options.add_int("FOLLOW_ROOT",1);
  }
  if(name == "STABILITY"|| options.read_globals()) {
     /*- MODULEDESCRIPTION Performs wavefunction stability analysis. Called when specifically requested
         by the user. -*/
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF", "RHF UHF ROHF");
    /*- -*/
    options.add_int("CACHELEVEL",2);
    /*- Do follow the most negative eigenvalue of the Hessian towards a lower
    energy HF solution? Follow a UHF $\rightarrow$ UHF instability of same symmetry? -*/
    options.add_bool("FOLLOW",false);
    /*- Number of lowest MO Hessian eigenvalues to print -*/
    options.add_int("NUM_VECS_PRINT",0);
    /*- Method for following eigenvectors, either 0 by angles or 1 by antisymmetric matrix. -*/
    options.add_int("ROTATION_SCHEME",0);
    /*- Scale factor (between 0 and 1) for orbital rotation step -*/
    options.add_double("SCALE",0.5);
  }
  if(name == "ADC" || options.read_globals()) {
     /*- MODULEDESCRIPTION Performs Algebraic-Diagrammatic Construction (ADC) propagator computations for excited states. -*/
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "");
    /*- How to cache quantities within the DPD library -*/
    options.add_int("CACHELEVEL", 2);
    /*- The amount of memory available (in Mb) -*/
    options.add_int("MEMORY", 1000);
    /*- The convergence criterion for pole searching step. -*/
    options.add_double("NEWTON_CONVERGENCE", 1e-7);
    /*- Maximum iteration number in pole searching -*/
    options.add_int("POLE_MAXITER", 20);
    /*- Maximum iteration number in simultaneous expansion method -*/
    options.add_int("SEM_MAXITER", 30);
    /*- The cutoff norm of residual vector in SEM step. -*/
    options.add_double("NORM_TOLERANCE", 1e-6);
    /*- The poles per irrep vector -*/
    options.add("ROOTS_PER_IRREP", new ArrayType());
    /*- Do use the partial renormalization scheme for the ground state wavefunction? -*/
    options.add_bool("PR", false);
    /*- Number of components of transition amplitudes printed -*/
    options.add_int("NUM_AMPS_PRINT", 5);
  }
  if(name == "CCHBAR"|| options.read_globals()) {
     /*- MODULEDESCRIPTION Assembles the coupled cluster effective Hamiltonian. Called whenever CC
         properties and/or gradients are required. -*/
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "SCF");
    /*- Reference wavefunction type for EOM computations -*/
    options.add_str("EOM_REFERENCE","RHF");
    /*- Do compute the Tamplitude equation matrix elements? -*/
    options.add_bool("T_AMPS",false);
    /*- Cacheing level for libdpd governing the storage of amplitudes,
    integrals, and intermediates in the CC procedure. A value of 0 retains
    no quantities in cache, while a level of 6 attempts to store all
    quantities in cache.  For particularly large calculations, a value of
    0 may help with certain types of memory problems.  The default is 2,
    which means that all four-index quantites with up to two virtual-orbital
    indices (e.g., $\langle ij | ab \rangle>$ integrals) may be held in the cache. -*/
    options.add_int("CACHELEVEL",2);
    /*- Do use the minimal-disk algorithm for Wabei? It's VERY slow! -*/
    options.add_bool("WABEI_LOWDISK", false);
  }
  if(name == "CCEOM"|| options.read_globals()) {
     /*- MODULEDESCRIPTION Performs equation-of-motion (EOM) coupled cluster excited state computations. -*/
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "EOM_CCSD", "EOM_CCSD EOM_CC2 EOM_CC3");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF", "RHF ROHF UHF");
    /*- Reference wavefunction type for EOM computations -*/
    options.add_str("EOM_REFERENCE","RHF", "RHF ROHF UHF");
    /*- Do use full effective Hamiltonian matrix? -*/
    options.add_bool("FULL_MATRIX",false);
    /*- Cacheing level for libdpd governing the storage of amplitudes,
    integrals, and intermediates in the CC procedure. A value of 0 retains
    no quantities in cache, while a level of 6 attempts to store all
    quantities in cache.  For particularly large calculations, a value of
    0 may help with certain types of memory problems.  The default is 2,
    which means that all four-index quantites with up to two virtual-orbital
    indices (e.g., $\langle ij | ab \rangle>$ integrals) may be held in the cache. -*/
    options.add_int("CACHELEVEL",2);
    /*- The criterion used to retain/release cached data -*/
    options.add_str("CACHETYPE", "LRU", "LOW LRU");
    /*- Number of threads -*/
    options.add_int("CC_NUM_THREADS", 1);
    /*- Type of ABCD algorithm will be used -*/
    options.add_str("ABCD", "NEW", "NEW OLD");
    /*- Do build W intermediates required for eom_cc3 in core memory? -*/
    options.add_bool("T3_WS_INCORE", false);
    /*- Do simulate the effects of local correlation techniques? -*/
    options.add_bool("LOCAL", false);
     /*- Value (always between one and zero) for the Broughton-Pulay completeness
    check used to contruct orbital domains for local-CC calculations. See
    J. Broughton and P. Pulay, J. Comp. Chem. 14, 736-740 (1993) and C. Hampel
    and H.-J. Werner, J. Chem. Phys. 104, 6286-6297 (1996). -*/
    options.add_double("LOCAL_CUTOFF", 0.02);
    /*- Type of local-CCSD scheme to be simulated. ``WERNER`` selects the method
    developed by H.-J. Werner and co-workers, and ``AOBASIS`` selects the method
    developed by G.E. Scuseria and co-workers (currently inoperative). -*/
    options.add_str("LOCAL_METHOD", "WERNER", "WERNER AOBASIS");
    /*- Desired treatment of "weak pairs" in the local-CCSD method. A value of
    ``NEGLECT`` ignores weak pairs entirely. A value of ``NONE`` treats weak pairs in
    the same manner as strong pairs. A value of MP2 uses second-order perturbation
    theory to correct the local-CCSD energy computed with weak pairs ignored. -*/
    options.add_str("LOCAL_WEAKP", "NONE", "NONE MP2 NEGLECT");
    /*- Preconditioner will be used in local CC computations -*/
    options.add_str("LOCAL_PRECONDITIONER", "HBAR", "HBAR FOCK");
    /*- Permit ghost atoms to hold projected atomic orbitals to include in the virtual space in local-EOM-CCSD calculations -*/
    options.add_int("LOCAL_GHOST", -1);
    /*- -*/
    options.add_bool("LOCAL_DO_SINGLES", true);
    /*- Do apply local filtering to singles amplitudes? -*/
    options.add_bool("LOCAL_FILTER_SINGLES", true);
    /*- Do use new triples? -*/
    options.add_bool("NEW_TRIPLES", true);
    /*- Number of excited states per irreducible representation for EOM-CC
    and CC-LR calculations. Irreps denote the final state symmetry, not the
    symmetry of the transition. -*/
    options.add("ROOTS_PER_IRREP", new ArrayType());
    /*- Maximum number of iterations -*/
    options.add_int("MAXITER", 80);
    /*- Symmetry of the state to compute properties. Defaults to last irrep
    for which states are requested. -*/
    options.add_int("PROP_SYM", 1);
    /*- Root number (within its irrep) for computing properties. Defaults to
    highest root requested. -*/
    options.add_int("PROP_ROOT", 0);
    /*- Do turn on root following for CC3 -*/
    options.add_bool("CC3_FOLLOW_ROOT", false);
    /*- Do form a triplet state from RHF reference? -*/
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
    /*- SS vectors stored per root -*/
    options.add_int("SS_VECS_PER_ROOT", 5);
    /*- Vectors stored per root -*/
    options.add_int("VECS_PER_ROOT", 12);
    /*- Vectors stored in CC3 computations -*/
    options.add_int("VECS_CC3", 10);
    /*- Do collapse with last vector? -*/
    options.add_bool("COLLAPSE_WITH_LAST", true);
    /*- Complex tolerance applied in CCEOM computations -*/
    options.add_double("COMPLEX_TOLERANCE", 1E-12);
    /*- Convergence criterion for norm of the residual vector in the Davidson algorithm for CC-EOM. -*/
    options.add_double("R_CONVERGENCE", 1E-6);
    /*- Convergence criterion for norm of the residual vector in the Davidson algorithm for the CIS guess to CC-EOM. -*/
    options.add_double("SS_R_CONVERGENCE", 1E-6);
    /*- Convergence criterion for excitation energy (change) in the Davidson algorithm for CC-EOM. -*/
    options.add_double("E_CONVERGENCE", 1E-8);
    /*- Convergence criterion for excitation energy (change) in the Davidson algorithm for the CIS guess to CC-EOM. -*/
    options.add_double("SS_E_CONVERGENCE", 1E-6);
    /*- Number of important CC amplitudes to print -*/
    options.add_int("NUM_AMPS_PRINT", 5);
    /*- Minimum absolute value above which a guess vector to a root is added
    to the Davidson algorithm in the EOM-CC iterative procedure. -*/
    options.add_double("SCHMIDT_ADD_RESIDUAL_TOLERANCE", 1E-3);
    /*- Do skip diagonalization of Hbar SS block? -*/
    options.add_bool("SS_SKIP_DIAG", false);
    /*- Do restart from on-disk? -*/
    options.add_bool("RESTART_EOM_CC3", false);
    /*- Specifies a set of single-excitation guess vectors for the EOM-CC
    procedure.  If EOM_GUESS = ``SINGLES``, the guess will be taken from
    the singles-singles block of the similarity-transformed Hamiltonian,
    Hbar.  If EOM_GUESS = ``DISK``, guess vectors from a previous computation
    will be read from disk.  If EOM_GUESS = ``INPUT``, guess vectors will be
    specified in user input.  The latter method is not currently available. -*/
    options.add_str("EOM_GUESS", "SINGLES", "SINGLES DISK INPUT");
    /*- Convert ROHF MOs to semicanonical MOs -*/
    options.add_bool("SEMICANONICAL", true);
  }
  if(name == "CCRESPONSE"|| options.read_globals()) {
     /*- MODULEDESCRIPTION Performs coupled cluster response property computations. -*/
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "SCF");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF");
    /*- Cacheing level for libdpd -*/
    options.add_int("CACHELEVEL",2);
    /*- Specifies the choice of representation of the electric dipole operator.
    Acceptable values are ``LENGTH`` for the usual length-gauge representation,
    ``VELOCITY`` for the modified velocity-gauge representation in which the
    static-limit optical rotation tensor is subtracted from the frequency-
    dependent tensor, or ``BOTH``. Note that, for optical rotation calculations,
    only the choices of ``VELOCITY`` or ``BOTH`` will yield origin-independent results. -*/
    options.add_str("GAUGE","LENGTH", "LENGTH VELOCITY BOTH");
    /*- Maximum number of iterations to converge perturbed amplitude equations -*/
    options.add_int("MAXITER",50);
    /*- Convergence criterion for wavefunction (change) in perturbed CC equations. -*/
    options.add_double("R_CONVERGENCE",1e-7);
    /*- Do use DIIS extrapolation to accelerate convergence? -*/
    options.add_bool("DIIS",1);
    /*- The response property desired.  Acceptable values are ``POLARIZABILITY``
    (default) for dipole-polarizabilities, ``ROTATION`` for specific rotations,
    ``ROA`` for Raman Optical Activity, and ``ALL`` for all of the above. -*/
    options.add_str("PROPERTY","POLARIZABILITY","POLARIZABILITY ROTATION ROA ALL");
    /*- Type of ABCD algorithm will be used -*/
    options.add_str("ABCD","NEW");
    /*- Do restart from on-disk amplitudes? -*/
    options.add_bool("RESTART",1);
    /*- Do simulate local correlation? -*/
    options.add_bool("LOCAL",0);
    /*- Value (always between one and zero) for the Broughton-Pulay completeness
    check used to contruct orbital domains for local-CC calculations. See
    J. Broughton and P. Pulay, J. Comp. Chem. 14, 736-740 (1993) and C. Hampel
    and H.-J. Werner, J. Chem. Phys. 104, 6286-6297 (1996). -*/
    options.add_double("LOCAL_CUTOFF",0.01);
    /*- Type of local-CCSD scheme to be simulated. ``WERNER`` (unique avaliable option) selects the method
    developed by H.-J. Werner and co-workers. -*/
    options.add_str("LOCAL_METHOD","WERNER");
    /*- Desired treatment of "weak pairs" in the local-CCSD method. The value of ``NONE`` (unique avaliable option) treats weak pairs in
    the same manner as strong pairs. -*/
    options.add_str("LOCAL_WEAKP","NONE");
    /*- Do apply local filtering to single excitation amplitudes? -*/
    options.add_bool("LOCAL_FILTER_SINGLES", false);
    /*- Cutoff value for local-coupled-perturbed-Hartree-Fock -*/
    options.add_double("LOCAL_CPHF_CUTOFF",0.10);
    /*- Definition of local pair domains -*/
    options.add_str("LOCAL_PAIRDEF","NONE");
    /*- Do analyze X2 amplitudes -*/
    options.add_bool("ANALYZE",0);
    /*- Number of important CC amplitudes per excitation level to print.
    CC analog to |detci__num_dets_print|. -*/
    options.add_int("NUM_AMPS_PRINT",5);
    /*- Do Sekino-Bartlett size-extensive model-III? -*/
    options.add_bool("SEKINO",0);
    /*- Do Bartlett size-extensive linear model? -*/
    options.add_bool("LINEAR",0);
    /*- Array that specifies the desired frequencies of the incident
    radiation field in CCLR calculations.  If only one element is
    given, the units will be assumed to be atomic units.  If more
    than one element is given, then the units must be specified as the final
    element of the array.  Acceptable units are ``HZ``, ``NM``, ``EV``, and ``AU``. -*/
    options.add("OMEGA",new ArrayType());
  }
  if(name == "RESPONSE"|| options.read_globals()){
     /*- MODULEDESCRIPTION Performs SCF linear response computations. -*/
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF");
    /*- Array that specifies the desired frequencies of the incident
    radiation field in CCLR calculations.  If only one element is
    given, the units will be assumed to be atomic units.  If more
    than one element is given, then the units must be specified as the final
    element of the array.  Acceptable units are ``HZ``, ``NM``, ``EV``, and ``AU``. -*/
    options.add("OMEGA", new ArrayType());
    /*- Array that specifies the desired frequencies of the incident
    radiation field in CCLR calculations.  If only one element is
    given, the units will be assumed to be atomic units.  If more
    than one element is given, then the units must be specified as the final
    element of the array.  Acceptable units are HZ, NM, EV, and AU. -*/
    /*- The response property desired.  Acceptable values are POLARIZABILITY
    (default) for dipole-polarizabilities, ROTATION for specific rotations,
    ROA for Raman Optical Activity, and ALL for all of the above.
    -*/
    options.add_str("PROPERTY","POLARIZABILITY","POLARIZABILITY ROTATION ROA ALL");
  }
  if(name == "MCSCF"|| options.read_globals()) {
     /*- MODULEDESCRIPTION Performs RHF/UHF/ROHF/TCSCF and more general MCSCF computations. Called
         as the starting point for multireference coupled cluster computations. -*/
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE","RHF","RHF ROHF UHF TWOCON MCSCF GENERAL");
    /*- Level shift to aid convergence -*/
    options.add_double("LEVEL_SHIFT",0.0);
    /*- Convergence criterion for energy. -*/
    options.add_double("E_CONVERGENCE", 1e-8);
    /*- Convergence criterion for density. -*/
    options.add_double("D_CONVERGENCE", 1e-6);
    /*- Maximum number of iterations -*/
    options.add_int("MAXITER",100);
    /*- Maximum number of error vectors stored for DIIS extrapolation -*/
    options.add_int("DIIS_MAX_VECS",7);
    /*- Which solution of the SCF equations to find, where 1 is the SCF ground state-*/
    options.add_int("FOLLOW_ROOT",1);
    /*- Iteration at which to begin using the averaged Fock matrix-*/
    options.add_int("FAVG_START",5);
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
    options.add_bool("FAVG",false);
    /*- Do canonicalize the active orbitals such that the average Fock matrix is diagonal? -*/
    options.add_bool("CANONICALIZE_ACTIVE_FAVG",false);
    /*- Do canonicalize the inactive (DOCC and Virtual) orbitals such that the average Fock matrix is diagonal? -*/
    options.add_bool("CANONICALIZE_INACTIVE_FAVG",false);
    /*- Do consider internal rotations? -*/
    options.add_bool("INTERNAL_ROTATIONS",true);
    /*- Do attempt to force a two configruation solution by starting with CI coefficents of $\pm \sqrt{\frac{1}{2}}$ ? -*/
    options.add_bool("FORCE_TWOCON",false);
    /*- The number of singly occupied orbitals, per irrep -*/
    options.add("SOCC", new ArrayType());
    /*- The number of doubly occupied orbitals, per irrep -*/
    options.add("DOCC", new ArrayType());
    /*- The symmetry of the SCF wavefunction.-*/
    options.add_str("WFN_SYM","1","A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");
  }
  if(name == "CCENERGY"|| options.read_globals()) {
    /*- MODULEDESCRIPTION Computes coupled cluster energies. Called as part of any coupled cluster computation. -*/

    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "NONE", "CCSD CCSD_T EOM_CCSD LEOM_CCSD BCCD BCCD_T CC2 CC3 EOM_CC2 EOM_CC3 CCSD_MVD");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF", "RHF ROHF UHF");
    /*- Do use new triples? -*/
    options.add_bool("NEW_TRIPLES", 1);
    /*- Do analyze T2 amplitudes -*/
    options.add_bool("ANALYZE", 0);
    /*- Maximum number of iterations to solve the CC equations -*/
    options.add_int("MAXITER", 50);
    /*- Convergence criterion for wavefunction (change) in CC amplitude equations. -*/
    options.add_double("R_CONVERGENCE", 1e-7);
    /*- Do restart the coupled-cluster iterations from old $t@@1$ and $t@@2$
    amplitudes?  For geometry optimizations, Brueckner
    calculations, etc. the iterative solution of the CC amplitude
    equations may benefit considerably by reusing old vectors as initial
    guesses.  Assuming that the MO phases remain the same between
    updates, the CC codes will, by default, re-use old vectors, unless
    the user sets RESTART = false. -*/
    options.add_bool("RESTART",1);
    /*- Do restart the coupled-cluster iterations even if MO phases are screwed up? !expert -*/
    options.add_bool("FORCE_RESTART", 0);
//#warning CCEnergy ao_basis keyword type was changed.
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms
    If AO_BASIS is ``NONE``, the MO-basis integrals will be used;
    if AO_BASIS is ``DISK``, the AO-basis integrals stored on disk will
    be used; if AO_BASIS is ``DIRECT``, the AO-basis integrals will be computed
    on the fly as necessary.  NB: The ``DIRECT`` option is not fully
    implemented and should only be used by experts.  Default is NONE.
    Note: The developers recommend use of this keyword only as a last
    resort because it significantly slows the calculation. The current
    algorithms for handling the MO-basis four-virtual-index integrals have
    been significantly improved and are preferable to the AO-based approach.
    !expert -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- Cacheing level for libdpd governing the storage of amplitudes,
    integrals, and intermediates in the CC procedure. A value of 0 retains
    no quantities in cache, while a level of 6 attempts to store all
    quantities in cache.  For particularly large calculations, a value of
    0 may help with certain types of memory problems.  The default is 2,
    which means that all four-index quantites with up to two virtual-orbital
    indices (e.g., $\langle ij | ab \rangle>$ integrals) may be held in the cache. -*/
    options.add_int("CACHELEVEL", 2);
    /*- Selects the priority type for maintaining the automatic memory
    cache used by the libdpd codes. A value of ``LOW`` selects a "low priority"
    scheme in which the deletion of items from the cache is based on
    pre-programmed priorities. A value of LRU selects a "least recently used"
    scheme in which the oldest item in the cache will be the first one deleted. -*/
    options.add_str("CACHETYPE", "LOW", "LOW LRU");
    /*- Number of threads -*/
    options.add_int("CC_NUM_THREADS",1);
    /*- Do use DIIS extrapolation to accelerate convergence? -*/
    options.add_bool("DIIS", true);
    /*- -*/
    options.add_bool("T2_COUPLED", false);
    /*- The response property desired.  Acceptable values are ``POLARIZABILITY``
    (default) for dipole-polarizabilities, ``ROTATION`` for specific rotations,
    ``ROA`` for Raman Optical Activity, and ``ALL`` for all of the above. -*/
    options.add_str("PROPERTY", "POLARIZABILITY", "POLARIZABILITY ROTATION MAGNETIZABILITY ROA ALL");
    /*- Type of ABCD algorithm will be used -*/
    options.add_str("ABCD", "NEW", "NEW OLD");
    /*- Do simulate the effects of local correlation techniques? -*/
    options.add_bool("LOCAL", 0);
    /*- Value (always between one and zero) for the Broughton-Pulay completeness
    check used to contruct orbital domains for local-CC calculations. See
    J. Broughton and P. Pulay, J. Comp. Chem. 14, 736-740 (1993) and C. Hampel
    and H.-J. Werner, J. Chem. Phys. 104, 6286-6297 (1996). -*/
    options.add_double("LOCAL_CUTOFF", 0.02);
    /*- Type of local-CCSD scheme to be simulated. ``WERNER`` selects the method
    developed by H.-J. Werner and co-workers, and ``AOBASIS`` selects the method
    developed by G.E. Scuseria and co-workers (currently inoperative). -*/
    options.add_str("LOCAL_METHOD", "WERNER", "WERNER AOBASIS");
    /*- Desired treatment of "weak pairs" in the local-CCSD method. A value of
    ``NEGLECT`` ignores weak pairs entirely. A value of ``NONE`` treats weak pairs in
    the same manner as strong pairs. A value of MP2 uses second-order perturbation
    theory to correct the local-CCSD energy computed with weak pairs ignored. -*/
    options.add_str("LOCAL_WEAKP", "NONE", "NONE NEGLECT MP2");
    //options.add_int("LOCAL_FILTER_SINGLES", 1);
    /*- Cutoff value for local-coupled-perturbed-Hartree-Fock -*/
    options.add_double("LOCAL_CPHF_CUTOFF", 0.10);
    /*- Definition of local pair domains, default is BP, Boughton-Pulay. -*/
    options.add_str("LOCAL_PAIRDEF", "BP", "BP RESPONSE");
    /*- Number of important $t@@1$ and $t@@2$ amplitudes to print -*/
    options.add_int("NUM_AMPS_PRINT", 10);
    /*- Convergence criterion for Breuckner orbitals. The convergence
       is determined based on the largest $T_1$ amplitude. -*/
    options.add_double("BRUECKNER_ORBS_R_CONVERGENCE", 1e-5);
    /*- Do print the MP2 amplitudes which are the starting guesses for RHF and UHF reference functions? -*/
    options.add_bool("MP2_AMPS_PRINT", 0);
    /*- Do print MP2 and CCSD pair energies for RHF references? -*/
    options.add_bool("PAIR_ENERGIES_PRINT", 0);
    /*- Do print spin-adapted pair energies? -*/
    options.add_bool("SPINADAPT_ENERGIES", false);
    /*- Do build W intermediates required for cc3 in core memory? -*/
    options.add_bool("T3_WS_INCORE", 0);
    /*- Do SCS-MP2 with parameters optimized for nucleic acids? -*/
    options.add_bool("SCSN_MP2", 0);
    /*- Do spin-component-scaled MP2 (SCS-MP2)? -*/
    options.add_bool("SCS_MP2", 0);
    /*- Do spin-component-scaled CCSD -*/
    options.add_bool("SCS_CCSD", 0);
    /*- MP2 opposite-spin scaling value -*/
    options.add_double("MP2_OS_SCALE",1.20);
    /*- MP2 same-spin scaling value -*/
    options.add_double("MP2_SS_SCALE",1.0/3.0);
    /*- Coupled-cluster opposite-spin scaling value -*/
    options.add_double("CC_OS_SCALE", 1.27);
    /*- Coupled-cluster same-spin scaling value -*/
    options.add_double("CC_SS_SCALE",1.13);
    /*- Convert ROHF MOs to semicanonical MOs -*/
    options.add_bool("SEMICANONICAL", true);
  }
  if(name == "CIS"|| options.read_globals()) {
    /*- MODULEDESCRIPTION Performs configuration interaction singles (CIS) computations. Currently unused in
        Psi4. -*/

    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "CIS", "CCSD CCSD_T EOM_CCSD CIS");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF", "RHF ROHF UHF");
    /*- Cutoff value for printing local amplitudes -*/
    options.add_double("LOCAL_AMPS_PRINT_CUTOFF", 0.60);
    /*- Maximum number of iterations -*/
    options.add_int("MAXITER", 500);
    /*- Convergence criterion for CIS wavefunction. -*/
    options.add_double("R_CONVERGENCE", 1e-7);
    /*- The number of electronic states to computed, per irreducible
    representation-*/
    options.add("ROOTS_PER_IRREP", new ArrayType());
    /*- Diagonalization method for the CI matrix -*/
    options.add_str("DIAG_METHOD", "DAVIDSON", "DAVIDSON FULL");
    /*- Do simulate the effects of local correlation techniques? -*/
    options.add_bool("LOCAL", false);
    /*- Value (always between one and zero) for the Broughton-Pulay completeness
    check used to contruct orbital domains for local-CC calculations. See
    J. Broughton and P. Pulay, J. Comp. Chem. 14, 736-740 (1993) and C. Hampel
    and H.-J. Werner, J. Chem. Phys. 104, 6286-6297 (1996). -*/
    options.add_double("LOCAL_CUTOFF", 0.02);
    /*- Type of local-CIS scheme to be simulated. ``WERNER`` selects the method
    developed by H.-J. Werner and co-workers, and ``AOBASIS`` selects the method
    developed by G.E. Scuseria and co-workers. -*/
    options.add_str("LOCAL_METHOD", "WERNER", "AOBASIS WERNER");
    /*- Desired treatment of "weak pairs" in the local-CIS method. A value of
    ``NEGLECT`` ignores weak pairs entirely. A value of ``NONE`` treats weak pairs in
    the same manner as strong pairs. A value of MP2 uses second-order perturbation
    theory to correct the local-CIS energy computed with weak pairs ignored. -*/
    options.add_str("LOCAL_WEAKP", "MP2", "MP2 NEGLECT NONE");
    /*- -*/
    options.add_int("LOCAL_GHOST", -1);
    /*- -*/
    options.add("DOMAINS", new ArrayType());
    /*- Do print the domains? -*/
    options.add_bool("DOMAIN_PRINT", 0);
  }
  if(name == "LMP2"|| options.read_globals()) {
    /*- MODULEDESCRIPTION Performs local MP2 computations for RHF reference functions. -*/

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
    /*- Convergence criterion for energy (change). -*/
    options.add_double("E_CONVERGENCE", 1e-7);
    /*- Convergence criterion for T2 amplitudes (RMS change). -*/
    options.add_double("R_CONVERGENCE", 1e-5);
    /*- Minimum absolute value below which parts of the Fock matrix are skipped. -*/
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
    /*- The amount of memory available (in Mb) -*/
    options.add_int("MEMORY", 2000);
    /*- Do spin-component-scaled MP2 (SCS-MP2)? -*/
    options.add_bool("SCS", false);
    /*- Do SCS-MP2 with parameters optimized for nucleic acids? -*/
    options.add_bool("SCS_N", false);
    /*- The scale factor used for opposite-spin pairs in SCS computations -*/
    options.add_double("MP2_OS_SCALE", 6.0/5.0);
    /*- The scale factor used for same-spin pairs in SCS computations-*/
    options.add_double("MP2_SS_SCALE", 1.0/3.0);
    /*- Do screen integrals? -*/
    options.add_bool("SCREEN_INTS", false);
    /*- Minimum absolute value below which integrals are neglected. -*/
    options.add_double("INTS_TOLERANCE", 1e-7);
    /*- Do exit after printing the domains? -*/
    options.add_bool("DOMAIN_PRINT_EXIT", 0);
   }
  if(name == "DFMP2"|| options.read_globals()) {
    /*- MODULEDESCRIPTION Performs density-fitted MP2 computations for RHF/UHF/ROHF reference wavefunctions. -*/

    /*- A helpful option, used only in debugging the MADNESS version !expert-*/
    options.add_int("MADMP2_SLEEP", 0);
    /*- Primary basis set -*/
    options.add_str("BASIS","NONE");
    /*- Auxiliary basis set for MP2 density fitting computations.
    :ref:`Defaults <apdx:basisFamily>` to a RI basis. -*/
    options.add_str("DF_BASIS_MP2","");
    /*- OS Scale -*/
    options.add_double("MP2_OS_SCALE", 6.0/5.0);
    /*- SS Scale  -*/
    options.add_double("MP2_SS_SCALE", 1.0/3.0);
    /*- \% of memory for DF-MP2 three-index buffers -*/
    options.add_double("DFMP2_MEM_FACTOR", 0.9);
    /*- Minimum absolute value below which integrals are neglected. -*/
    options.add_double("INTS_TOLERANCE", 0.0);
    /*- Minimum error in the 2-norm of the P(2) matrix for corrections to Lia and P. -*/
    options.add_double("DFMP2_P2_TOLERANCE", 0.0);
    /*- Minimum error in the 2-norm of the P matrix for skeleton-core Fock matrix derivatives. -*/
    options.add_double("DFMP2_P_TOLERANCE", 0.0);
    /*- Number of threads to compute integrals with. 0 is wild card -*/
    options.add_int("DF_INTS_NUM_THREADS", 0);
    /*- IO caching for CP corrections, etc !expert -*/
    options.add_str("DF_INTS_IO", "NONE", "NONE SAVE LOAD");
  }
  if(name == "PSIMRCC"|| options.read_globals()) {
    /*- MODULEDESCRIPTION Performs multireference coupled cluster computations.  This theory should be used only by
        advanced users with a good working knowledge of multireference techniques. -*/

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
    options.add_int("CC_NUM_THREADS",1);
    /*- Which root of the effective hamiltonian is the target state? -*/
    options.add_int("FOLLOW_ROOT",1);
    /*- Convergence criterion for energy. -*/
    options.add_double("E_CONVERGENCE",1e-9);
    /*- Convergence criterion for amplitudes (residuals). -*/
    options.add_double("R_CONVERGENCE",1e-9);
    /*- Maximum number of iterations to determine the amplitudes -*/
    options.add_int("MAXITER",100);
    /*- The number of DIIS vectors needed before extrapolation is performed -*/
    options.add_int("DIIS_START",2);
    /*- The shift to apply to the denominators, {\it c.f.} Taube and Bartlett, JCP, 130, 144112 (2009) -*/
    options.add_double("TIKHONOW_OMEGA",0.0);  // Omega = TIKHONOW_OMEGA / 1000
    /*- The cycle after which Tikhonow regularization is stopped.
    Set to zero to allow regularization in all iterations -*/
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
    /*- The ansatz to use for MRCC computations -*/
    options.add_str("CORR_ANSATZ","MK","SR MK BW APBW");
    /*- The order of coupling terms to include in MRCCSDT computations -*/
    options.add_str("COUPLING","CUBIC","NONE LINEAR QUADRATIC CUBIC");
    /*- The symmetry of the target wavefunction, specified either by Sch\ |o_dots|\ nflies symbol,
        or irrep number (in Cotton ordering) -*/
    options.add_str("WFN_SYM","1","A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");
    /*- The type of algorithm to use for (T) computations -*/
    options.add_str("TRIPLES_ALGORITHM","RESTRICTED","SPIN_ADAPTED RESTRICTED UNRESTRICTED");
    /*- How to perform MP2_CCSD computations -*/
    options.add_str("MP2_CCSD_METHOD","II","I IA II");
    /*- Whether to use spin symmetry to map equivalent configurations onto each other, for efficiency !expert -*/
    options.add_bool("USE_SPIN_SYMMETRY", true);
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
    /*- Do disregard updating single excitation amplitudes? -*/
    options.add_bool("NO_SINGLES", false);
  }
  if(name == "OPTKING"|| options.read_globals()) {
    /*- MODULEDESCRIPTION Performs geometry optimizations and vibrational frequency analyses. -*/

      /*- SUBSECTION Optimization Algorithm -*/

      /*- Maximum number of geometry optimization steps -*/
      options.add_int("GEOM_MAXITER", 20);
      /*- Specifies minimum search, transition-state search, or IRC following -*/
      options.add_str("OPT_TYPE", "MIN", "MIN TS IRC");
      /*- Geometry optimization step type, either Newton-Raphson or Rational Function Optimization -*/
      options.add_str("STEP_TYPE", "RFO", "RFO NR SD");
      /*- Do follow the initial RFO vector after the first step? -*/
      options.add_bool("RFO_FOLLOW_ROOT", false);
      /*- Root for RFO to follow, 0 being lowest (for a minimum) -*/
      options.add_int("RFO_ROOT", 0);
      /*- IRC step size in bohr(amu)\ $^{1/2}$. -*/
      options.add_double("IRC_STEP_SIZE", 0.2);
      /*- IRC mapping direction -*/
      options.add_str("IRC_DIRECTION", "FORWARD", "FORWARD BACKWARD");
      /*- Decide when to stop IRC calculations -*/
      options.add_str("IRC_STOP", "STOP", "ASK STOP GO");
      /*- Initial maximum step size in bohr or radian along an internal coordinate -*/
      options.add_double("INTRAFRAG_STEP_LIMIT", 0.4);
      /*- Lower bound for dynamic trust radius [au] -*/
      options.add_double("INTRAFRAG_STEP_LIMIT_MIN", 0.001);
      /*- Upper bound for dynamic trust radius [au] -*/
      options.add_double("INTRAFRAG_STEP_LIMIT_MAX", 1.0);
      /*- Maximum step size in bohr or radian along an interfragment coordinate -*/
      options.add_double("INTERFRAG_STEP_LIMIT", 0.4);
      /*- Set number of consecutive backward steps allowed in optimization -*/
      options.add_int("CONSECUTIVE_BACKSTEPS", 0);

      /*- SUBSECTION Convergence Control -*/

      /*- Set of optimization criteria. Specification of any MAX_*_G_CONVERGENCE 
      or RMS_*_G_CONVERGENCE options will append to overwrite the criteria set here
      unless |optking__flexible_g_convergence| is also on. 
      See Table :ref:`Geometry Convergence <table:optkingconv>` for details. -*/
      options.add_str("G_CONVERGENCE", "QCHEM", "QCHEM MOLPRO GAU GAU_LOOSE GAU_TIGHT GAU_VERYTIGHT TURBOMOLE CFOUR NWCHEM_LOOSE");
      /*- Convergence criterion for geometry optmization: maximum force 
      (internal coordinates, atomic units). -*/
      options.add_double("MAX_FORCE_G_CONVERGENCE", 3.0e-4);
      /*- Convergence criterion for geometry optmization: rms force 
      (internal coordinates, atomic units). -*/
      options.add_double("RMS_FORCE_G_CONVERGENCE", 3.0e-4);
      /*- Convergence criterion for geometry optmization: maximum energy change. -*/
      options.add_double("MAX_ENERGY_G_CONVERGENCE", 1.0e-6);
      /*- Convergence criterion for geometry optmization: maximum displacement 
      (internal coordinates, atomic units). -*/
      options.add_double("MAX_DISP_G_CONVERGENCE", 1.2e-3);
      /*- Convergence criterion for geometry optmization: rms displacement 
      (internal coordinates, atomic units). -*/
      options.add_double("RMS_DISP_G_CONVERGENCE", 1.2e-3);
      /*- Even if a user-defined threshold is set, allow for normal, flexible convergence criteria -*/
      options.add_bool("FLEXIBLE_G_CONVERGENCE", false);

      /*- SUBSECTION Hessian Update -*/

      /*- Hessian update scheme -*/
      options.add_str("HESS_UPDATE", "BFGS", "NONE BFGS MS POWELL BOFILL");
      /*- Number of previous steps to use in Hessian update, 0 uses all -*/
      options.add_int("HESS_UPDATE_USE_LAST", 1);
      /*- Do limit the magnitude of changes caused by the Hessian update? -*/
      options.add_bool("HESS_UPDATE_LIMIT", true);
      /*- If |optking__hess_update_limit| is true, changes to the Hessian
      from the update are limited to the larger of
      |optking__hess_update_limit_scale| * (the previous value) and
      HESS_UPDATE_LIMIT_MAX [au]. -*/
      options.add_double("HESS_UPDATE_LIMIT_MAX", 1.00);
      /*- If |optking__hess_update_limit| is true, changes to the Hessian
      from the update are limited to the larger of HESS_UPDATE_LIMIT_SCALE
      * (the previous value) and |optking__hess_update_limit_max| [au]. -*/
      options.add_double("HESS_UPDATE_LIMIT_SCALE", 0.50);
      /*- Do read Cartesian Hessian?  Only for experts - use
      |optking__full_hess_every| instead. -*/
      options.add_bool("CART_HESS_READ", false);
      /*- Frequency with which to compute the full Hessian in the course
      of a geometry optimization. 0 means to compute the initial Hessian only, 1
      means recompute every step, and N means recompute every N steps. The
      default (-1) is to never compute the full Hessian. -*/
      options.add_int("FULL_HESS_EVERY", -1);
      /*- Model Hessian to guess intrafragment force constants -*/
      options.add_str("INTRAFRAG_HESS", "SCHLEGEL", "FISCHER SCHLEGEL SIMPLE LINDH");

      /*- SUBSECTION Fragment/Internal Coordinate Control -*/

      /*- For multi-fragment molecules, treat as single bonded molecule
      or via interfragment coordinates. A primary difference is that in ``MULTI`` mode,
      the interfragment coordinates are not redundant. -*/
      options.add_str("FRAG_MODE", "SINGLE", "SINGLE MULTI");
      /*- Do freeze all fragments rigid? -*/
      options.add_bool("FREEZE_INTRAFRAG", false);
      /*- Do freeze all interfragment modes? -*/
      options.add_bool("FREEZE_INTERFRAG", false);
      /*- When interfragment coordinates are present, use as reference points either
      principal axes or fixed linear combinations of atoms. -*/
      options.add_str("INTERFRAG_MODE", "FIXED", "FIXED INTERFRAGMENT");
      /*- Do add bond coordinates at nearby atoms for non-bonded systems? -*/
      options.add_bool("ADD_AUXILIARY_BONDS", false);
      /*- Do use $\frac{1}{R@@{AB}}$ for the stretching coordinate between fragments? 
      Otherwise, use $R@@{AB}$. -*/
      options.add_bool("INTERFRAG_DIST_INV", false);
      /*- Model Hessian to guess interfragment force constants -*/
      options.add_str("INTERFRAG_HESS", "DEFAULT", "DEFAULT FISCHER_LIKE");
      /*- When determining connectivity, a bond is assigned if interatomic distance
      is less than (this number) * sum of covalent radii. -*/
      options.add_double("COVALENT_CONNECT", 1.3);
      /*- For now, this is a general maximum distance for the definition of H-bonds -*/
      options.add_double("H_BOND_CONNECT", 4.3);
      /*- Do only generate the internal coordinates and then stop? -*/
      options.add_bool("INTCOS_GENERATE_EXIT", false);

      /*- SUBSECTION Misc. -*/

      /*- Do save and print the geometry from the last projected step at the end
      of a geometry optimization? Otherwise (and by default), save and print
      the previous geometry at which was computed the gradient that satisfied
      the convergence criteria. -*/
      options.add_bool("FINAL_GEOM_WRITE", false);
      /*- Do test B matrix? -*/
      options.add_bool("TEST_B", false);
      /*- Do test derivative B matrix? -*/
      options.add_bool("TEST_DERIVATIVE_B", false);
      /*- Keep internal coordinate definition file. -*/
      options.add_bool("KEEP_INTCOS", false);
      /*- In constrained optimizations, for internal coordinates with user-specified
      equilibrium values, this is the force constant (in au) used to apply an additional
      force to each coordinate.  If the user is only concerned to satify the desired constraint,
      then the user need only ensure that this value is sufficiently large.  Alternatively,
      the user may specify this value to apply a force of a particular magnitude, in which case the
      given equilibrium value may or may not be reached by the optimization. -*/
      options.add_double("INTCO_FIXED_EQ_FORCE_CONSTANT", 2.0);
  }
  if(name == "FINDIF"|| options.read_globals()) {
    /*- MODULEDESCRIPTION Performs finite difference computations of energy derivative, with respect to nuclear displacements
        for geometry optimizations and vibrational frequency analyses, where the required analytical derivatives are not
        available. -*/

      /*- Number of points for finite-differences (3 or 5) -*/
      options.add_int("POINTS", 3); // Can we error check integers?
      /*- Displacement size in au for finite-differences. -*/
      options.add_double("DISP_SIZE", 0.005);
  }
  if (name == "OMP2"|| options.read_globals()) {
    /*- MODULEDESCRIPTION Performs orbital-optimized MP2 computations. -*/

    //options.add_int("MEMORY", 256);
    //options.add_str("REFERENCE", "UHF", "UHF");

    /*- Convergence criterion for energy. -*/
    options.add_double("E_CONVERGENCE",1e-8);
    /*- Convergence criterion for amplitudes (residuals). -*/
    options.add_double("R_CONVERGENCE",1e-5);
    /*- Convergence criterion for RMS orbital gradient. -*/
    options.add_double("RMS_MOGRAD_CONVERGENCE",1e-5);
    /*- Convergence criterion for maximum orbital gradient -*/
    options.add_double("MAX_MOGRAD_CONVERGENCE",1e-4);
    /*- Maximum number of iterations to determine the amplitudes -*/
    options.add_int("CC_MAXITER",50);
    /*- Maximum number of iterations to determine the orbitals -*/
    options.add_int("MO_MAXITER",50);
    /*- Cacheing level for libdpd governing the storage of amplitudes,
    integrals, and intermediates in the CC procedure. A value of 0 retains
    no quantities in cache, while a level of 6 attempts to store all
    quantities in cache.  For particularly large calculations, a value of
    0 may help with certain types of memory problems.  The default is 2,
    which means that all four-index quantites with up to two virtual-orbital
    indices (e.g., $\langle ij | ab \rangle>$ integrals) may be held in the cache. -*/
    options.add_int("CACHELEVEL",2);
    /*- Number of vectors used in DIIS -*/
    options.add_int("DIIS_MAX_VECS",4);
    /*- Cutoff value for numerical procedures -*/
    options.add_int("CUTOFF",14);

    /*- Maximum step size in orbital-optimization procedure -*/
    options.add_double("MO_STEP_MAX",0.5);
    /*- Level shift to aid convergence -*/
    options.add_double("LEVEL_SHIFT",0.02);
    /*- MP2 opposite-spin scaling value -*/
    options.add_double("MP2_OS_SCALE",6.0/5.0);
    /*- MP2 same-spin scaling value -*/
    options.add_double("MP2_SS_SCALE",1.0/3.0);
    /*- Spin-opposite scaling (SOS) value for SCF orbitals -*/
    options.add_double("SOS_SCALE",1.3); 
    /*- Spin-opposite scaling (SOS) value for optimized-MP2 orbitals -*/
    options.add_double("SOS_SCALE2",1.2); 
    //options.add_str("LINEQ","CDGESV","CDGESV FLIN POPLE");
    /*- The algorithm for orthogonalization of MOs -*/
    options.add_str("ORTH_TYPE","MGS","GS MGS");
    //options.add_str("STABILITY","FALSE","TRUE FALSE");
    /*- Do compute natural orbitals? -*/
    options.add_bool("NAT_ORBS",false);
    /*- The optimization algorithm -*/
    options.add_str("OPT_METHOD","DIIS","SD DIIS");
    /*- Type Hessian matrix will be used in orbital optimization procedure -*/
    options.add_str("HESS_TYPE","NONE","NONE");
    /*- Do print OMP2 orbital energies? -*/
    options.add_bool("OMP2_ORBS_PRINT",false);
    /*- Do perform spin-component-scaled OMP2 (SCS-OMP2)? In all computation, SCS-OMP2 energy is computed automatically. 
     However, in order to perform geometry optimizations and frequency computations with SCS-OMP2, one needs to set 
     'DO_SCS' to true -*/
    options.add_bool("DO_SCS",false);
    /*- Do perform spin-opposite-scaled OMP2 (SOS-OMP2)? In all computation, SOS-OMP2 energy is computed automatically. 
     However, in order to perform geometry optimizations and frequency computations with SOS-OMP2, one needs to set 
     'DO_SOS' to true -*/
    options.add_bool("DO_SOS",false);
    /*- Do write coefficient matrices to external files for direct reading MOs in a subsequent job? -*/
    options.add_bool("MO_WRITE",false);
    /*- Do read coefficient matrices from external files of a previous OMP2 or OMP3 computation? -*/
    options.add_bool("MO_READ",false);
  }
  if (name == "OMP3"|| options.read_globals()) {
    /*- MODULEDESCRIPTION Performs orbital-optimized MP3 computations. -*/

    /*- Convergence criterion for energy. -*/
    options.add_double("E_CONVERGENCE",1e-8);
    /*- Convergence criterion for amplitudes (residuals). -*/
    options.add_double("R_CONVERGENCE",1e-5);
    /*- Convergence criterion for RMS orbital gradient. -*/
    options.add_double("RMS_MOGRAD_CONVERGENCE",1e-5);
    /*- Convergence criterion for maximum orbital gradient -*/
    options.add_double("MAX_MOGRAD_CONVERGENCE",1e-4);
    /*- Maximum number of iterations to determine the amplitudes -*/
    options.add_int("CC_MAXITER",50);
    /*- Maximum number of iterations to determine the orbitals -*/
    options.add_int("MO_MAXITER",50);
    /*- Cacheing level for libdpd governing the storage of amplitudes,
    integrals, and intermediates in the CC procedure. A value of 0 retains
    no quantities in cache, while a level of 6 attempts to store all
    quantities in cache.  For particularly large calculations, a value of
    0 may help with certain types of memory problems.  The default is 2,
    which means that all four-index quantites with up to two virtual-orbital
    indices (e.g., $\langle ij | ab \rangle>$ integrals) may be held in the cache. -*/
    options.add_int("CACHELEVEL",2);
    /*- Number of vectors used in DIIS -*/
    options.add_int("DIIS_MAX_VECS",4);
    /*- Cutoff value for numerical procedures -*/
    options.add_int("CUTOFF",14);

    /*- Maximum step size in orbital-optimization procedure -*/
    options.add_double("MO_STEP_MAX",0.5);
    /*- Level shift to aid convergence -*/
    options.add_double("LEVEL_SHIFT",0.02);
    /*- MP2 opposite-spin scaling value -*/
    options.add_double("MP2_OS_SCALE",6.0/5.0);
    /*- MP2 same-spin scaling value -*/
    options.add_double("MP2_SS_SCALE",1.0/3.0);
    /*- Spin-opposite scaling (SOS) value for SCF orbitals -*/
    options.add_double("SOS_SCALE",1.3);  
    /*- Spin-opposite scaling (SOS) value for optimized-MP2 orbitals -*/
    options.add_double("SOS_SCALE2",1.2);  
    /*- Scaling value for 3rd order energy correction (S. Grimme, Vol. 24, pp. 1529, J. Comput. Chem.) -*/
    options.add_double("E3_SCALE",0.25); 
    /*- The algorithm for orthogonalization of MOs -*/
    options.add_str("ORTH_TYPE","MGS","GS MGS");

    /*- Do compute natural orbitals? -*/
    options.add_bool("NAT_ORBS",false);
    /*- The optimization algorithm -*/
    options.add_str("OPT_METHOD","DIIS","SD DIIS");
    /*- Type Hessian matrix will be used in orbital optimization procedure -*/
    options.add_str("HESS_TYPE","NONE","NONE");
    /*- Do print OMP3 orbital energies? -*/
    options.add_bool("OMP3_ORBS_PRINT",false);
    /*- Do perform spin-component-scaled OMP3 (SCS-OMP3)? In all computation, SCS-OMP3 energy is computed automatically. 
     However, in order to perform geometry optimizations and frequency computations with SCS-OMP3, one needs to set 
     'DO_SCS' to true -*/
    options.add_bool("DO_SCS",false);
    /*- Do perform spin-opposite-scaled OMP3 (SOS-OMP3)? In all computation, SOS-OMP3 energy is computed automatically. 
     However, in order to perform geometry optimizations and frequency computations with SOS-OMP3, one needs to set 
     'DO_SOS' to true -*/
    options.add_bool("DO_SOS",false);
    /*- Do write coefficient matrices to external files for direct reading MOs in a subsequent job? -*/
    options.add_bool("MO_WRITE",false);
    /*- Do read coefficient matrices from external files of a previous OMP2 or OMP3 computation? -*/
    options.add_bool("MO_READ",false);
  }
  if (name == "MRCC"|| options.read_globals()) {
      /*- MODULEDESCRIPTION Interface to MRCC program written by Mih\ |a_acute|\ ly K\ |a_acute|\ llay. -*/

      /*- Sets the OMP_NUM_THREADS environment variable before calling MRCC.
          If the environment variable :envvar:`OMP_NUM_THREADS` is set prior to calling PSI4 then
          that value is used. When set, this option overrides everything. Be aware
          the ``-n`` command-line option described in section :ref:`sec:threading`
          does not affect MRCC.
          !expert -*/
      options.add_int("MRCC_OMP_NUM_THREADS", 1);

      /*- This becomes ``tol`` (option \#16) in fort.56. -*/
      options.add_double("E_CONVERGENCE",1e-8);

      /*- Minimum absolute value below which integrals are neglected. -*/
      options.add_double("INTS_TOLERANCE",1.0E-12);

      /*- Maximum excitation level. This is used ONLY if it is explicity set by the user.
        Single-reference case: all excitations up to this level are
        included, e.g., 2 for CCSD, 3 for CCSDT, 4 for CCSDTQ, etc.
        This becomes ``ex.lev`` (option \#1) in fort.56. -*/
      options.add_int("MRCC_LEVEL", 2);

      /*- Number of singlet roots. (Strictly speaking number of
          of roots with M_s=0 and S is even.) Use this option only with
          closed shell reference determinant, it must be zero otherwise.
          This becomes ``nsing`` (option \#2) in fort.56. -*/
      options.add_int("MRCC_NUM_SINGLET_ROOTS", 1);

      /*- Number of triplet roots. (Strictly speaking number of of roots
          with $M_s=0$ and S is odd.) See notes at option
          |mrcc__mrcc_num_singlet_roots|. This becomes ``ntrip`` (option \#3)
          in fort.56. -*/
      options.add_int("MRCC_NUM_TRIPLET_ROOTS", 0);

      /*- The program restarts from the previously
          calculated parameters if it is 1. In case it is 2, the program
          executes automatically the lower-level calculations of the same
          type consecutively (e.g., CCSD, CCSDT, and CCSDTQ if CCSDTQ is
          requested) and restarts each calculation from the previous one
          (rest=2 is available only for energy calculations).
          Currently, only a value of 0 and 2 are supported.
          This becomes ``rest`` (option \#4) in fort.56. !expert -*/
      options.add_int("MRCC_RESTART", 0);

      /*- If more than one root is requested and calc=1, LR-CC (EOM-CC)
          calculation is performed automatically for the excited states.
          This overrides all automatic determination of method
          and will only work with :py:func:`~driver.energy`.
          This becomes CC/CI (option \#5) in fort.56
          |  \begin{tabular}{ccc}
          |         Value  &  Method      &  Description  \\
          |   \hline
          |              1 & CC           & \\
          |              2 & CC(n-1)[n]   & \\
          |              3 & CC(n-1)(n)   &  (CC(n-1)[n] energy is also calculated) \\
          |              4 & CC(n-1)(n)_L & (CC(n-1)[n] and CC(n-1)(n) energies are also calculated) \\
          |              5 & CC(n)-1a     & \\
          |              6 & CC(n)-1b     & \\
          |              7 & CCn          & \\
          |              8 & CC(n)-3      & \\
          |  \end{tabular}
            !expert
          -*/
      options.add_int("MRCC_METHOD", 1);
  }
  if (name == "CEPA"|| options.read_globals()) {
      /*- Desired convergence for the t1 and t2 amplitudes, defined as
      the norm of the change in the amplitudes between iterations.-*/
      options.add_double("R_CONVERGENCE", 1.0e-7);
      /*- Maximum number of iterations to converge the t1 and t2
      amplitudes. -*/
      options.add_int("MAXITER", 100);
      /*- Number of vectors to store for DIIS extrapolation. -*/
      options.add_int("DIIS_MAX_VECS", 8);
      /*- Opposite-spin scaling factor for SCS-MP2. -*/
      options.add_double("MP2_SCALE_OS",1.20);
      /*- Same-spin scaling factor for SCS-MP2-*/
      options.add_double("MP2_SCALE_SS",1.0/3.0);
      /*- Perform SCS-CEPA? If true, note that the
      default values for the spin component scaling factors
      are optimized for the CCSD method. -*/
      options.add_bool("SCS_CEPA", false);
      /*- Oppposite-spin scaling factor for SCS-CEPA. -*/
      options.add_double("CEPA_SCALE_OS", 1.27);
      /*- Same-spin scaling factor for SCS-CEPA. -*/
      options.add_double("CEPA_SCALE_SS",1.13);
      /*- Which coupled-pair method is called?  This parameter is
      used internally by the python driver.  Changing its value 
      won't have any effect on the procedure. -*/
      options.add_str("CEPA_LEVEL","CEPA0");
      /*- Compute the dipole moment? Note that quadrupole moments
      will also be computed if PRINT >= 2. -*/
      options.add_bool("DIPMOM",false);
      /*- Flag to exclude singly excited configurations from the
      computation. Note that this algorithm is not optimized for
      doubles-only computations. -*/
      options.add_bool("CEPA_NO_SINGLES",false);
      /*- Use integral-direct implementation of the (ac|bd) t(ij,cd)
      contraction? AO integrals will be generated on the fly. The 
      CEPA iterations will be slower, but the AO->MO integral 
      transform will be faster, and the out-of-core sort of the 
      (AC|BD) integrals will be avoided. -*/
      options.add_bool("CEPA_VABCD_DIRECT",false);
  }
  if (name == "THERMO"|| options.read_globals()) {
      /*- Temperature in Kelvin for thermodynamic analysis. -*/
      options.add_double("T", 298.15);
      /*- Pressure in Pascal for thermodynamic analysis. -*/
      options.add_double("P", 101325);
  }
  return true;
}

} //end ::psi

