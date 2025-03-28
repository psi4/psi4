/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file read_calculation_options
    \defgroup PSI4
*/

#include "psi4/physconst.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"

#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"

// clang-format off

namespace psi {

/**
 * This is called immediately before a module is run.  Any options
 * expected by that module must be added here
 *
 * @param name    - the name of the module.
 * @param options - the liboptions module used in the computations.
 * @param suppress_printing - boolean to specify whether to print to output file [false]
 */
int read_options(const std::string &name, Options &options, bool suppress_printing) {
    //  options.clear();

    // dodoc == "GLOBALS" fake line to make document_options_and_tests.pl generate a GLOBALS doc section

    /*- An array containing the number of doubly-occupied orbitals per irrep
    (in Cotton order) -*/
    options.add("DOCC", new ArrayType());
    /*- An array containing the number of singly-occupied orbitals per irrep
    (in Cotton order).  The value of |globals__docc| should also be set. -*/
    options.add("SOCC", new ArrayType());
    /*- An array containing the number of frozen doubly-occupied orbitals per
    irrep (these are not excited in a correlated wavefunction, nor can they be
    optimized in MCSCF. This trumps |globals__num_frozen_docc| and
    |globals__freeze_core|. -*/
    options.add("FROZEN_DOCC", new ArrayType());
    /*- An array containing the number of frozen unoccupied orbitals per
    irrep (these are not populated in a correlated wavefunction, nor can they be
    optimized in MCSCF.  This trumps |globals__num_frozen_uocc|. -*/
    options.add("FROZEN_UOCC", new ArrayType());
    /*- The number of core orbitals to freeze in later correlated computations.
    This trumps |globals__freeze_core|.  -*/
    options.add_int("NUM_FROZEN_DOCC", 0);
    /*- The number of virtual orbitals to freeze in later correlated computations. -*/
    options.add_int("NUM_FROZEN_UOCC", 0);

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

    /*- Specifies how many core orbitals to freeze in correlated computations.
    ``TRUE`` or ``1`` will default to freezing the previous noble gas shell
    on each atom. In case of positive charges on fragments, an additional
    shell may be unfrozen, to ensure there are valence electrons in each
    fragment. With ``FALSE`` or ``0``, no electrons are frozen (with the
    exception of electrons treated by an ECP). With ``-1``, ``-2``, and ``-3``,
    the user might request strict freezing of the previous first/second/third
    noble gas shell on every atom. In this case, when there are no valence
    electrons, the code raises an exception. More precise control over the
    number of frozen orbitals can be attained by using the keywords
    |globals__num_frozen_docc| (gives the total number of orbitals to freeze,
    program picks the lowest-energy orbitals) or |globals__frozen_docc| (gives
    the number of orbitals to freeze per irreducible representation) or by
    the option ``POLICY`` in combination with appropriate inputs to
    |globals__freeze_core_policy|. At present, ``POLICY`` is an experimental
    option and is subject to change.-*/
    options.add_str("FREEZE_CORE", "FALSE", "FALSE TRUE 1 0 -1 -2 -3 POLICY");

    /*- NOTE: This is an experimental feature and subject to change! Specifies
    a custom frozen-core policy on a per-element basis. Input should be a list
    of integers representing the number of orbitals to freeze for each atomic
    number MINUS one (so H is 0, He is 1, etc). For example, to specify that
    elements H-Be should have 0 frozen orbitals, B-Mg should have 1, and Al
    should have 2, you would provide the input ``[0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
    1, 1, 2]``. Please make sure to fill in the list up to the highest atomic
    number included in any calculations. This option is only used if
    |globals__freeze_core| is set to ``POLICY``. -*/
    options.add("FREEZE_CORE_POLICY", new ArrayType());

    options.add_int("NUM_GPUS", 1);
    /*- Do use pure angular momentum basis functions?
    If not explicitly set, the default comes from the basis set.
    **Cfour Interface:** Keyword translates into |cfour__cfour_spherical|. -*/
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
    /*- For displacements, symmetry (Schoenflies symbol) of 'parent' (undisplaced)
    reference molecule. Internal use only for finite difference. !expert -*/
    options.add_str("PARENT_SYMMETRY", "");
    /*- Number of columns to print in calls to ``Matrix::print_mat``. !expert -*/
    options.add_int("MAT_NUM_COLUMN_PRINT", 5);
    /*- List of properties to compute -*/
    options.add("PROPERTIES", new ArrayType());
    /*- Either :ref:`a set of 3 coordinates or a string <table:oe_origin>`
    describing the origin about which one-electron properties are computed. -*/
    options.add("PROPERTIES_ORIGIN", new ArrayType());

    /*- Psi4 dies if energy does not converge. !expert -*/
    options.add_bool("DIE_IF_NOT_CONVERGED", true);
    /*- Integral package to use. If compiled with Simint support, change this option to use them; LibInt2 is used
       otherwise. -*/
    options.add_str("INTEGRAL_PACKAGE", "LIBINT2", "LIBINT2 SIMINT");
#ifdef USING_BrianQC
    /*- Whether to enable using the BrianQC GPU module -*/
    options.add_bool("BRIANQC_ENABLE", false);
#endif

    // Note that case-insensitive options are only functional as
    //   globals, not as module-level, and should be defined sparingly

    /*- Base filename for text files written by PSI, such as the
    MOLDEN output file, the Hessian file, the internal coordinate file,
    etc. Use the add_str_i function to make this string case sensitive. -*/
    options.add_str_i("WRITER_FILE_LABEL", "");
    /*- The density fitting basis to use in coupled cluster computations. -*/
    options.add_str("DF_BASIS_CC", "");
    /*- Assume external fields are arranged so that they have symmetry. It is up to the user to know what to do here.
       The code does NOT help you out in any way! !expert -*/
    options.add_bool("EXTERNAL_POTENTIAL_SYMMETRY", false);
    /*- Text to be passed directly into CFOUR input files. May contain
    molecule, options, percent blocks, etc. Access through ``cfour {...}``
    block. -*/
    options.add_str_i("LITERAL_CFOUR", "");
    /*- When several modules can compute the same methods and the default
    routing is not suitable, this targets a module. ``CCENERGY`` covers
    CCHBAR, etc. ``OCC`` covers OCC and DFOCC. -*/
    options.add_str("QC_MODULE", "", "CCENERGY DETCI DFMP2 FNOCC OCC CCT3 BUILTIN MRCC");
    /*- What algorithm to use for the SCF computation. See Table :ref:`SCF
    Convergence & Algorithm <table:conv_scf>` for default algorithm for
    different calculation types. -*/
    options.add_str("SCF_TYPE", "PK", "DIRECT DF MEM_DF DISK_DF PK OUT_OF_CORE CD GTFOCK DFDIRJ DFDIRJ+COSX DFDIRJ+LINK DFDIRJ+SNLINK");
    /*- Algorithm to use for MP2 computation.
    See :ref:`Cross-module Redundancies <table:managedmethods>` for details. -*/
    options.add_str("MP2_TYPE", "DF", "DF CONV CD");
    /*- Algorithm to use for MPn ( $n>2$ ) computation (e.g., MP3 or MP2.5 or MP4(SDQ)).
    See :ref:`Cross-module Redundancies <table:managedmethods>` for details.
    Since v1.4, default for non-orbital-optimized MP2.5 and MP3 is DF. -*/
    options.add_str("MP_TYPE", "CONV", "DF CONV CD");
    // The type of integrals to use in coupled cluster computations. DF activates density fitting for the largest
    // integral files, while CONV results in no approximations being made.
    /*- Algorithm to use for CC or CEPA computation (e.g., CCD, CCSD(T), CEPA(3), ACPF, REMP).
    See :ref:`Cross-module Redundancies <table:managedmethods>` for details. -*/
    options.add_str("CC_TYPE", "CONV", "DF CONV CD");
    /*- Algorithm to use for CI computation (e.g., CID or CISD).
    See :ref:`Cross-module Redundancies <table:managedmethods>` for details. -*/
    options.add_str("CI_TYPE", "CONV", "CONV");
    /*- Write all the MOs to the MOLDEN file (true) or discard the unoccupied MOs (false). -*/
    options.add_bool("MOLDEN_WITH_VIRTUAL", true);

    /*- The type of screening used when computing two-electron integrals. -*/
    options.add_str("SCREENING", "CSAM", "SCHWARZ CSAM DENSITY NONE");

    // CDS-TODO: We should go through and check that the user hasn't done
    // something silly like specify frozen_docc in DETCI but not in TRANSQT.
    // That would create problems.  (This was formerly checked in DETCI
    // itself, but I don't think DETCI will have the info available to check
    // this anymore).  This problem has affected users in the past.
    // Same goes for restricted_docc, restricted_uocc, ras1, ras2, ras3,
    // frozen_uocc.

#ifdef USING_dkh
    /*- Relativistic Hamiltonian type !expert -*/
    options.add_str("RELATIVISTIC", "NO", "NO X2C DKH");
#else
    /*- Relativistic Hamiltonian type !expert -*/
    options.add_str("RELATIVISTIC", "NO", "NO X2C");
#endif
    /*- Auxiliary basis set for solving Dirac equation in X2C and DKH
        calculations. Defaults to decontracted orbital basis. -*/
    options.add_str("BASIS_RELATIVISTIC", "");
    /*- Order of Douglas-Kroll-Hess !expert -*/
    options.add_int("DKH_ORDER", 2);

    /*- Directory to which to write cube files. Default is the input file
    directory. -*/
    options.add_str_i("CUBEPROP_FILEPATH", ".");

    /*- Properties to compute. Valid tasks include:
        ``DENSITY`` - Da, Db, Dt, Ds;
        ``ESP`` - Dt, ESP;
        ``ORBITALS`` - Psi_a_N, Psi_b_N;
        ``BASIS_FUNCTIONS`` - Phi_N;
        ``LOL`` - LOLa, LOLb;
        ``ELF`` - ELFa, ELFb;
        ``FRONTIER_ORBITALS`` - Psi_a_N_HOMO + Psi_a_N_LUMO;
        ``DUAL_DESCRIPTOR`` - DUAL_N_HOMO-M_LUMO.
    -*/
    options.add("CUBEPROP_TASKS", new ArrayType());
    /*- List of orbital indices for which cube files are generated (1-based,
    $+$ for alpha, $-$ for beta). All orbitals computed if empty. -*/
    options.add("CUBEPROP_ORBITALS", new ArrayType());
    /*- List of basis function indices for which cube files are generated
    (1-based). All basis functions computed if empty.-*/
    options.add("CUBEPROP_BASIS_FUNCTIONS", new ArrayType());
    /*- Fraction of density captured by adaptive isocontour values -*/
    options.add_double("CUBEPROP_ISOCONTOUR_THRESHOLD", 0.85);
    /*- CubicScalarGrid basis cutoff. !expert -*/
    options.add_double("CUBIC_BASIS_TOLERANCE", 1.0E-12);
    /*- CubicScalarGrid maximum number of grid points per evaluation block. !expert -*/
    options.add_int("CUBIC_BLOCK_MAX_POINTS", 1000);
    /*- CubicScalarGrid spatial extent in bohr [O_X, O_Y, O_Z]. Defaults to 4.0 bohr each. -*/
    options.add("CUBIC_GRID_OVERAGE", new ArrayType());
    /*- CubicScalarGrid grid spacing in bohr [D_X, D_Y, D_Z]. Defaults to 0.2 bohr each. -*/
    options.add("CUBIC_GRID_SPACING", new ArrayType());
    /*- How many NOONS to print -- used in libscf_solver/uhf.cc and libmints/oeprop.cc -*/
    options.add_str("PRINT_NOONS", "3");

    /// Tensor Hypercontration (THC) Options (libmints/thc_eri.cc)

    /*- Use DF approximation when computing LS-THC factorization? -*/
    options.add_bool("LS_THC_DF", true);
    /*- Number of spherical points in LS-THC grid -*/
    options.add_int("LS_THC_SPHERICAL_POINTS", 50);
    /*- Number of radial points in LS-THC grid -*/
    options.add_int("LS_THC_RADIAL_POINTS", 10);
    /*- Screening criteria for basis function values on LS-THC grids !expert -*/
    options.add_double("LS_THC_BASIS_TOLERANCE", 1.0E-10);
    /*- Grid weights cutoff for LS-THC grids !expert -*/
    options.add_double("LS_THC_WEIGHTS_TOLERANCE", 1.0E-12);
    /*- Pruning scheme for LS-THC grids !expert -*/
    options.add_str("LS_THC_PRUNING_SCHEME", "ROBUST", 
                        "ROBUST TREUTLER NONE FLAT P_GAUSSIAN D_GAUSSIAN P_SLATER D_SLATER LOG_GAUSSIAN LOG_SLATER NONE");
    /*- Tolerance for pseudoinversion of grid point overlap matrix (Parrish 2012 eq. 30) !expert -*/
    options.add_double("LS_THC_S_EPSILON", 1.0E-10);

    /// MBIS Options (libmints/oeprop.cc)

    /*- Maximum Number of MBIS Iterations -*/
    options.add_int("MBIS_MAXITER", 500);
    /*- MBIS Convergence Criteria -*/
    options.add_double("MBIS_D_CONVERGENCE", 1.0e-8);
    /*- MBIS Number of Radial Points -*/
    /*- Additional Radial and/or Spherical Points may be needed for Heavier Atoms (200-300) like Zinc -*/
    options.add_int("MBIS_RADIAL_POINTS", 75);
    /*- MBIS Number of Spherical Points -*/
    options.add_int("MBIS_SPHERICAL_POINTS", 302);
    /*- Pruning scheme for MBIS Grid -*/
    options.add_str("MBIS_PRUNING_SCHEME", "ROBUST",
                    "ROBUST TREUTLER NONE FLAT P_GAUSSIAN D_GAUSSIAN P_SLATER D_SLATER LOG_GAUSSIAN LOG_SLATER NONE");
    /*- Maximum Radial Moment to Calculate -*/
    options.add_int("MAX_RADIAL_MOMENT", 4);

    /*- PCM boolean for pcmsolver module -*/
    options.add_bool("PCM", false);
    /*- PE boolean for polarizable embedding module -*/
    options.add_bool("PE", false);
    /*- DDX boolean for ddx module -*/
    options.add_bool("DDX", false);

    if (name == "PCM" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs polarizable continuum model (PCM) computations. -*/

        /*- Use total or separate potentials and charges in the PCM-SCF step. !expert -*/
        options.add_str("PCM_SCF_TYPE", "TOTAL", "TOTAL SEPARATE");
        /*- Name of the PCMSolver input file as parsed by pcmsolver.py !expert -*/
        options.add_str_i("PCMSOLVER_PARSED_FNAME", "");
        /*- PCM-CCSD algorithm type. -*/
        options.add_str("PCM_CC_TYPE", "PTE", "PTE");
    }

    if (name == "DDX" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs continuum solvation model computations using
            the domain-decomposition paradigm. -*/

        /*- Switch available solvation models -*/
        options.add_str("DDX_MODEL", "PCM", "PCM COSMO LPB");

        /*- Radius set for cavity spheres. Ignored if RADII is set. -*/
        options.add_str("DDX_RADII_SET", "UFF", "UFF BONDI");

        /*- Scaling factor for cavity spheres. Ignored if RADII is set. The default depends on
            the radii set chosen. -*/
        options.add_double("DDX_RADII_SCALING", 1.1);

        /*- Custom cavity radii. One per atom, uses the unit of the molecule. -*/
        options.add("DDX_RADII", new ArrayType());

        /*- Solvent to use. Not case sensitive. Ignored if SOLVENT_EPSILON is set. -*/
        options.add_str("DDX_SOLVENT", "");

        /*- Dielectric constant of the solvent to use -*/
        options.add_double("DDX_SOLVENT_EPSILON", 0);

        /*- Optical dielectric constant of the solvent to use for non-equilibrium solvation -*/
        options.add_double("DDX_SOLVENT_EPSILON_OPTICAL", 0);

        /*- Debye-Hückel parameter of the solvent to use. Ignored if DDX_MODEL is not LPB;
            mandatory for LPB. Uses the unit of the molecule (i.e. either ang^{-1} or bohr^{-1}). -*/
        options.add_double("DDX_SOLVENT_KAPPA", 0);

        /*- Maximal degree of modelling spherical harmonics -*/
        options.add_int("DDX_LMAX", 9);

        /*- Number of Lebedev grid points to use.
            (A :ref:`Lebedev Points <table:lebedevorder>` number) -*/
        options.add_int("DDX_N_LEBEDEV", 302);

        /*- Maximal number of iterations used inside DDX -*/
        options.add_int("DDX_MAXITER", 100);

        /*- Number of previous iterates to use in DIIS acceleration inside DDX -*/
        options.add_int("DDX_DIIS_MAX_VECS", 20);

        /*- Tolerance to which DDX linear systems are solved -*/
        options.add_double("DDX_SOLVATION_CONVERGENCE", 1e-8);

        /*- Number of spherical points used to compute the solute electric potential/field
            integrals for DDX calculations (A :ref:`Lebedev Points <table:lebedevorder>` number) -*/
        options.add_int("DDX_SOLUTE_SPHERICAL_POINTS", 110);

        /*- Number of radial points used to compute the integrals for DDX calculations -*/
        options.add_int("DDX_SOLUTE_RADIAL_POINTS", 35);

        /*- Use an in-core version, which uses more memory, but is generally faster -*/
        options.add_bool("DDX_INCORE", false);

        /*- Use the fast multipole method to accelerate the solver -*/
        options.add_bool("DDX_FMM", true);

        /*- Maximal degree of multipole spherical harmonics (far-field FMM interactions).
            Using the same value as |ddx__ddx_lmax| is recommended and done by default. -*/
        options.add_int("DDX_FMM_MULTIPOLE_LMAX", 9);

        /*- Maximal degree of local spherical harmonics (near-field FMM interations). -*/
        options.add_int("DDX_FMM_LOCAL_LMAX", 6);

        /*- Logfile to dump a full trace of the DDX solver history for debugging. !expert -*/
        options.add_str("DDX_LOGFILE", "");

        /*- Regularization parameter for characteristic function of sphere overlap.
            Advanced parameter, which usually does not need to be modified. Valid
            values are within the range [0, 1]. !expert -*/
        options.add_double("DDX_ETA", 0.1);

        /*- Shift for characteristic function of sphere overlap.
            Advanced parameter, which usually does not need to be modified. Valid values
            are within the range [-1, 1] with -100 denoting an automatic selection of the
            best shift. !expert -*/
        options.add_double("DDX_SHIFT", -100.0);
    }

    if (name == "PE" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs polarizable embedding model (PE) computations. -*/

        /*- Name of the potential file OR contents of potential file to be written anonymously on-the-fly. -*/
        options.add_str_i("POTFILE", "potfile.pot");
        /*- Threshold for induced moments convergence -*/
        options.add_double("INDUCED_CONVERGENCE", 1e-8);
        /*- Maximum number of iterations for induced moments -*/
        options.add_int("MAXITER", 50);
        /*- Make polarizabilities isotropic -*/
        options.add_bool("ISOTROPIC_POL", false);
        /*- Enable Thole damping for induced moments -*/
        options.add_bool("DAMP_INDUCED", false);
        /*- Enable Thole damping for multipole fields -*/
        options.add_bool("DAMP_MULTIPOLE", false);
        /*- Thole damping factor for induced moments -*/
        options.add_double("DAMPING_FACTOR_INDUCED", 2.1304);
        /*- Thole damping factor for multipole fields -*/
        options.add_double("DAMPING_FACTOR_MULTIPOLE", 2.1304);

        /*- Summation scheme for field computations, can be direct or fmm -*/
        options.add_str("SUMMATION_FIELDS", "DIRECT", "DIRECT FMM");
        /*- Expansion order of the multipoles for FMM -*/
        options.add_int("TREE_EXPANSION_ORDER", 5);
        /*- Opening angle theta -*/
        options.add_double("TREE_THETA", 0.5);

        /*- Activate border options for sites in proximity to the QM/MM border -*/
        options.add_bool("BORDER", false);
        /*- border type, either remove or redistribute moments/polarizabilities -*/
        options.add_str("BORDER_TYPE", "REMOVE", "REMOVE REDIST");
        /*- minimum radius from QM atoms to MM sites to be taken into account
        for removal/redistribution -*/
        options.add_double("BORDER_RMIN", 2.2);
        /*- unit of BORDER_RMIN, default is atomic units (AU) -*/
        options.add_str("BORDER_RMIN_UNIT", "AU", "AU AA");
        /*- order from which moments are removed, e.g.,
        if set to 1 (default), only charges are redistributed and
        all higher order moments are removed -*/
        options.add_int("BORDER_REDIST_ORDER", 1);
        /*- number of neighbor sites to redistribute to.
        The default (-1) redistributes to all sites which are not in the border region -*/
        options.add_int("BORDER_N_REDIST", -1);
        /*- redistribute polarizabilities? If false, polarizabilities are removed (default) -*/
        options.add_bool("BORDER_REDIST_POL", false);

        /*- use PE(ECP) repulsive potentials -*/
        options.add_bool("PE_ECP", false);
    }

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
        options.add_str("REFERENCE", "RHF", "RHF ROHF");

        /*- Convergence criterion for CI residual vector in the Davidson
        algorithm (RMS error).
        The default is 1e-4 for energies and 1e-7 for gradients. -*/
        options.add_double("R_CONVERGENCE", 1e-4);

        /*- Convergence criterion for energy. See Table :ref:`Post-SCF
        Convergence <table:conv_corl>` for default convergence criteria for
        different calculation types. -*/
        options.add_double("E_CONVERGENCE", 1e-6);

        /*- Maximum number of iterations to diagonalize the Hamiltonian -*/
        options.add_int("CI_MAXITER", 24);

        /*- Do a full CI (FCI)? If TRUE, overrides the value of |detci__ex_level|. -*/
        options.add_bool("FCI", false);

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
        options.add_bool("ISTOP", false);

        /*- Do print a summary of the CI blocks? -*/
        options.add_bool("CIBLKS_PRINT", false);

        /*- Number of important determinants to print -*/
        options.add_int("NUM_DETS_PRINT", 20);

        /*- Do freeze core orbitals? -*/
        // CDS-TODO: Need to make DETCI compatible with normal FREEZE_CORE
        options.add_bool("DETCI_FREEZE_CORE", true);

        /*- Do calculate the value of $\langle S^2\rangle$ for each root?
        Only supported for |detci__icore| = 1. -*/
        options.add_bool("CALC_S_SQUARED", false);

        /*- Specifies how to handle buffering of CI vectors.  A value of 0
        makes the program perform I/O one RAS subblock at a time; 1
        uses entire CI vectors at a time; and 2 uses one irrep block
        at a time.  Values of 0 or 2 cause some inefficiency in the I/O
        (requiring multiple reads of the C vector when constructing
        H in the iterative subspace if |detci__diag_method| = SEM), but require
        less core memory. -*/
        options.add_int("ICORE", 1);

        /*- Number of threads for DETCI. !expert -*/
        options.add_int("CI_NUM_THREADS", 1);

        /*- Do print the sigma overlap matrix?  Not generally useful.  !expert -*/
        options.add_bool("SIGMA_OVERLAP", false);

        /*- Array giving the root numbers of the states to average in a
        state-averaged procedure such as SA-CASSCF. Root numbering starts
        from 0. -*/
        options.add("AVG_STATES", new ArrayType());

        /*- Array giving the weights for each state in a state-averaged
        procedure -*/
        // CDS:TODO - Does this work for doubles??
        options.add("AVG_WEIGHTS", new ArrayType());

        /*- The value of the spin quantum number $S$ is given by this option.
        The default is determined by the value of the multiplicity.  This is used
        for two things: (1) determining the phase of the redundant half of the CI
        vector when the $M@@s = 0$ component is used (i.e., |detci__ms0| = ``TRUE``), and (2) making
        sure the guess vector has the desired value of $\langle S^2\rangle$
        (if |detci__calc_s_squared| is ``TRUE`` and |detci__icore| = ``1``). -*/
        options.add_double("S", 0.0);

        /*- Do use the $M@@s = 0$ component of the state? Defaults to TRUE
        if closed-shell and FALSE otherwise. Related to the |detci__s| option. -*/
        options.add_bool("MS0", false);

        /*- An array of length |detci__ex_level| specifying whether each excitation type
        (S,D,T, etc.) is allowed (1 is allowed, 0 is disallowed).  Used to
        specify non-standard CI spaces such as CIST.  !expert -*/
        options.add("EX_ALLOW", new ArrayType());

        /*- Do eliminate determinants not valid for spin-complete spin-flip CI's?
        [see J. S. Sears et al, J. Chem. Phys. 118, 9084-9094 (2003)] !expert -*/
        options.add_bool("SF_RESTRICT", false);

        /*- maximum number of alpha electrons in RAS III -*/
        options.add_int("A_RAS3_MAX", -1);

        /*- maximum number of beta electrons in RAS III -*/
        options.add_int("B_RAS3_MAX", -1);

        /*- maximum number of electrons in RAS III -*/
        options.add_int("RAS3_MAX", -1);

        /*- maximum number of electrons in RAS IV -*/
        options.add_int("RAS4_MAX", -1);

        /*- maximum number of electrons in RAS III + IV -*/
        options.add_int("RAS34_MAX", -1);

        /*- Do allow "mixed" RAS II/RAS III excitations into the CI space?
        If FALSE, then if there are any electrons
        in RAS III, then the number of holes in RAS I cannot exceed the given
        excitation level |detci__ex_level|. !expert -*/
        options.add_bool("MIXED", true);

        /*- Do allow "mixed" excitations involving RAS IV into the CI space.
        Useful to specify a split-virtual
        CISD[TQ] computation.  If FALSE, then if there are any electrons
        in RAS IV, then the number of holes in RAS I cannot exceed the given
        excitation level |detci__ex_level|.  !expert -*/
        options.add_bool("MIXED4", true);

        /*- Do restrict strings with $e-$ in RAS IV?  Useful to reduce the number
        of strings required if MIXED4=true, as in a split-virutal CISD[TQ]
        computation.  If more than one electron is in RAS IV, then the
        holes in RAS I cannot exceed the number of particles in
        RAS III + RAS IV (i.e., |detci__ex_level|), or else the string is discarded.
        !expert -*/
        options.add_bool("R4S", false);

        /*- SUBSECTION Diagonalization Methods -*/

        /*- This specifies which method is to be used in diagonalizing the
        Hamiltonian.  The valid options are: ``RSP``, to form the entire H
        matrix and diagonalize using libciomr to obtain all eigenvalues
        (n.b. requires HUGE memory); ``OLSEN``, to use Olsen's preconditioned
        inverse subspace method (1990); ``MITRUSHENKOV``, to use a 2x2
        Olsen/Davidson method; and ``DAVIDSON`` (or ``SEM``) to use Liu's
        Simultaneous Expansion Method, which is identical to the Davidson method
        if only one root is to be found.
        The ``SEM`` method is the most robust, but it also
        requires $2NM+1$ CI vectors on disk, where $N$ is the maximum number of
        iterations and $M$ is the number of roots. -*/
        options.add_str("DIAG_METHOD", "SEM", "RSP DAVIDSON SEM");

        /*- This specifies the type of preconditioner to use in the selected
        diagonalization method.  The valid options are: ``DAVIDSON`` which
        approximates the Hamiltonian matrix by the diagonal elements;
        ``H0BLOCK_INV`` which uses an exact Hamiltonian of |detci__h0_blocksize| and
        explicitly inverts it; ``GEN_DAVIDSON`` which does a spectral
        decomposition of H0BLOCK; ``ITER_INV`` using an iterative approach
        to obtain the correction vector of H0BLOCK.  The ``H0BLOCK_INV``, ``GEN_DAVIDSON``,
        and ``ITER_INV`` approaches are all formally equivalent but the ``ITER_INV`` is
        less computationally expensive.  Default is ``DAVIDSON``. -*/
        options.add_str("PRECONDITIONER", "DAVIDSON", "LANCZOS DAVIDSON GEN_DAVIDSON H0BLOCK ITER_INV EVANGELISTI");
        // options.add_str("PRECONDITIONER", "DAVIDSON", "LANCZOS DAVIDSON GEN_DAVIDSON H0BLOCK H0BLOCK_INV ITER_INV
        // H0BLOCK_COUPLING EVANGELISTI"); // Failures

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
        SEM method if |detci__guess_vector| is ``H0_BLOCK``.  Defaults to 1000.
        Note that the program may change the given size for Ms=0 cases
        (|detci__ms0| is TRUE) if it determines that the H0 block includes only
        one member of a pair of determinants related by time reversal symmetry.
        For very small block sizes, this could conceivably eliminate the entire
        H0 block; the program should print warnings if this occurs. !expert -*/
        options.add_int("H0_BLOCKSIZE", 1000);

        /*- size of H0 block for initial guess !expert -*/
        options.add_int("H0_GUESS_SIZE", 1000);

        /*- Do use coupling block in preconditioner? !expert -*/
        options.add_bool("H0_BLOCK_COUPLING", false);

        /*- Parameters which specifies the size of the coupling block
        within the generalized davidson preconditioner. !expert -*/
        options.add_int("H0_BLOCK_COUPLING_SIZE", 0);

        /*- Do use least-squares extrapolation in iterative solution of CI
        vector? -*/
        options.add_bool("LSE", false);

        /*- Number of iterations between least-squares extrapolations -*/
        options.add_int("LSE_COLLAPSE", 3);

        /*- Minimum converged energy for least-squares
        extrapolation to be performed -*/
        options.add_double("LSE_TOLERANCE", 3);

        /*- SUBSECTION Density Matrices -*/

        /*- Do compute one-particle density matrix if not otherwise required? -*/
        options.add_bool("OPDM", false);

        /*- Do compute two-particle density matrix if not otherwise required?
            Warning: This will hold 4 dense active TPDM's in memory !expert -*/
        options.add_bool("TPDM", false);

        /*- Do compute the transition density?  Note: only transition densities
        between roots of the same symmetry will be evaluated.  DETCI
        does not compute states of different irreps within the same
        computation; to do this, lower the symmetry of the computation.-*/
        options.add_bool("TDM", false);

        /*- Do compute the dipole moment? -*/
        options.add_bool("DIPMOM", false);

        /*- Do compute natural orbitals? -*/
        options.add_bool("NAT_ORBS", false);

        /*- SUBSECTION Root Following -*/

        /*- The root to write out the two-particle density matrix for
        (the one-particle density matrices are written for all roots).
        Useful for a state-specific CASSCF or CI optimization on an
        excited state. -*/
        options.add_int("FOLLOW_ROOT", 0);

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

        /*- What file do we start at for hd/c/s/d CIvects? Should be 350 for normal
        CI calculations and 354 if we are going to do a second monomer. !expert -*/
        options.add_int("CI_FILE_START", 350);

        /*- Guess vector type.  Accepted values are ``UNIT`` for a unit vector
        guess (|detci__num_roots| and |detci__num_init_vecs| must both be 1); ``H0_BLOCK`` to use
        eigenvectors from the H0 BLOCK submatrix (default); ``DFILE`` to use
        NUM_ROOTS previously converged vectors in the D file; !expert -*/
        options.add_str("GUESS_VECTOR", "H0_BLOCK", "UNIT H0_BLOCK DFILE");

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
        disk. -*/
        options.add_bool("RESTART", false);

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
        helpful for saving disk space.  Defaults to |detci__ci_maxiter| * |detci__num_roots|
        + |detci__num_init_vecs|. -*/
        options.add_int("MAX_NUM_VECS", 0);

        /*- Gives the number of vectors to retain when the Davidson subspace is
        collapsed (see |detci__max_num_vecs|).  If greater than one, the
        collapsed subspace retains the best estimate of the CI vector for
        the previous n iterations.   Defaults to 1. -*/
        options.add_int("COLLAPSE_SIZE", 1);

        /*- Do compute the diagonal elements of the Hamiltonian matrix
        on-the-fly? Otherwise, a diagonal element vector is written
        to a separate file on disk. !expert -*/
        options.add_bool("HD_OTF", true);

        /*- Do use the last vector space in the BVEC file to write
        scratch DVEC rather than using a separate DVEC file? (Only
        possible if |detci__num_roots| = 1.) !expert -*/
        options.add_bool("NO_DFILE", false);

        /*- SUBSECTION General-Order Perturbation Theory -*/

        /*- Do compute the MPn series out to
        kth order where k is determined by |detci__max_num_vecs| ?  For open-shell systems
        (|detci__reference| is ROHF, |detci__wfn| is ZAPTN), DETCI will compute the ZAPTn series.
        |detci__guess_vector| must be set to UNIT, |detci__hd_otf| must be set to TRUE, and
        |detci__hd_avg| must be set to orb_ener; these should happen by default for
        MPN = TRUE. -*/
        options.add_bool("MPN", false);

        /*- If 0, save the MPn energy; if 1, save the MP(2n-1) energy (if
        available from |detci__mpn_wigner| = true); if 2, save the MP(2n-2) energy (if
        available from |detci__mpn_wigner| = true). !expert -*/
        options.add_int("MPN_ORDER_SAVE", 0);

        /*- Do employ an orthonormal vector space rather than
          storing the kth order wavefunction? !expert -*/
        options.add_bool("MPN_SCHMIDT", false);

        /*- Do use Wigner formulas in the $E_{text{mp}n}$ series? !expert -*/
        options.add_bool("MPN_WIGNER", true);

        /*- The magnitude of perturbation $z$ in $H = H@@0 + z H@@1$ !expert -*/
        options.add_double("PERTURB_MAGNITUDE", 1.0);

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
        options.add_int("NUM_AMPS_PRINT", 10);

        /*- maximum number of alpha electrons in RAS III, for CC -*/
        options.add_int("CC_A_RAS3_MAX", -1);

        /*- maximum number of beta electrons in RAS III, for CC -*/
        options.add_int("CC_B_RAS3_MAX", -1);

        /*- maximum number of electrons in RAS III, for CC -*/
        options.add_int("CC_RAS3_MAX", -1);

        /*- maximum number of electrons in RAS IV, for CC -*/
        options.add_int("CC_RAS4_MAX", -1);

        /*- maximum number of electrons in RAS III + IV, for CC -*/
        options.add_int("CC_RAS34_MAX", -1);

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
        Optional additional restrictions on allowed excitations in
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
        options.add_bool("FCI_STRINGS", false);

        /*- Do string replacements on the fly in DETCI? Can
        save a gigantic amount of memory (especially for truncated CI's) but
        is somewhat flaky and hasn't been tested for a while.  It may work
        only works for certain classes of RAS calculations.  The current
        code is very slow with this option turned on. !expert -*/
        options.add_bool("REPL_OTF", false);

        /*- Do use some routines based on the papers of Bendazzoli et al.
        to calculate sigma?  Seems to be slower and not worthwhile; may disappear
        eventually.  Works only for full CI and I don't remember if I could see
        how their clever scheme might be extended to RAS in general. !expert -*/
        options.add_bool("BENDAZZOLI", false);

        /*- SUBSECTION MCSCF -*/

        /*- Convergence criterion for the RMS of the orbital gradient -*/
        options.add_double("MCSCF_R_CONVERGENCE", 1e-5);

        /*- Convergence criterion for energy. See Table :ref:`Post-SCF
        Convergence <table:conv_corl>` for default convergence criteria for
        different calculation types. -*/
        options.add_double("MCSCF_E_CONVERGENCE", 1e-7);

        /*- Maximum number MCSCF of iterations -*/
        options.add_int("MCSCF_MAXITER", 30);

        /*- Maximum value in the rotation matrix. If a value is greater than this number
        all values are scaled. -*/
        options.add_double("MCSCF_MAX_ROT", 0.5);

        /*- Method to handle the two-electron integrals -*/
        options.add_str("MCSCF_TYPE", "CONV", "DF CONV AO");

        /*- Initial MCSCF starting guess, MP2 natural orbitals only available for DF-RHF reference -*/
        options.add_str("MCSCF_GUESS", "SCF", "MP2 SCF");

        /*- Apply a list of 2x2 rotation matrices to the orbitals in the form of
        [irrep, orbital1, orbital2, theta] where an angle of 0 would do nothing and an angle
        of 90 would switch the two orbitals. -*/
        options.add("MCSCF_ROTATE", new ArrayType());

        /*- Convergence algorithm to utilize. Two-Step, Augmented Hessian. Defaults
        to TS for RASSCF. -*/
        options.add_str("MCSCF_ALGORITHM", "TS", "TS AH");

        /*- Start second-order (AH or OS) orbital-orbital MCSCF based on RMS of orbital gradient -*/
        options.add_double("MCSCF_SO_START_GRAD", 1e-4);

        /*- Start second-order (AH or OS) orbital-orbital MCSCF based on energy convergence -*/
        options.add_double("MCSCF_SO_START_E", 1e-4);

        /*- Iteration to turn on DIIS for TS convergence -*/
        options.add_int("MCSCF_DIIS_START", 3);

        /*- How often to do a DIIS extrapolation for TS convergence -*/
        options.add_int("MCSCF_DIIS_FREQ", 1);

        /*- Maximum number of DIIS vectors for TS convergence -*/
        options.add_int("MCSCF_DIIS_MAX_VECS", 8);

        /*- DIIS error vector type either, the AO orbital gradient or the orbital rotation update matrix -*/
        options.add_str("MCSCF_DIIS_ERROR_TYPE", "GRAD", "GRAD UPDATE");

        /*- Auxiliary basis set for MCSCF density fitted ERI computations.
        This only effects the "Q" matrix in Helgaker's language.
        :ref:`Defaults <apdx:basisFamily>` to a JKFIT basis. -*/
        options.add_str("DF_BASIS_MCSCF", "");

        /*- Cleanup the CI info at the end of a run? -*/
        options.add_bool("MCSCF_CI_CLEANUP", true);

        /*- Cleanup the DPD MCSCF object at the end of a run? -*/
        options.add_bool("MCSCF_DPD_CLEANUP", true);
    }

    if (name == "SAPT" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs symmetry adapted perturbation theory (SAPT)
        analysis to quantitatively analyze non-covalent interactions. -*/

        /*- SUBSECTION SAPT(HF) -*/

        /*- The level of theory for SAPT -*/
        options.add_str("SAPT_LEVEL", "SAPT0", "SAPT0 SAPT2 SAPT2+ SAPT2+3");

        /*- Whether or not to perform exchange scaling for SAPT exchange components.
        Default is false, i.e. no scaling. If set to true, performs scaling with
        $Exch10 / Exch10(S^2)$. If set to a value $\alpha$, performs scaling with
        $(Exch10 / Exch10(S^2))^{\alpha}$. -*/
        options.add_str("EXCH_SCALE_ALPHA", "FALSE", "");
        /*- For SAPT0 only, compute only first-order electrostatics and exchange.
        The integrals are computed before any terms, so all integrals will
        be computed even if they are not needed for the requested term !expert -*/
        options.add_bool("SAPT0_E10", false);
        /*- For SAPT0 only, compute only second-order induction
        The integrals are computed before any terms, so all integrals will
        be computed even if they are not needed for the requested term !expert -*/
        options.add_bool("SAPT0_E20IND", false);
        /*- For SAPT0 only, compute only second-order induction
        The integrals are computed before any terms, so all integrals will
        be computed even if they are not needed for the requested term !expert -*/
        options.add_bool("SAPT0_E20DISP", false);

        /*- Convergence criterion for residual of the CPHF/CPKS coefficients
          in the SAPT $E@@{ind,resp}^{(20)}$ term. This applies to
          wavefunction-based SAPT or SAPT(DFT). See |fisapt__cphf_r_convergence| for 
          fragment-partitioned or intramolecular SAPT. -*/
        options.add_double("CPHF_R_CONVERGENCE", 1e-8);

        /*- Solve the CPHF equations to compute coupled induction and
            exchange-induction. These are not available for ROHF, and
            the option is automatically false in this case. In all other cases,
            coupled induction is strongly recommended. Only turn it off if the
            induction energy is not going to be used.
            !expert -*/
        options.add_bool("COUPLED_INDUCTION", true);

        /*- For SAPT0 or SAPT(DFT), compute the non-approximated second-order exchange-induction term. !expert -*/
        options.add_bool("DO_IND_EXCH_SINF", false);

        /*- For SAPT0 or SAPT(DFT), compute the non-approximated second-order exchange-dispersion term. !expert -*/
        options.add_bool("DO_DISP_EXCH_SINF", false);

        /*- For SAPT2+3, compute the non-approximated third-order exchange-induction term. !expert -*/
        options.add_bool("DO_IND30_EXCH_SINF", false);

        /*- Do use asynchronous disk I/O in the solution of the CPHF equations?
        Use may speed up the computation slightly at the cost of spawning an
        additional thread. -*/
        options.add_bool("AIO_CPHF", false);

        /*- Do use asynchronous disk I/O in the formation of the DF integrals?
        Use may speed up the computation slightly at the cost of spawning an
        additional thread. -*/
        options.add_bool("AIO_DF_INTS", false);

        /*- Maximum number of CPHF iterations -*/
        options.add_int("MAXITER", 50);
        /*- Do CCD dispersion correction in SAPT2+, SAPT2+(3) or SAPT2+3? !expert -*/
        options.add_bool("DO_CCD_DISP", false);
        /*- Do MBPT dispersion correction in SAPT2+, SAPT2+(3) or SAPT2+3, if also doing CCD? !expert -*/
        options.add_bool("DO_MBPT_DISP", true);
        /*- E converge value for CCD -*/
        options.add_double("CCD_E_CONVERGENCE", 1E-8);
        /*- Convergence tolerance for CCD amplitudes -*/
        options.add_double("CCD_T_CONVERGENCE", 1E-8);
        /*- Maximum number of vectors used in CCD-DIIS -*/
        options.add_int("MAX_CCD_DIISVECS", 10);
        /*- Minimum number of vectors used in CCD-DIIS -*/
        options.add_int("MIN_CCD_DIISVECS", 4);
        /*- Max CCD iterations -*/
        options.add_int("CCD_MAXITER", 50);
        /*- Do compute third-order corrections? !expert -*/
        options.add_bool("DO_THIRD_ORDER", false);
        /*- Do natural orbitals to speed up evaluation of the triples
        contribution to dispersion by truncating the virtual orbital space?
        Recommended true for all SAPT computations. -*/
        options.add_bool("NAT_ORBS_T3", true);
        /*- Do use MP2 natural orbital approximations for the $v^4$ block of
        two-electron integrals in the evaluation of second-order T2 amplitudes?
        Recommended true for all SAPT computations. -*/
        options.add_bool("NAT_ORBS_T2", true);
        /*- Do use MP2 natural orbital approximations for the $v^4$ block of
        two-electron integrals in the evaluation of CCD T2 amplitudes?
        Recommended true for all SAPT computations. -*/
        options.add_bool("NAT_ORBS_V4", true);

        /*- Minimum occupation (eigenvalues of the MP2 OPDM) below which virtual
        natural orbitals are discarded for in each of the above three truncations
        -*/
        options.add_double("OCC_TOLERANCE", 1.0E-6);
        /*- Schwarz screening threshold.
        Minimum absolute value below which all three-index DF integrals
        and those contributing to four-index integrals are neglected. The
        default is conservative, but there isn't much to be gained from
        loosening it, especially for higher-order SAPT. -*/
        options.add_double("INTS_TOLERANCE", 1.0E-12);
        /*- Memory safety -*/
        options.add_double("SAPT_MEM_SAFETY", 0.9);
        /*- Do force SAPT2 and higher to die if it thinks there isn't enough
        memory?  Turning this off is ill-advised. -*/
        options.add_bool("SAPT_MEM_CHECK", true);
        /*- Primary basis set, describes the monomer molecular orbitals -*/
        options.add_str("BASIS", "");
        /*- Auxiliary basis set for SAPT density fitting computations.
        :ref:`Defaults <apdx:basisFamily>` to a RI basis. -*/
        options.add_str("DF_BASIS_SAPT", "");
        /*- Auxiliary basis set for SAPT Elst10 and Exch10 density fitting
        computations, may be important if heavier elements are involved.
        :ref:`Defaults <apdx:basisFamily>` to a JKFIT basis.
        Previous to v1.6, defaulted to |sapt__df_basis_sapt|. See :ref:`fitting notes <sec:saptfitA>` . -*/
        options.add_str("DF_BASIS_ELST", "");
        /*- Maximum error allowed (Max error norm in Delta tensor)
        in the approximate energy denominators employed for most of the
        $E@@{disp}^{(20)}$ and $E@@{exch-disp}^{(20)}$ evaluation. -*/
        options.add_double("DENOMINATOR_DELTA", 1.0E-6);
        /*- Denominator algorithm for PT methods. Laplace transformations
        are slightly more efficient. -*/
        options.add_str("DENOMINATOR_ALGORITHM", "LAPLACE", "LAPLACE CHOLESKY");
        /*- The scope of core orbitals to freeze in evaluation of SAPT
        $E@@{disp}^{(20)}$ and $E@@{exch-disp}^{(20)}$ terms. Recommended true
        for all SAPT computations -*/
        options.add_str("FREEZE_CORE", "FALSE", "FALSE TRUE");
        /*- The amount of information to print to the output file for the sapt
        module. For 0, only the header and final results are printed. For 1,
        (recommended for large calculations) some intermediate quantities are also
        printed. -*/
        options.add_int("PRINT", 1);
        /*- Proportion of memory available for the DF-MP2 three-index integral
            buffers used to evaluate dispersion. !expert -*/
        options.add_double("SAPT_MEM_FACTOR", 0.9);

        /*- SUBSECTION SAPT(DFT) -*/

        /*- Monomer A GRAC shift in Hartree -*/
        options.add_double("SAPT_DFT_GRAC_SHIFT_A", 0.0);
        /*- Monomer B GRAC shift in Hartree -*/
        options.add_double("SAPT_DFT_GRAC_SHIFT_B", 0.0);
        /*- Compute the Delta-HF correction? -*/
        options.add_bool("SAPT_DFT_DO_DHF", true);
        /*- How is the GRAC correction determined? !expert -*/
        options.add_str("SAPT_DFT_GRAC_DETERMINATION", "INPUT", "INPUT");
        /*- Enables the hybrid xc kernel in dispersion? !expert -*/
        options.add_bool("SAPT_DFT_DO_HYBRID", true);
        /*- Scheme for approximating exchange-dispersion for SAPT-DFT.
        Previous to Nov 2022, default was ``FIXED`` with Hesselmann value.
        ``NONE`` Use unscaled ``Exch-Disp2,u`` .
        ``FIXED`` Use a fixed factor |sapt__sapt_dft_exch_disp_fixed_scale| to scale ``Exch-Disp2,u`` .
        ``DISP`` Use the ratio of ``Disp2,r`` and ``Disp2,u`` to scale ``Exch-Disp2,u`` . -*/
        options.add_str("SAPT_DFT_EXCH_DISP_SCALE_SCHEME", "FIXED", "NONE FIXED DISP");
        /*- Exch-disp scaling factor for FIXED scheme for |sapt__sapt_dft_exch_disp_scale_scheme|.
        Default value of 0.770 suggested in Y. Xie, D. G. A. Smith and C. D. Sherrill, 2022 (submitted).
        Previous to Nov 2022, default value was 0.686 suggested by Hesselmann and Korona, J. Chem. Phys. 141, 094107 (2014). !expert -*/
        options.add_double("SAPT_DFT_EXCH_DISP_FIXED_SCALE", 0.770);
        /*- Underlying funcitonal to use for SAPT(DFT) !expert -*/
        options.add_str("SAPT_DFT_FUNCTIONAL", "PBE0", "");
        /*- Number of points in the Legendre FDDS Dispersion time integration !expert -*/
        options.add_int("SAPT_FDDS_DISP_NUM_POINTS", 10);
        /*- Lambda shift in the space morphing for the FDDS Dispersion time integration !expert -*/
        options.add_double("SAPT_FDDS_DISP_LEG_LAMBDA", 0.3);
        /*- Minimum rho cutoff for the in the LDA response for FDDS !expert -*/
        options.add_double("SAPT_FDDS_V2_RHO_CUTOFF", 1.e-6);
        /*- Which MP2 Exch-Disp module to use? !expert -*/
        options.add_str("SAPT_DFT_MP2_DISP_ALG", "SAPT", "FISAPT SAPT");
        /*- Interior option to clean up printing !expert -*/
        options.add_bool("SAPT_QUIET", false);
    }

    if (name == "FISAPT" || options.read_globals()) {
        // ==> FISAPT Options <== //

        // => Overall Options <= //

        /*- Memory safety factor for heavy FISAPT operations !expert -*/
        options.add_double("FISAPT_MEM_SAFETY_FACTOR", 0.9);
        /*- Convergence criterion for residual of the CPHF coefficients in the SAPT
        $E@@{ind,resp}^{(20)}$ term. -*/
        options.add_double("CPHF_R_CONVERGENCE", 1E-8);
        /*- Maximum number of iterations for CPHF -*/
        options.add_int("MAXITER", 50);
        /*- Schwarz screening threshold. Mininum absolute value below which TEI are neglected. -*/
        options.add_double("INTS_TOLERANCE", 0.0);

        // => ISAPT Zero-th Order Wavefunction Options <= //

        /*- Specification algorithm for link bonds in ISAPT -*/
        options.add_str("FISAPT_LINK_SELECTION", "AUTOMATIC", "AUTOMATIC MANUAL");
        /*- Amount of fragment charge completeness to distinguish link bonds -*/
        options.add_double("FISAPT_CHARGE_COMPLETENESS", 0.8);
        /*- Manual link bond specification [[Atom1, Atom2], ...] -*/
        options.add("FISAPT_MANUAL_LINKS", new ArrayType());
        /*- Where do sigma links go (to C, AB, or split into IHOs)? -*/
        options.add_str("FISAPT_LINK_ASSIGNMENT", "C", "C AB SAO0 SAO1 SAO2 SIAO0 SIAO1 SIAO2");
        /*- Orthogonalization of link orbitals for FISAPT_LINK_ASSIGNMENT=SAOx/SIAOx 
            Link A orthogonalized to A in whole (interacting) molecule or in the (noninteracting) fragment? -*/
        options.add_str("FISAPT_LINK_ORTHO", "FRAGMENT", "FRAGMENT WHOLE NONE");
        /*- Calculate separate exchange corrections for parallel and perpendicular spin coupling of link orbitals? 
            When false, only the averaged out exchange corrections are computed. -*/
        options.add_bool("FISAPT_EXCH_PARPERP", false);
        /*- Generate cube files for unsplit link orbitals (IBOs)? -*/
        options.add_bool("FISAPT_CUBE_LINKIBOS", false);
        /*- Generate cube files for split link orbitals (IHOs)? -*/
        options.add_bool("FISAPT_CUBE_LINKIHOS", false);
        /*- Generate cube files for fragment density matrices? -*/
        options.add_bool("FISAPT_CUBE_DENSMAT", false);

        // => F-SAPT Options <= //

        /*- Do an F-SAPT analysis? -*/
        options.add_bool("FISAPT_DO_FSAPT", true);
        /*- Do F-SAPT Dispersion? -*/
        options.add_bool("FISAPT_DO_FSAPT_DISP", true);
        /*- Filepath to drop F-SAPT data within input file directory -*/
        options.add_str_i("FISAPT_FSAPT_FILEPATH", "fsapt/");
        /*- Do F-SAPT exchange scaling? (ratio of S^\infty to S^2) -*/
        options.add_bool("FISAPT_FSAPT_EXCH_SCALE", true);
        /*- Do F-SAPT induction scaling? (ratio of HF induction to F-SAPT induction) -*/
        options.add_bool("FISAPT_FSAPT_IND_SCALE", true);
        /*- Do F-SAPT coupled response? (not recommended) -*/
        options.add_bool("FISAPT_FSAPT_IND_RESPONSE", false);
        /*- Do sSAPT0 exchange-scaling with F-SAPT -*/
        options.add_bool("SSAPT0_SCALE", false);
        /*- Filepath to drop sSAPT0 exchange-scaling F-SAPT data within input file directory -*/
        options.add_str_i("FISAPT_FSSAPT_FILEPATH", "s-fsapt/");

        // => CubicScalarGrid options <= //

        /*- CubicScalarGrid spatial extent in bohr [O_X, O_Y, O_Z]. Defaults to 4.0 bohr each. -*/
        options.add("CUBIC_GRID_OVERAGE", new ArrayType());
        /*- CubicScalarGrid grid spacing in bohr [D_X, D_Y, D_Z]. Defaults to 0.2 bohr each. -*/
        options.add("CUBIC_GRID_SPACING", new ArrayType());
        /*- CubicScalarGrid basis cutoff. !expert -*/
        options.add_double("CUBIC_BASIS_TOLERANCE", 1.0E-12);
        /*- CubicScalarGrid maximum number of grid points per evaluation block. !expert -*/
        options.add_int("CUBIC_BLOCK_MAX_POINTS", 1000);

        // => Scalar Field Plotting Options <= //

        /*- Plot a scalar-field analysis -*/
        options.add_bool("FISAPT_DO_PLOT", false);
        /*- Filepath to drop scalar data within input file directory -*/
        options.add_str_i("FISAPT_PLOT_FILEPATH", "plot/");

        // => Localization Tech <= //

        /*- Relative convergence in orbital localization -*/
        options.add_double("LOCAL_CONVERGENCE", 1.0E-12);
        /*- Maximum iterations in localization -*/
        options.add_int("LOCAL_MAXITER", 1000);
        /*- Use ghost atoms in Pipek-Mezey or IBO metric !expert -*/
        options.add_bool("LOCAL_USE_GHOSTS", false);
        /*- Condition number to use in IBO metric inversions !expert -*/
        options.add_double("LOCAL_IBO_CONDITION", 1.0E-7);
        /*- IBO localization metric power -*/
        options.add_int("LOCAL_IBO_POWER", 4);
        /*- MinAO Basis for IBO !expert -*/
        options.add_str("MINAO_BASIS", "CC-PVTZ-MINAO");
        /*- IBO Stars procedure -*/
        options.add_bool("LOCAL_IBO_USE_STARS", false);
        /*- IBO Charge metric for classification as Pi -*/
        options.add_double("LOCAL_IBO_STARS_COMPLETENESS", 0.90);
        /*- IBO Centers for Pi Degeneracy -*/
        options.add("LOCAL_IBO_STARS", new ArrayType());
    }
    if (name == "DCT" || options.read_globals()) {
        /*-MODULEDESCRIPTION Performs density cumulant (functional) theory
        computations -*/

        /*- Reference wavefunction type -*/
        options.add_str("REFERENCE", "RHF", "UHF RHF ROHF");
        /*- Algorithm to use for the density cumulant and orbital updates in the DCT energy computation.
        Two-step algorithm is usually more efficient for small
        systems, but for large systems simultaneous algorithm (default) is recommended.
        If convergence problems are encountered (especially
        for highly symmetric systems) QC algorithm can be used. -*/
        options.add_str("ALGORITHM", "SIMULTANEOUS", "TWOSTEP SIMULTANEOUS QC");
        /*- Algorithm to use for the solution of DC-06 response equations in computation of analytic gradients and
         * properties-*/
        options.add_str("RESPONSE_ALGORITHM", "TWOSTEP", "TWOSTEP SIMULTANEOUS");
        /*- Controls the type of the quadratically-convergent algorithm (effective for ALGORITHM = QC).
        If set to TWOSTEP the Newton-Raphson equations are only solved for the orbital updates,
        the cumulant is updated using the standard Jacobi algorithm. If set to SIMULTANEOUS both cumulant
        and orbitals are updated in a single Newton-Raphson step. -*/
        options.add_str("QC_TYPE", "SIMULTANEOUS", "TWOSTEP SIMULTANEOUS");
        /*- Convergence criterion for the RMS of the residual vector in density cumulant updates, as well as
        the solution of the density cumulant and orbital response equations. In the orbital updates controls
        the RMS of the SCF error vector -*/
        options.add_double("R_CONVERGENCE", 1e-10);
        /*- Convergence criterion for energy. See Table :ref:`Post-SCF
        Convergence <table:conv_corl>` for default convergence criteria for
        different calculation types. -*/
        options.add_double("E_CONVERGENCE", 1e-10);
        /*- Convergence criterion for the density cumulant and orbital guess for the
        variationally orbital-optimized DFT methods. Currently only available for ALGORITHM = SIMULTANEOUS. -*/
        options.add_double("GUESS_R_CONVERGENCE", 1e-3);
        /*- Maximum number of macro- or micro-iterations for both energy and response equations -*/
        options.add_int("MAXITER", 40);
        /*- Value of RMS of the density cumulant residual and SCF error vector below which DIIS extrapolation starts.
        Same keyword controls the DIIS extrapolation for the solution of the response equations. -*/
        options.add_double("DIIS_START_CONVERGENCE", 1e-3);
        /*- Maximum number of error vectors stored for DIIS extrapolation !expert-*/
        options.add_int("DIIS_MAX_VECS", 6);
        /*- Minimum number of error vectors stored for DIIS extrapolation !expert-*/
        options.add_int("DIIS_MIN_VECS", 3);
        /*- Controls whether to avoid the AO->MO transformation of the
        two-electron integrals for the four-virtual case ($\langle VV||
        VV \rangle$) by computing the corresponding terms in the AO
        basis. AO_BASIS = DISK algorithm reduces the memory requirements
        and can significantly reduce the cost of the energy computation
        if SIMULTANEOUS algorithm is used. For the TWOSTEP algorithm,
        however, AO_BASIS = DISK option is not recommended due to extra
        I/O. -*/
        options.add_str("AO_BASIS", "DISK", "NONE DISK");
        /*- The amount (percentage) of damping to apply to the orbital update procedure:
        0 will result in a full update, 100 will completely stall the
        update. A value around 20 (which corresponds to 20\% of the previous
        iteration's density being mixed into the current iteration)
        can help in cases where oscillatory convergence is observed. !expert-*/
        options.add_double("DAMPING_PERCENTAGE", 0.0);
        /*- The shift applied to the denominator in the density cumulant update iterations !expert-*/
        options.add_double("TIKHONOW_OMEGA", 0.0);
        /*- The shift applied to the denominator in the orbital update iterations !expert-*/
        options.add_double("ORBITAL_LEVEL_SHIFT", 0.0);
        /*- Controls how to cache quantities within the DPD library !expert-*/
        options.add_int("CACHELEVEL", 2);
        /*- Schwarz screening threshold. Mininum absolute value below which TEI are neglected. !expert -*/
        options.add_double("INTS_TOLERANCE", 1e-14);
        /*- Whether to read the orbitals from a previous computation, or to compute
            an MP2 guess. !expert -*/
        options.add_str("DCT_GUESS", "MP2", "CC BCC MP2 DCT");
        /*- Whether to perform a guess DC-06 or DC-12 computation for ODC-06 or ODC-12 methods, respectively.
            Currently only available for ALGORITHM = SIMULTANEOUS. -*/
        options.add_bool("ODC_GUESS", false);
        /*- Controls whether to relax the guess orbitals by taking the guess density cumulant
        and performing orbital update on the first macroiteration (for ALOGRITHM = TWOSTEP only) !expert-*/
        options.add_bool("RELAX_GUESS_ORBITALS", false);
        /*- Controls whether to include the coupling terms in the DCT electronic Hessian (for ALOGRITHM = QC
        with QC_TYPE = SIMULTANEOUS only) -*/
        options.add_bool("QC_COUPLING", false);
        /*- Performs stability analysis of the DCT energy !expert-*/
        options.add_bool("STABILITY_CHECK", false);
        /*- The value of the rms of the residual in Schmidt orthogonalization which is used as a threshold
            for augmenting the vector subspace in stability check !expert-*/
        options.add_double("STABILITY_AUGMENT_SPACE_TOL", 0.1);
        /*- Controls the convergence of the Davidson's diagonalization in stability check !expert-*/
        options.add_double("STABILITY_CONVERGENCE", 1e-4);
        /*- The number of vectors that can be added simultaneously into the subspace for Davidson's diagonalization in
           stability check !expert-*/
        options.add_int("STABILITY_ADD_VECTORS", 20);
        /*- The number of guess vectors used for Davidson's diagonalization in stability check !expert-*/
        options.add_int("STABILITY_N_GUESS_VECTORS", 20);
        /*- The number of Hessian eigenvalues computed during the stability check !expert-*/
        options.add_int("STABILITY_N_EIGENVALUES", 3);
        /*- The maximum size of the subspace for the stability check. The program will terminate if this parameter is
           exceeded and the convergence (STABILITY_CONVERGENCE) is not satisfied !expert-*/
        options.add_int("STABILITY_MAX_SPACE_SIZE", 200);
        /*- Chooses appropriate DCT method -*/
        options.add_str("DCT_FUNCTIONAL", "ODC-12", "DC-06 DC-12 ODC-06 ODC-12 ODC-13 CEPA0");
        /*- Whether to compute three-particle energy correction or not -*/
        options.add_str("THREE_PARTICLE", "NONE", "NONE PERTURBATIVE");
        /*- Level shift applied to the diagonal of the density-weighted Fock operator. While this shift can improve
           convergence, it does change the DCT energy. !expert-*/
        options.add_double("ENERGY_LEVEL_SHIFT", 0.0);
        /*- What algorithm to use for the DCT computation -*/
        options.add_str("DCT_TYPE", "CONV", "CONV DF");
        /*- Auxiliary basis set for DCT density fitting computations.
        :ref:`Defaults <apdx:basisFamily>` to a RI basis. -*/
        options.add_str("DF_BASIS_DCT", "");
        /*- Compute a (relaxed) one-particle density matrix? Can be set manually. Set internally for
         property and gradient computations. -*/
        options.add_bool("OPDM", false);
    }
    if (name == "GDMA" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs distributed multipole analysis (DMA), using
        Anthony Stone's GDMA program. See :ref:`GDMA <sec:gdma>` for more details. -*/

        /*- The order of multipole expansion on each site.  Currently limited to the same
            order for all sites; for more advanced usage a user-provided GDMA data file
            should be provided. -*/
        options.add_int("GDMA_LIMIT", 2);
        /*- The radii to be used, overriding the defaults.  Specified as an array
            [ n1, r1, n2, r2, ... ] where n1,n2,n3... are atom type strings and
            r1,r2,r3 are radii in Angstrom. -*/
        options.add("GDMA_RADIUS", new ArrayType());
        /*- The origin (in Angstrom, expressed as an [x, y, z] array) about which the total multipoles
            will be computed during DMA.  Useful for determining single site expansions at an arbitrary point. -*/
        options.add("GDMA_ORIGIN", new ArrayType());
        /*- Whether to print DMA results in atomic units or SI. -*/
        options.add_str("GDMA_MULTIPOLE_UNITS", "AU SI", "AU");
        /*- The value to switch between the older standard DMA and the new grid-based approach.
            Pairs of primitives whose exponents sum is above this value will be treated using
            standard DMA.  Set to 0 to force all pairs to be treated with standard DMA. -*/
        options.add_double("GDMA_SWITCH", 4.0);
    }

    if (name == "MINTS" || options.read_globals()) {
        /*- MODULEDESCRIPTION Called at the beginning of SCF computations,
        whenever disk-based molecular integrals are required. -*/

        /*- Primary basis set. :ref:`Available basis sets <apdx:basisElement>` -*/
        options.add_str("BASIS", "");
        /*- Omega scaling for Erf and Erfc.-*/
        options.add_double("OMEGA_ERF", 0.20);
    }
    if (name == "SCF" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs self consistent field (Hartree-Fock and
        Density Functional Theory) computations.  These are the starting
        points for most computations, so this code is called in most cases. -*/

        /*- SUBSECTION General Wavefunction Info -*/

        /*- Wavefunction type !expert -*/
        options.add_str("WFN", "SCF", "SCF");
        /*- Reference wavefunction type.
        **Cfour Interface:** Keyword translates into |cfour__cfour_reference|. -*/
        options.add_str("REFERENCE", "RHF", "RHF ROHF UHF CUHF RKS UKS");
        /*- Primary basis set -*/
        options.add_str("BASIS", "");
        /*- Auxiliary basis set for SCF density fitting computations.
        :ref:`Defaults <apdx:basisFamily>` to a JKFIT basis. -*/
        options.add_str("DF_BASIS_SCF", "");
        /*- Maximum numbers of batches to read PK supermatrix. !expert -*/
        options.add_int("PK_MAX_BUCKETS", 500);
        /*- All densities are considered non symmetric, debug only. !expert -*/
        options.add_bool("PK_ALL_NONSYM", false);
        /*- Max memory per buf for PK algo REORDER, for debug and tuning -*/
        options.add_int("MAX_MEM_BUF", 0);
        /*- Tolerance for Cholesky decomposition of the ERI tensor -*/
        options.add_double("CHOLESKY_TOLERANCE", 1e-4);
        /*- Do a density fitting SCF calculation to converge the
            orbitals before switching to the use of exact integrals in
            a |globals__scf_type| ``DIRECT`` calculation -*/
        options.add_bool("DF_SCF_GUESS", true);
        /*- For certain |globals__scf_type| algorithms that have internal sub-algorithms
            depending on available memory or other hardware constraints, allow the best
            sub-algorithm for the molecule and conditions (``AUTO`` ; usual mode) or
            forcibly select a sub-algorithm (usually only for debugging or profiling).
            Presently, ``SCF_SUBTYPE=DF``, ``SCF_SUBTYPE=MEM_DF``, and ``SCF_SUBTYPE=DISK_DF`` 
	        can have ``INCORE`` and ``OUT_OF_CORE`` selected; and ``SCF_TYPE=PK``  can have ``INCORE``,
	        ``OUT_OF_CORE``, ``YOSHIMINE_OUT_OF_CORE``, and ``REORDER_OUT_OF_CORE`` selected. !expert -*/
	    options.add_str("SCF_SUBTYPE", "AUTO", "AUTO INCORE OUT_OF_CORE YOSHIMINE_OUT_OF_CORE REORDER_OUT_OF_CORE");
        /*- Keep JK object for later use? -*/
        options.add_bool("SAVE_JK", false);
        /*- Memory safety factor for allocating JK -*/
        options.add_double("SCF_MEM_SAFETY_FACTOR", 0.75);
        /*- SO orthogonalization: automatic, symmetric, or canonical? -*/
        options.add_str("S_ORTHOGONALIZATION", "AUTO", "AUTO SYMMETRIC CANONICAL PARTIALCHOLESKY");
        /*- Minimum S matrix eigenvalue to allow before linear dependencies are removed. -*/
        options.add_double("S_TOLERANCE", 1E-7);
        /*- Tolerance for partial Cholesky decomposition of overlap matrix. -*/
        options.add_double("S_CHOLESKY_TOLERANCE", 1E-8);
        /*- Screening threshold for the chosen screening method (SCHWARZ, CSAM, DENSITY)
          Absolute value below which TEI are neglected. -*/
        options.add_double("INTS_TOLERANCE", 1E-12);
        /*- The type of guess orbitals. See :ref:`sec:scfguess` for what the options mean and
         what the defaults are. -*/
        options.add_str("GUESS", "AUTO", "AUTO CORE GWH SAD SADNO SAP SAPGAU HUCKEL MODHUCKEL READ");
        /*- The potential basis set used for the SAPGAU guess -*/
        options.add_str("SAPGAU_BASIS", "sap_helfem_large");
        /*- Mix the HOMO/LUMO in UHF or UKS to break alpha/beta spatial symmetry.
        Useful to produce broken-symmetry unrestricted solutions.
        Notice that this procedure is defined only for calculations in C1 symmetry. -*/
        options.add_bool("GUESS_MIX", false);
        /*- Do write a MOLDEN output file?  If so, the filename will end in
        .molden, and the prefix is determined by |globals__writer_file_label|
        (if set), or else by the name of the output file plus the name of
        the current molecule. -*/
        options.add_bool("MOLDEN_WRITE", false);
        /*- If true, then repeat the specified guess procedure for the orbitals every time -
        even during a geometry optimization. -*/
        options.add_bool("GUESS_PERSIST", false);
        /*- File name (case sensitive) to which to serialize Wavefunction orbital data. -*/
        options.add_str_i("ORBITALS_WRITE", "");

        /*- Do print the molecular orbitals? -*/
        options.add_bool("PRINT_MOS", false);
        /*- Do print the basis set? -*/
        options.add_bool("PRINT_BASIS", false);
        /*- Do perform a QCHF computation?  -*/
        options.add_bool("QCHF", false);

        /*- SCF Properties to calculate after an energy evaluation. Note, this
        keyword is not used for property evaluations. -*/
        options.add("SCF_PROPERTIES", new ArrayType());

        /*- SUBSECTION Convergence Control/Stabilization -*/

        /*- Maximum number of iterations.
        **Cfour Interface:** Keyword translates into |cfour__cfour_scf_maxcyc|. -*/
        options.add_int("MAXITER", 100);
        /*- Fail if we reach maxiter without converging? -*/
        options.add_bool("FAIL_ON_MAXITER", true);
        /*- Convergence criterion for SCF energy. See Table :ref:`SCF
        Convergence & Algorithm <table:conv_scf>` for default convergence
        criteria for different calculation types. -*/
        options.add_double("E_CONVERGENCE", 1e-6);
        /*- Convergence criterion for SCF density, defined as the RMS
        or maximum absolute value of the orbital gradient.  See Table
        :ref:`SCF Convergence & Algorithm <table:conv_scf>` for
        default convergence criteria for different calculation types.
        **Cfour Interface:** Keyword translates into
        |cfour__cfour_scf_conv|. -*/
        options.add_double("D_CONVERGENCE", 1e-6);
        /*- The amount (percentage) of damping to apply to the early density updates.
            0 will result in a full update, 100 will completely stall the update.  A
            value around 20 (which corresponds to 20\% of the previous iteration's
            density being mixed into the current density)
            could help to solve problems with oscillatory convergence. -*/
        options.add_double("DAMPING_PERCENTAGE", 0.0);
        /*- The density convergence threshold after which damping is no longer performed, if it is enabled.
            It is recommended to leave damping on until convergence, which is the default.
        **Cfour Interface:** Keyword translates into |cfour__cfour_scf_damping|. -*/
        options.add_double("DAMPING_CONVERGENCE", 1.0E-18);
        /*- Accelerate convergence by performing a preliminary SCF with
        this small basis set followed by projection into the full target
        basis. A value of ``TRUE`` turns on projection using the
        :ref:`Defaults <apdx:basisFamily>` small basis set 3-21G, pcseg-0, or def2-SV(P). -*/
        options.add_str("BASIS_GUESS", "FALSE", "");
        /*- When |scf__basis_guess| is active, run the preliminary scf in
        density-fitted mode with this as fitting basis for the small basis
        set. A value of ``TRUE`` turns on density fitting with the
        default basis, otherwise the specified basis is used. -*/
        options.add_str("DF_BASIS_GUESS", "FALSE", "");
        /*- Use RMS error instead of the more robust absolute error? -*/
        options.add_bool("DIIS_RMS_ERROR", true);
        /*- The minimum iteration to start storing DIIS vectors and performing ADIIS/EDIIS. -*/
        options.add_int("DIIS_START", 1);
        /*- Minimum number of error vectors stored for DIIS extrapolation. Will be removed in v1.7. -*/
        options.add_int("DIIS_MIN_VECS", 2);
        /*- Maximum number of error vectors stored for DIIS extrapolation -*/
        options.add_int("DIIS_MAX_VECS", 10);
        /*- Do use DIIS extrapolation to accelerate convergence? -*/
        options.add_bool("DIIS", true);
        /*- Do use a level shift? -*/
        options.add_double("LEVEL_SHIFT", 0.0);
        /*- DIIS error at which to stop applying the level shift -*/
        options.add_double("LEVEL_SHIFT_CUTOFF", 1e-2);
        /*- The iteration to start MOM on (or 0 for no MOM) -*/
        options.add_int("MOM_START", 0);
        /*- The absolute indices of orbitals to excite from in MOM (+/- for alpha/beta) -*/
        options.add("MOM_OCC", new ArrayType());
        /*- The absolute indices of orbitals to excite to in MOM (+/- for alpha/beta) -*/
        options.add("MOM_VIR", new ArrayType());
        /*- Convergence threshold (max 2-norm) for numerical solvers (instability analysis and CPHF/CPKS). -*/
        options.add_double("SOLVER_CONVERGENCE", 1.0E-6);
        /*- Maximum iterations for numerical solvers (instability analysis and CPHF/CPKS).  -*/
        options.add_int("SOLVER_MAXITER", 100);
        /*- Number of guess vectors per root for instability analysis. -*/
        options.add_int("SOLVER_N_GUESS", 1);
        /*- Number of roots to converge for all irreps during instability analysis. (Overridden by SOLVER_ROOTS_PER_IRREP.) -*/
        options.add_int("SOLVER_N_ROOT", 1);
        /*- Number of roots to converge, per irrep, during instability analysis. (Overrides SOLVER_N_ROOT.) -*/
        options.add("SOLVER_ROOTS_PER_IRREP", new ArrayType());
        /*- Do use second-order SCF convergence methods? -*/
        options.add_bool("SOSCF", false);
        /*- When to start second-order SCF iterations based on gradient RMS. -*/
        options.add_double("SOSCF_START_CONVERGENCE", 1.0E-2);
        /*- Minimum number of second-order microiterations to perform. -*/
        options.add_int("SOSCF_MIN_ITER", 1);
        /*- Maximum number of second-order microiterations to perform. -*/
        options.add_int("SOSCF_MAX_ITER", 5);
        /*- Second order convergence threshold. Cease microiterating at this value. -*/
        options.add_double("SOSCF_CONV", 5.0E-3);
        /*- Do we print the SOSCF microiterations?. -*/
        options.add_bool("SOSCF_PRINT", false);
        /*- Whether to perform stability analysis after convergence.  NONE prevents analysis being
            performed. CHECK will print out the analysis of the wavefunction stability at the end of
            the computation.  FOLLOW will perform the analysis and, if a totally symmetric instability
            is found, will attempt to follow the eigenvector and re-run the computations to find a stable
            solution. -*/
        options.add_str("STABILITY_ANALYSIS", "NONE", "NONE CHECK FOLLOW");
        /*- When using |scf__stability_analysis| ``FOLLOW``, how much to scale the step along the eigenvector
            by. A full step of $pi/2$ corresponds to a value of 1.0. !expert -*/
        options.add_double("FOLLOW_STEP_SCALE", 0.5);
        /*- When using STABILITY_ANALYSIS = FOLLOW, the increment to modify |scf__follow_step_scale| value
            if we end up in the same SCF solution. !expert -*/
        options.add_double("FOLLOW_STEP_INCREMENT", 0.2);
        /*- When using |scf__stability_analysis| ``FOLLOW``, maximum number of orbital optimization attempts
            to make the wavefunction stable. !expert -*/
        options.add_int("MAX_ATTEMPTS", 1);
        /*- Use a method to accelerate initial SCF convergence? Use ``NONE`` for DIIS alone (if enabled) and ``EDIIS`` or ``ADIIS``
            to have both the chosen accelerator and DIIS (if enabled). For restricted-open references, ``EDIIS`` and ``ADIIS`` have no effect. -*/
        options.add_str("SCF_INITIAL_ACCELERATOR", "ADIIS", "NONE EDIIS ADIIS");
        /*- SCF error at which to start the linear interpolation between DIIS steps and steps of the initial SCF accelerator.
            Value taken from Garza and Scuseria, DOI: 10.1063/1.4740249 -*/
        options.add_double("SCF_INITIAL_START_DIIS_TRANSITION", 1.0E-1);
        /*- SCF error at which to complete the linear interpolation between DIIS steps and steps of the initial SCF accelerator
            Value taken from Garza and Scuseria, DOI: 10.1063/1.4740249 -*/
        options.add_double("SCF_INITIAL_FINISH_DIIS_TRANSITION", 1.0E-4);
        /*- Do perform incremental Fock build? -*/
        options.add_bool("INCFOCK", false);
        /*- Frequency with which to compute the full Fock matrix if using |scf__incfock| . 
        N means rebuild every N SCF iterations to avoid accumulating error from the incremental procedure. -*/
        options.add_int("INCFOCK_FULL_FOCK_EVERY", 5);
        /*- The density threshold at which to stop building the Fock matrix incrementally -*/
        options.add_double("INCFOCK_CONVERGENCE", 1.0e-5);

        /*- The screening tolerance used for ERI/Density sparsity in the LinK algorithm -*/
        options.add_double("LINK_INTS_TOLERANCE", 1.0e-12);

        /*- SUBSECTION Fractional Occupation UHF/UKS -*/

        /*- The iteration to start fractionally occupying orbitals (or 0 for no fractional occupation) -*/
        options.add_int("FRAC_START", 0);
        /*- The absolute indices of occupied orbitals to fractionally occupy (+/- for alpha/beta) -*/
        options.add("FRAC_OCC", new ArrayType());
        /*- The occupations of the orbital indices specified above ($0.0\le {\rm occ} \le 1.0$) -*/
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
        /*- Size of the perturbation (applies only to dipole perturbations).  Deprecated - use PERTURB_DIPOLE instead
           -*/
        options.add_double("PERTURB_MAGNITUDE", 0.0);
        /*- An array of length three describing the magnitude (atomic units) of the dipole field in the {x,y,z}
           directions -*/
        options.add("PERTURB_DIPOLE", new ArrayType());
        /*- The operator used to perturb the Hamiltonian, if requested.  DIPOLE_X, DIPOLE_Y and DIPOLE_Z will be
            removed in favor of the DIPOLE option in the future -*/
        options.add_str("PERTURB_WITH", "DIPOLE", "DIPOLE DIPOLE_X DIPOLE_Y DIPOLE_Z EMBPOT SPHERE DX");
        /*- An ExternalPotential (built by Python or nullptr/None) -*/
        options.add_bool("EXTERN", false);

        /*- Radius (bohr) of a hard-sphere external potential -*/
        options.add_double("RADIUS", 10.0);  // bohr
        /*- Thickness (bohr) of a hard-sphere external potential -*/
        options.add_double("THICKNESS", 20.0);  // bohr
        /*- Number of radial grid points for spherical potential integration -*/
        options.add_int("R_POINTS", 100);
        /*- Number of colatitude grid points for spherical potential integration -*/
        options.add_int("THETA_POINTS", 360);
        /*- Number of azimuthal grid points for spherical potential integration -*/
        options.add_int("PHI_POINTS", 360);
        /*- Read an external potential from the .dx file? -*/
        options.add_bool("ONEPOT_GRID_READ", false);

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
        options.add_str("SAPT", "FALSE",
                        "FALSE 2-DIMER 2-MONOMER_A 2-MONOMER_B 3-TRIMER 3-DIMER_AB 3-DIMER_BC 3-DIMER_AC 3-MONOMER_A "
                        "3-MONOMER_B 3-MONOMER_C");

        /*- SUBSECTION DFSCF Algorithm -*/

        /*- Number of threads for integrals (may be turned down if memory is an issue). 0 is blank -*/
        options.add_int("DF_INTS_NUM_THREADS", 0);
        /*- IO caching for CP corrections, etc. Previous to v1.10, changing this selected Disk_DF over Mem_DF.
	That is, setting this forced DiskDFJK when SCF_TYPE=DF. Starting with v1.10, changing this affects 
 	|globals__scf_type| = ``CD`` or ``DISK_DF`` but does not force ``DISK_DF`` when given ``DF``. !expert -*/
        options.add_str("DF_INTS_IO", "NONE", "NONE SAVE LOAD");
        /*- Fitting Condition, i.e. eigenvalue threshold for RI basis. Analogous to S_TOLERANCE !expert -*/
        options.add_double("DF_FITTING_CONDITION", 1.0E-10);
        /*- FastDF Fitting Metric -*/
        options.add_str("DF_METRIC", "COULOMB", "COULOMB EWALD OVERLAP");
        /*- FastDF SR Ewald metric range separation parameter -*/
        options.add_double("DF_THETA", 1.0);
        /*- FastDF geometric fitting domain selection algorithm -*/
        options.add_str("DF_DOMAINS", "DIATOMIC", "DIATOMIC SPHERES");
        /*- Bump function min radius -*/
        options.add_double("DF_BUMP_R0", 0.0);
        /*- Bump function max radius -*/
        options.add_double("DF_BUMP_R1", 0.0);

        /*- SUBSECTION COSX Algorithm -*/

        /*- Number of spherical points in initial COSX grid. -*/
        options.add_int("COSX_SPHERICAL_POINTS_INITIAL", 50);
        /*- Number of radial points in initial COSX grid. -*/
        options.add_int("COSX_RADIAL_POINTS_INITIAL", 25);
        /*- Number of spherical points in final COSX grid. -*/
        options.add_int("COSX_SPHERICAL_POINTS_FINAL", 110);
        /*- Number of radial points in final COSX grid. -*/
        options.add_int("COSX_RADIAL_POINTS_FINAL", 35);
        /*- Screening criteria for integrals and intermediates in COSX -*/
        options.add_double("COSX_INTS_TOLERANCE", 1.0E-11);
        /*- Controls SCF iteration behavior for the larger (i.e., final) COSX grid.
        -1 fully converges the SCF on the final grid if possible, ending early if |scf__maxiter| total SCF iterations are reached (failure).
        0 disables the final COSX grid entirely.
        n runs up to n iterations on the final COSX grid, ending early if SCF convergence is reached (success) or if |scf__maxiter| total SCF iterations are reached (failure). -*/
        options.add_int("COSX_MAXITER_FINAL", 1);
        /*- Screening criteria for shell-pair densities in COSX !expert -*/
        options.add_double("COSX_DENSITY_TOLERANCE", 1.0E-10);
        /*- Screening criteria for basis function values on COSX grids !expert -*/
        options.add_double("COSX_BASIS_TOLERANCE", 1.0E-10);
        /*- Pruning scheme for COSX grids !expert -*/
        options.add_str("COSX_PRUNING_SCHEME", "ROBUST", 
                        "ROBUST TREUTLER NONE FLAT P_GAUSSIAN D_GAUSSIAN P_SLATER D_SLATER LOG_GAUSSIAN LOG_SLATER NONE");
        /*- Do reduce numerical COSX errors with overlap fitting? !expert -*/
        options.add_bool("COSX_OVERLAP_FITTING", true);

        /*- SUBSECTION snLinK Algorithm -*/

        /*- Number of spherical points in snLinK grid. -*/
        options.add_int("SNLINK_SPHERICAL_POINTS", 302);
        /*- Number of radial points in snLinK grid. -*/
        options.add_int("SNLINK_RADIAL_POINTS", 70);
        /*- Radial Scheme for snLinK grid. 
        MURA is default here as it matches the GauXC default option -*/
        options.add_str("SNLINK_RADIAL_SCHEME", "MURA", "MURA TREUTLER EM");
        /*- Use GPU for GauXC? -*/
        options.add_bool("SNLINK_USE_GPU", false);
        /*- Proportion (in %) of available GPU memory to allocate to snLinK. !expert-*/
        options.add_bool("SNLINK_GPU_MEM", 90);
        /*- Screening criteria for integrals and intermediates in snLinK -*/
        options.add_double("SNLINK_INTS_TOLERANCE", 1.0E-11);
        /*- Screening criteria for shell-pair densities in snLinK !expert -*/
        options.add_double("SNLINK_DENSITY_TOLERANCE", 1.0E-10);
        /*- Screening criteria for basis function values on snLinK grids !expert -*/
        options.add_double("SNLINK_BASIS_TOLERANCE", 1.0E-10);
        /*- Force snLinK to use cartesian coordinates !expert -*/
        options.add_bool("SNLINK_FORCE_CARTESIAN", false);
        /*- Pruning scheme for snLinK grids !expert -*/
        options.add_str("SNLINK_PRUNING_SCHEME", "ROBUST", "ROBUST TREUTLER NONE");
        /*- Maximum number of grid points per grid block for GauXC !expert -*/ 
        options.add_int("SNLINK_GRID_BATCH_SIZE", 2048);
        /*- Load Balancer kernel for snLinK !expert -*/
        options.add_str("SNLINK_LOAD_BALANCER_KERNEL", "DEFAULT", "DEFAULT REPLICATED REPLICATED-PETITE REPLICATED-FILLIN");
        /*- Molecular Weights kernel for snLinK !expert -*/
        options.add_str("SNLINK_MOL_WEIGHTS_KERNEL", "DEFAULT", "DEFAULT");
        /*- Integrator execution kernel for snLinK !expert 
        GauXC also has SHELLBATCHED, but it is incompatible with Psi4 due to not 
        being yet implemented with sn-LinK. -*/
        options.add_str("SNLINK_INTEGRATOR_KERNEL", "DEFAULT", "DEFAULT INCORE");
        /*- Integrator reduction kernel for snLinK !expert 
        GauXC also has NCCL, but it is incompatible with Psi4 due to requiring MPI. -*/
        options.add_str("SNLINK_REDUCTION_KERNEL", "DEFAULT", "DEFAULT BASICMPI");
        /*- Integrator local work driver kernel for snLinK !expert 
        GauXC also has SCHEME1-CUTLASS, but it is disabled in Psi4 for now 
        due to compile-time issues and requiring very modern CUDA CCs (>=80) -*/
        options.add_str("SNLINK_LWD_KERNEL", "DEFAULT", "DEFAULT REFERENCE SCHEME1 SCHEME1-MAGMA"); 
        /*- Overwrite sn-LinK grid options with debug grid matching GauXC's Ultrafine grid spec !expert -*/
        options.add_bool("SNLINK_USE_DEBUG_GRID", false);

        /*- SUBSECTION SAD Guess Algorithm -*/

        /*- The amount of SAD information to print to the output !expert -*/
        options.add_int("SAD_PRINT", 0);
        /*- Convergence criterion for SCF energy in the SAD guess, analogous to |scf__e_convergence|. -*/
        options.add_double("SAD_E_CONVERGENCE", 1E-5);
        /*- Convergence criterion for SCF density in the SAD guess, analogous to |scf__d_convergence|. -*/
        options.add_double("SAD_D_CONVERGENCE", 1E-5);
        /*- Density fitting basis used in SAD !expert -*/
        options.add_str("DF_BASIS_SAD", "SAD-FIT");
        /*- Maximum number of atomic SCF iterations within SAD !expert -*/
        options.add_int("SAD_MAXITER", 50);
        /*- SCF type used for atomic calculations in SAD guess !expert -*/
        options.add_str("SAD_SCF_TYPE", "DF", "DIRECT DF MEM_DF DISK_DF PK OUT_OF_CORE CD GTFOCK");
        /*- Do force an even distribution of occupations across the last partially occupied orbital shell? !expert -*/
        options.add_bool("SAD_FRAC_OCC", true);
        /*- Do use spin-averaged occupations instead of atomic ground spin state in fractional SAD? !expert -*/
        options.add_bool("SAD_SPIN_AVERAGE", true);
        /*- SAD guess density decomposition threshold !expert -*/
        options.add_double("SAD_CHOL_TOLERANCE", 1E-7);
#ifdef USING_OpenOrbitalOptimizer
        /*- Orbital optimizer package to use for SAD guess. If compiled with OpenOrbitalOptimizer support, change this option to use the internal code. -*/
        options.add_str("ORBITAL_OPTIMIZER_PACKAGE", "OPENORBITALOPTIMIZER", "INTERNAL OPENORBITALOPTIMIZER");
#else
        /*- Orbital optimizer package to use for SAD guess. -*/
        options.add_str("ORBITAL_OPTIMIZER_PACKAGE", "INTERNAL", "INTERNAL");
#endif

        /*- SUBSECTION DFT -*/

        /*- The DFT Range-separation parameter -*/
        options.add_double("DFT_OMEGA", 0.0);
        /*- The DFT Exact-exchange parameter -*/
        options.add_double("DFT_ALPHA", 0.0);
        /*- The DFT Correlation Range-separation parameter -*/
        options.add_double("DFT_OMEGA_C", 0.0);
        /*- The DFT Correlation hybrid parameter -*/
        options.add_double("DFT_ALPHA_C", 0.0);
        /*- Minima spin-summed density cutoff for the second derivative. Defaults to the density tolerance. -*/
        options.add_double("DFT_V2_RHO_CUTOFF", -1.0);
        /*- The gradient regularized asymptotic correction shift value -*/
        options.add_double("DFT_GRAC_SHIFT", 0.0);
        /*- The gradient regularized asymptotic correction alpha value -*/
        options.add_double("DFT_GRAC_ALPHA", 0.5);
        /*- The gradient regularized asymptotic correction beta value -*/
        options.add_double("DFT_GRAC_BETA", 40.0);
        /*- The gradient regularized asymptotic correction functional exch form. !expert -*/
        options.add_str("DFT_GRAC_X_FUNC", "XC_GGA_X_LB");
        /*- The gradient regularized asymptotic correction functional corr form. !expert -*/
        options.add_str("DFT_GRAC_C_FUNC", "XC_LDA_C_VWN");
        /*- Number of spherical points (A :ref:`Lebedev Points <table:lebedevorder>` number). -*/
        options.add_int("DFT_SPHERICAL_POINTS", 302);
        /*- Number of radial points. -*/
        options.add_int("DFT_RADIAL_POINTS", 75);
        /*- Spherical Scheme. -*/
        options.add_str("DFT_SPHERICAL_SCHEME", "LEBEDEV", "LEBEDEV");
        /*- Radial Scheme. -*/
        options.add_str("DFT_RADIAL_SCHEME", "TREUTLER", "TREUTLER BECKE MULTIEXP EM MURA");
        /*- Nuclear Scheme. -*/
        options.add_str("DFT_NUCLEAR_SCHEME", "TREUTLER", "TREUTLER BECKE NAIVE STRATMANN SBECKE");
        /*- Factor for effective BS radius in radial grid. -*/
        options.add_double("DFT_BS_RADIUS_ALPHA", 1.0);
        /*- DFT basis cutoff. -*/
        options.add_double("DFT_BASIS_TOLERANCE", 1.0E-12);
        /*- grid weight cutoff. Disable with -1.0. !expert -*/
        options.add_double("DFT_WEIGHTS_TOLERANCE", 1.0E-15);
        /*- density cutoff for LibXC. A negative value turns the feature off and LibXC defaults are used. !expert -*/
        options.add_double("DFT_DENSITY_TOLERANCE", -1.0);
        /*- The DFT grid specification, such as SG1.!expert -*/
        options.add_str("DFT_GRID_NAME", "", "SG0 SG1");
        /*- Select approach for pruning. Options ``ROBUST`` and ``TREUTLER`` prune based on regions (proximity to nucleus) while
        ``FLAT`` ``P_GAUSSIAN`` ``D_GAUSSIAN`` ``P_SLATER`` ``D_SLATER`` ``LOG_GAUSSIAN`` ``LOG_SLATER`` prune based on decaying functions (experts only!).
        The recommended scheme is ``ROBUST``. -*/
        options.add_str("DFT_PRUNING_SCHEME", "NONE",
                        "ROBUST TREUTLER NONE FLAT P_GAUSSIAN D_GAUSSIAN P_SLATER D_SLATER LOG_GAUSSIAN LOG_SLATER NONE");
        /*- Spread alpha for logarithmic pruning. !expert -*/
        options.add_double("DFT_PRUNING_ALPHA", 1.0);
        /*- The maximum number of grid points per evaluation block. !expert -*/
        options.add_int("DFT_BLOCK_MAX_POINTS", 256);
        /*- The minimum number of grid points per evaluation block. !expert -*/
        options.add_int("DFT_BLOCK_MIN_POINTS", 100);
        /*- The maximum radius to terminate subdivision of an octree block [au]. !expert -*/
        options.add_double("DFT_BLOCK_MAX_RADIUS", 3.0);
        /*- Remove points from the quadrature grid that exceed the spatial extend of the basis functions. !expert -*/
        options.add_bool("DFT_REMOVE_DISTANT_POINTS",true);
        /*- The blocking scheme for DFT. !expert -*/
        options.add_str("DFT_BLOCK_SCHEME", "OCTREE", "NAIVE OCTREE ATOMIC");
        /*- Parameters defining the dispersion correction. See Table
        :ref:`-D Functionals <table:dft_disp>` for default values and Table
        :ref:`Dispersion Corrections <table:dashd>` for the order in which
        parameters are to be specified in this array option.
        Unused for functionals constructed by user. -*/
        options.add("DFT_DISPERSION_PARAMETERS", new ArrayType());
        /*- Parameters defining the -NL/-V dispersion correction. First b, then C -*/
        options.add("NL_DISPERSION_PARAMETERS", new ArrayType());
        /*- Number of spherical points (A :ref:`Lebedev Points <table:lebedevorder>` number) for VV10 NL integration.
           -*/
        options.add_int("DFT_VV10_SPHERICAL_POINTS", 146);
        /*- Number of radial points for VV10 NL integration. -*/
        options.add_int("DFT_VV10_RADIAL_POINTS", 50);
        /*- Rho cutoff for VV10 NL integration. !expert -*/
        options.add_double("DFT_VV10_RHO_CUTOFF", 1.e-8);
        /*- Define VV10 parameter b -*/
        options.add_double("DFT_VV10_B", 0.0);
        /*- Define VV10 parameter C -*/
        options.add_double("DFT_VV10_C", 0.0);
        /*- post-scf VV10 correction -*/
        options.add_bool("DFT_VV10_POSTSCF", false);
        /*- The convergence on the orbital localization procedure -*/
        options.add_double("LOCAL_CONVERGENCE", 1E-12);
        /*- The maxiter on the orbital localization procedure -*/
        options.add_int("LOCAL_MAXITER", 200);
        /*- The number of NOONs to print in a UHF calc -*/
        options.add_str("UHF_NOONS", "3");
        /*- Save the UHF NOs -*/
        options.add_bool("SAVE_UHF_NOS", false);

        /*- SUBSECTION TDSCF -*/
        /*- Number of roots (excited states) we should seek to converge. This
        can be either an integer (total number of states to seek) or a list
        (number of states per irrep). The latter is only valid if the system has
        symmetry. Furthermore, the total number of states will be redistributed
        among irreps when symmetry is used.-*/
        options.add("TDSCF_STATES", new ArrayType());
        /*- Controls inclusion of triplet states, which is only valid for restricted references. Valid options:
            - none : No triplets computed (default)
            - also : lowest-energy triplets and singlets included, in 50-50
              ratio. Note that singlets are privileged, i.e. if seeking to
              converge 5 states in total, 3 will be singlets and 2 will be
              triplets.
            - only : Only triplet states computed
             -*/
        options.add_str("TDSCF_TRIPLETS", "NONE", "NONE ALSO ONLY");
        /*- Run with Tamm-Dancoff approximation (TDA), uses random-phase approximation (RPA) when false -*/
        options.add_bool("TDSCF_TDA", false);
        /*- Convergence threshold for the norm of the residual vector. If unset,
        default based on |scf__d_convergence|. -*/
        options.add_double("TDSCF_R_CONVERGENCE", 1E-4);
        /*- Guess type, only 'denominators' currently supported -*/
        options.add_str("TDSCF_GUESS", "DENOMINATORS");
        /*- Maximum number of TDSCF solver iterations -*/
        options.add_int("TDSCF_MAXITER", 60);
        /*- Verbosity level in TDSCF -*/
        options.add_int("TDSCF_PRINT", 1);
        /*- Cutoff for printing excitations and de-excitations contributing to each excited state -*/
        options.add_double("TDSCF_COEFF_CUTOFF", 0.1);
        /*- Which transition dipole moments to print out:
            - E_TDM_LEN : electric transition dipole moments, length representation
            - E_TDM_VEL : electric transition dipole moments, velocity representation
            - M_TDM : magnetic transition dipole moments -*/
        options.add("TDSCF_TDM_PRINT", new ArrayType());

        /*- combine omega exchange and Hartree--Fock exchange into
              one matrix for efficiency?
              Disabled until fixed.-*/
            //NOTE: Re-enable with below doc string:
            // Default is True for MemDFJK
            //   (itself the default for |globals__scf_type| DF),
            // False otherwise as not yet implemented. -*/
        options.add_bool("WCOMBINE", false);
    }
    if (name == "CPHF" || options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- The amount of debug information printed
            to the output file -*/
        options.add_int("DEBUG", 0);
        /*- What app to test?
          -*/
        options.add_str("MODULE", "RCPHF", "RCPHF");
        /*- Do explicit hamiltonian only? -*/
        options.add_bool("EXPLICIT_HAMILTONIAN", false);
        /*- Which tasks to run CPHF For
         *  Valid choices:
         *  -Polarizability
         * -*/
        options.add("CPHF_TASKS", new ArrayType());
        /*- Memory safety factor for allocating JK
        -*/
        options.add_double("CPHF_MEM_SAFETY_FACTOR", 0.75);
        /*- SCF Type
         -*/
        options.add_str("SCF_TYPE", "DIRECT", "DIRECT DF PK OUT_OF_CORE PS INDEPENDENT GTFOCK DFDIRJ+SNLINK DFDIRJ+LINK DFDIRJ+COSX");
        /*- Auxiliary basis for SCF
         -*/
        options.add_str("DF_BASIS_SCF", "");
        /*- Solver maximum iterations  -*/
        options.add_int("SOLVER_MAXITER", 100);
        /*- Solver convergence threshold (max 2-norm). -*/
        options.add_double("SOLVER_CONVERGENCE", 1.0E-6);
        /*- DL Solver number of guesses  -*/
        options.add_int("SOLVER_N_GUESS", 1);
        /*- Solver precondition type
         -*/
        options.add_str("SOLVER_PRECONDITION", "JACOBI", "SUBSPACE JACOBI NONE");
    }
    if (name == "CCTRANSORT" || options.read_globals()) {
        /*- MODULEDESCRIPTION Transforms and sorts integrals for CC codes. Called before (non-density-fitted) MP2 and
           coupled cluster computations. -*/
        /*- Wavefunction type !expert -*/
        options.add_str("WFN", "");
        /*- Reference wavefunction type -*/
        options.add_str("REFERENCE", "RHF");
        /*- The algorithm to use for the $\left\langle VV||VV \right\rangle$ terms -*/
        options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
        /*- Delete the SO two-electron integrals after the transformation? -*/
        options.add_bool("DELETE_TEI", true);
        /*- Caching level for libdpd -*/
        options.add_int("CACHELEVEL", 2);
        /*- Force conversion of ROHF MOs to semicanonical MOs to run UHF-based energies -*/
        options.add_bool("SEMICANONICAL", false);
        /*- Use cctransort module NOTE: Turning this option off requires separate
           installation of  ccsort and transqt2 modules, see http://github.com/psi4/psi4pasture -*/
        options.add_bool("RUN_CCTRANSORT", true);
    }
    if (name == "CCTRIPLES" || options.read_globals()) {
        /*- MODULEDESCRIPTION Computes the triples component of CCSD(T) energies (and gradients, if necessary). -*/
        /*- Wavefunction type !expert -*/
        options.add_str("WFN", "SCF");
        /*- Reference wavefunction type -*/
        options.add_str("REFERENCE", "RHF");
        /*- Number of threads -*/
        options.add_int("CC_NUM_THREADS", 1);
        /*- Convert ROHF MOs to semicanonical MOs -*/
        options.add_bool("SEMICANONICAL", true);
    }
    if (name == "CCDENSITY" || options.read_globals()) {
        /*- MODULEDESCRIPTION Computes the coupled cluster density matrices. Called whenever CC properties and/or
            gradients are required. -*/
        /*- Wavefunction type !expert -*/
        options.add_str("WFN", "SCF");
        /*- Reference wavefunction type -*/
        options.add_str("REFERENCE", "RHF");
        /*- Schwarz screening threshold. Mininum absolute value below which TEI are neglected. -*/
        options.add_double("INTS_TOLERANCE", 1e-14);
        /*- The amount of caching of data to perform -*/
        options.add_int("CACHELEVEL", 2);
        /*- The algorithm to use for the $\left\langle VV||VV\right \rangle$ terms -*/
        options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
        /*- The type of gauge to use for properties -*/
        options.add_str("GAUGE", "LENGTH");
        /*- Do relax the one-particle density matrix? -*/
        options.add_bool("OPDM_RELAX", false);
        /*- Do require $\bar{H}$ and $R$ to be connected? !expert -*/
        options.add_bool("XI_CONNECT", false);
        /*- The number of electronic states to computed, per irreducible
        representation -*/
        options.add("ROOTS_PER_IRREP", new ArrayType());
        /*- Compute non-relaxed properties for all excited states. -*/
        options.add_bool("PROP_ALL", true);
        /*- The symmetry of states -*/
        options.add_int("PROP_SYM", 1);
        /*- Root number (within its irrep) for computing properties -*/
        options.add_int("PROP_ROOT", 1);
        /*- Do compute Xi? -*/
        options.add_bool("XI", false);
        /*- Do use zeta?  -*/
        options.add_bool("ZETA", false);
        /*- For internal use only! Compute the one-particle density matrix, but not the two-particle density matrix. !expert -*/
        options.add_bool("OPDM_ONLY", false);
        /*- Do write natural orbitals (molden) -*/
        options.add_bool("WRITE_NOS", false);
        /*- Reproducing energies from densities ? -*/
        options.add_int("DEBUG", 0);
    }
    if (name == "CCLAMBDA" || options.read_globals()) {
        /*- MODULEDESCRIPTION Solves for the Lagrange multipliers, which are needed whenever coupled cluster properties
            or gradients are requested. -*/
        /*- Wavefunction type !expert -*/
        options.add_str("WFN", "SCF");
        /*- Convergence criterion for wavefunction (change) in CC lambda-amplitude equations. -*/
        options.add_double("R_CONVERGENCE", 1e-7);
        /*- Do restart the coupled-cluster iterations from old $\lambda@@1$ and $\lambda@@2$
        amplitudes? -*/
        options.add_bool("RESTART", false);
        /*- Caching level for libdpd governing the storage of amplitudes,
        integrals, and intermediates in the CC procedure. A value of 0 retains
        no quantities in cache, while a level of 6 attempts to store all
        quantities in cache.  For particularly large calculations, a value of
        0 may help with certain types of memory problems.  The default is 2,
        which means that all four-index quantities with up to two virtual-orbital
        indices (e.g., $\left\langle ij | ab \right\rangle$ integrals) may be held in the cache. -*/
        options.add_int("CACHELEVEL", 2);
        /*- Do Sekino-Bartlett size-extensive model-III? -*/
        options.add_bool("SEKINO", false);
        /*- Do use DIIS extrapolation to accelerate convergence? -*/
        options.add_bool("DIIS", true);
        /*- The algorithm to use for the $\left\langle VV||VV \right\rangle$ terms -*/
        options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
        /*- Type of ABCD algorithm will be used -*/
        options.add_str("ABCD", "NEW");
        /*- Number of important CC amplitudes per excitation level to print.
        CC analog to |detci__num_dets_print|. -*/
        options.add_int("NUM_AMPS_PRINT", 10);
        /*- Type of job being performed !expert -*/
        options.add_str("JOBTYPE", "");
        /*- Do simulate the effects of local correlation techniques? -*/
        options.add_bool("LOCAL", false);
        /*- Desired treatment of "weak pairs" in the local-CCSD method. The value of ``NONE`` (unique available option)
        treats weak pairs in the same manner as strong pairs. -*/
        options.add_str("LOCAL_WEAKP", "NONE");
        /*- Value (always between one and zero) for the Broughton-Pulay completeness
        check used to contruct orbital domains for local-CC calculations. See
        J. Broughton and P. Pulay, J. Comp. Chem. 14, 736-740 (1993) and C. Hampel
        and H.-J. Werner, J. Chem. Phys. 104, 6286-6297 (1996). -*/
        options.add_double("LOCAL_CUTOFF", 0.02);
        /*- Type of local-CCSD scheme to be simulated. ``WERNER`` (unique available option) selects the method
        developed by H.-J. Werner and co-workers. -*/
        options.add_str("LOCAL_METHOD", "WERNER");
        /*- Do apply local filtering to single de-excitation ($\lambda 1$ amplitudes? -*/
        options.add_bool("LOCAL_FILTER_SINGLES", true);
        /*- Cutoff value for local-coupled-perturbed-Hartree-Fock -*/
        options.add_double("LOCAL_CPHF_CUTOFF", 0.10);
        /*- Definition of local pair domains -*/
        options.add_str("LOCAL_PAIRDEF", "");
        /*- The number of electronic states to computed, per irreducible
        representation -*/
        options.add("ROOTS_PER_IRREP", new ArrayType());
        /*- Compute unrelaxed properties for all excited states. -*/
        options.add_bool("PROP_ALL", true);
        /*- The symmetry of states -*/
        options.add_int("PROP_SYM", 1);
        /*- Root number (within its irrep) for computing properties -*/
        options.add_int("PROP_ROOT", 1);
        /*- Maximum number of iterations -*/
        options.add_int("MAXITER", 50);
        /*- Do use zeta?  -*/
        options.add_bool("ZETA", false);
    }
    if (name == "ADC" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs Algebraic-Diagrammatic Construction (ADC) propagator computations for excited
           states. -*/
        /*- Reference wavefunction type -*/
        options.add_str("REFERENCE", "RHF", "RHF UHF");
        /*- The number of poles / excited states to obtain per irrep vector -*/
        options.add("ROOTS_PER_IRREP", new ArrayType());
        /*- Tolerance for extracted or printed amplitudes. This option is only available for the adcc backend. -*/
        options.add_double("CUTOFF_AMPS_PRINT", 0.01);
        /*- Convergence threshold for ADC matrix diagonalisation. Negative values keep the
         *   adcc default (1e-6) -*/
        options.add_double("R_CONVERGENCE", -1);
        /*- Number of guess vectors to generate and use. Negative values keep
         *  the adcc default (currently 2 * ROOTS_PER_IRREP). This option is only available for the adcc backend. -*/
        options.add_int("NUM_GUESSES", -1);
        /*- Number of orbitals to place in the core. This option is only available for the adcc backend.  -*/
        options.add_int("NUM_CORE_ORBITALS", 0);
        /*- The kind of states to compute. -*/
        options.add_str("KIND", "SINGLET", "SINGLET TRIPLET SPIN_FLIP ANY");
        /*- Maximum number of iterations -*/
        options.add_int("MAXITER", 50);
        /*- Maximum number of subspace vectors. A negative value uses
         *  the adcc default (roughly between 20 and 5 * N_GUESSES). This option is only available for the adcc backend. -*/
        options.add_int("MAX_NUM_VECS", -1);
        /*- Specifies the choice of representation of the electric dipole operator.
         *  Acceptable values are ``LENGTH`` (default) and ``VELOCITY``. -*/
        options.add_str("GAUGE", "LENGTH", "LENGTH VELOCITY");
    }
    if (name == "CCHBAR" || options.read_globals()) {
        /*- MODULEDESCRIPTION Assembles the coupled cluster effective Hamiltonian. Called whenever CC
            properties and/or gradients are required. -*/
        /*- Wavefunction type !expert -*/
        options.add_str("WFN", "SCF");
        /*- Reference wavefunction type for EOM computations -*/
        options.add_str("EOM_REFERENCE", "RHF");
        /*- Do compute the T amplitude equation matrix elements? -*/
        options.add_bool("T_AMPS", false);
        /*- Caching level for libdpd governing the storage of amplitudes,
        integrals, and intermediates in the CC procedure. A value of 0 retains
        no quantities in cache, while a level of 6 attempts to store all
        quantities in cache.  For particularly large calculations, a value of
        0 may help with certain types of memory problems.  The default is 2,
        which means that all four-index quantities with up to two virtual-orbital
        indices (e.g., $\langle ij | ab \rangle$ integrals) may be held in the cache. -*/
        options.add_int("CACHELEVEL", 2);
        /*- Do use the minimal-disk algorithm for Wabei? It's VERY slow! -*/
        options.add_bool("WABEI_LOWDISK", false);
    }
    if (name == "CCEOM" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs equation-of-motion (EOM) coupled cluster excited state computations. -*/
        /*- Wavefunction type !expert -*/
        options.add_str("WFN", "EOM_CCSD", "EOM_CCSD EOM_CC2 EOM_CC3");
        /*- Reference wavefunction type -*/
        options.add_str("REFERENCE", "RHF", "RHF ROHF UHF");
        /*- Reference wavefunction type for EOM computations -*/
        options.add_str("EOM_REFERENCE", "RHF", "RHF ROHF UHF");
        /*- Do use full effective Hamiltonian matrix? -*/
        options.add_bool("FULL_MATRIX", false);
        /*- Caching level for libdpd governing the storage of amplitudes,
        integrals, and intermediates in the CC procedure. A value of 0 retains
        no quantities in cache, while a level of 6 attempts to store all
        quantities in cache.  For particularly large calculations, a value of
        0 may help with certain types of memory problems.  The default is 2,
        which means that all four-index quantities with up to two virtual-orbital
        indices (e.g., $\left\langle ij | ab \right\rangle$ integrals) may be held in the cache. -*/
        options.add_int("CACHELEVEL", 2);
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
        /*- When collapsing Davidson subspace, whether to also include the
        previous approximate solution (for each root)? This doubles the
        number of resulting vectors but generally improves convergence. -*/
        options.add_bool("COLLAPSE_WITH_LAST", true);
        /*- Has the same effect as "COLLAPSE_WITH_LAST" but only in
        CC3 computations and after the initial solution of EOM CCSD.
        May help efficiency, but hazardous when solving for higher roots. -*/
        options.add_bool("COLLAPSE_WITH_LAST_CC3", false);
        /*- Complex tolerance applied in CCEOM computations -*/
        options.add_double("COMPLEX_TOLERANCE", 1E-12);
        /*- Convergence criterion for norm of the residual vector in the Davidson algorithm for CC-EOM. -*/
        options.add_double("R_CONVERGENCE", 1E-6);
        /*- Convergence criterion for norm of the residual vector in the Davidson algorithm for the CIS guess to CC-EOM.
           -*/
        options.add_double("SS_R_CONVERGENCE", 1E-6);
        /*- Convergence criterion for excitation energy (change) in the
        Davidson algorithm for CC-EOM. See Table :ref:`Post-SCF Convergence
        <table:conv_corl>` for default convergence criteria for different
        calculation types. -*/
        options.add_double("E_CONVERGENCE", 1E-6);
        /*- Convergence criterion for excitation energy (change) in the Davidson algorithm for the CIS guess to CC-EOM.
           -*/
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
        /*- Report overlaps with old excited-state wave functions, if
           available, and store current wave functions for later use. -*/
        options.add_bool("OVERLAP_CHECK", false);
    }
    if (name == "CCRESPONSE" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs coupled cluster response property computations. -*/
        /*- Wavefunction type !expert -*/
        options.add_str("WFN", "SCF");
        /*- Reference wavefunction type -*/
        options.add_str("REFERENCE", "RHF");
        /*- Caching level for libdpd -*/
        options.add_int("CACHELEVEL", 2);
        /*- Specifies the choice of representation of the electric dipole operator.
        For polarizability, this keyword is ignored and ``LENGTH`` gauge is computed.
        For optical rotation and raman optical activity, this keyword is active, and
        acceptable values are ``LENGTH`` for the usual length-gauge representation,
        ``VELOCITY``(default) for the modified velocity-gauge representation in which the
        static-limit optical rotation tensor is subtracted from the frequency-
        dependent tensor, or ``BOTH``. Note that, for optical rotation and raman optical
        activity calculations, only the choices of ``VELOCITY`` or ``BOTH`` will yield
        origin-independent results. -*/
        options.add_str("GAUGE", "VELOCITY", "LENGTH VELOCITY BOTH");
        /*- Maximum number of iterations to converge perturbed amplitude equations -*/
        options.add_int("MAXITER", 50);
        /*- Convergence criterion for wavefunction (change) in perturbed CC equations. -*/
        options.add_double("R_CONVERGENCE", 1e-7);
        /*- Do use DIIS extrapolation to accelerate convergence? -*/
        options.add_bool("DIIS", 1);
        /*- The response property desired.  Acceptable values are ``POLARIZABILITY``
        (default) for dipole polarizabilities, ``ROTATION`` for specific rotations,
        ``ROA`` for Raman Optical Activity (``ROA_TENSOR`` for each displacement),
        and ``ALL`` for all of the above. -*/
        options.add_str("PROPERTY", "POLARIZABILITY", "POLARIZABILITY ROTATION ROA ROA_TENSOR ALL");
        /*- Type of ABCD algorithm will be used -*/
        options.add_str("ABCD", "NEW");
        /*- Do restart from on-disk amplitudes? -*/
        options.add_bool("RESTART", 1);
        /*- Do simulate local correlation? -*/
        options.add_bool("LOCAL", 0);
        /*- Value (always between one and zero) for the Broughton-Pulay completeness
        check used to contruct orbital domains for local-CC calculations. See
        J. Broughton and P. Pulay, J. Comp. Chem. 14, 736-740 (1993) and C. Hampel
        and H.-J. Werner, J. Chem. Phys. 104, 6286-6297 (1996). -*/
        options.add_double("LOCAL_CUTOFF", 0.01);
        /*- Type of local-CCSD scheme to be simulated. ``WERNER`` (unique available option) selects the method
        developed by H.-J. Werner and co-workers. -*/
        options.add_str("LOCAL_METHOD", "WERNER");
        /*- Desired treatment of "weak pairs" in the local-CCSD method. The value of ``NONE`` (unique available option)
        treats weak pairs in the same manner as strong pairs. -*/
        options.add_str("LOCAL_WEAKP", "NONE");
        /*- Do apply local filtering to single excitation amplitudes? -*/
        options.add_bool("LOCAL_FILTER_SINGLES", false);
        /*- Cutoff value for local-coupled-perturbed-Hartree-Fock -*/
        options.add_double("LOCAL_CPHF_CUTOFF", 0.10);
        /*- Definition of local pair domains -*/
        options.add_str("LOCAL_PAIRDEF", "NONE");
        /*- Do analyze X2 amplitudes -*/
        options.add_bool("ANALYZE", 0);
        /*- Number of important CC amplitudes per excitation level to print.
        CC analog to |detci__num_dets_print|. -*/
        options.add_int("NUM_AMPS_PRINT", 5);
        /*- Do Sekino-Bartlett size-extensive model-III? -*/
        options.add_bool("SEKINO", 0);
        /*- Do Bartlett size-extensive linear model? -*/
        options.add_bool("LINEAR", 0);
        /*- Array that specifies the desired frequencies of the incident
        radiation field in CCLR calculations.  If only one element is
        given, the units will be assumed to be atomic units.  If more
        than one element is given, then the units must be specified as the final
        element of the array.  Acceptable units are ``HZ``, ``NM``, ``EV``, and ``AU``. -*/
        options.add("OMEGA", new ArrayType());
    }
    //  if(name == "RESPONSE"|| options.read_globals()){
    //     /*- MODULEDESCRIPTION Performs SCF linear response computations. -*/
    //    /*- Reference wavefunction type -*/
    //    options.add_str("REFERENCE", "RHF");
    //    /*- Array that specifies the desired frequencies of the incident
    //    radiation field in CCLR calculations.  If only one element is
    //    given, the units will be assumed to be atomic units.  If more
    //    than one element is given, then the units must be specified as the final
    //    element of the array.  Acceptable units are ``HZ``, ``NM``, ``EV``, and ``AU``. -*/
    //    options.add("OMEGA", new ArrayType());
    //    /*- Array that specifies the desired frequencies of the incident
    //    radiation field in CCLR calculations.  If only one element is
    //    given, the units will be assumed to be atomic units.  If more
    //    than one element is given, then the units must be specified as the final
    //    element of the array.  Acceptable units are HZ, NM, EV, and AU. -*/
    //    /*- The response property desired.  Acceptable values are POLARIZABILITY
    //    (default) for dipole-polarizabilities, ROTATION for specific rotations,
    //    ROA for Raman Optical Activity, and ALL for all of the above.
    //    -*/
    //    options.add_str("PROPERTY","POLARIZABILITY","POLARIZABILITY ROTATION ROA ALL");
    //  }
    if (name == "MCSCF" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs RHF/UHF/ROHF/TCSCF and more general MCSCF computations. Called
            as the starting point for multireference coupled cluster computations. -*/
        /*- Reference wavefunction type -*/
        options.add_str("REFERENCE", "RHF", "RHF ROHF UHF TWOCON MCSCF GENERAL");
        /*- Level shift to aid convergence -*/
        options.add_double("LEVEL_SHIFT", 0.0);
        /*- Convergence criterion for energy. -*/
        options.add_double("E_CONVERGENCE", 1e-6);
        /*- Convergence criterion for density, as measured by the orbital gradient. -*/
        options.add_double("D_CONVERGENCE", 1e-6);
        /*- Maximum number of iterations -*/
        options.add_int("MAXITER", 100);
        /*- Maximum number of error vectors stored for DIIS extrapolation -*/
        options.add_int("DIIS_MAX_VECS", 7);
        /*- Which solution of the SCF equations to find, where 1 is the SCF ground state-*/
        options.add_int("FOLLOW_ROOT", 1);
        /*- Iteration at which to begin using the averaged Fock matrix-*/
        options.add_int("FAVG_START", 5);
        /*- -*/
        options.add_int("TURN_ON_ACTV", 0);
        /*- For orbital rotations after convergence, the angle (in degrees) by which to rotate. !expert -*/
        options.add_double("ROTATE_MO_ANGLE", 0.0);
        /*- For orbital rotations after convergence, irrep (1-based, Cotton order) of the orbitals to rotate. !expert
           -*/
        options.add_int("ROTATE_MO_IRREP", 1);
        /*- For orbital rotations after convergence, number of the first orbital (1-based) to rotate. !expert -*/
        options.add_int("ROTATE_MO_P", 1);
        /*- For orbital rotations after convergence, number of the second orbital (1-based) to rotate. !expert -*/
        options.add_int("ROTATE_MO_Q", 2);
        /*- Do use DIIS extrapolation to accelerate convergence of the CI coefficients? -*/
        options.add_bool("CI_DIIS", false);
        /*- Do use DIIS extrapolation to accelerate convergence of the SCF energy (MO coefficients only)? -*/
        options.add_bool("DIIS", true);
        /*- Do read in from file the MOs from a previous computation? -*/
        options.add_bool("MO_READ", true);
        /*- Do use the average Fock matrix during the SCF optimization? -*/
        options.add_bool("FAVG", false);
        /*- Do canonicalize the active orbitals such that the average Fock matrix is diagonal? -*/
        options.add_bool("CANONICALIZE_ACTIVE_FAVG", false);
        /*- Do canonicalize the inactive (DOCC and Virtual) orbitals such that the average Fock matrix is diagonal? -*/
        options.add_bool("CANONICALIZE_INACTIVE_FAVG", false);
        /*- Do consider internal rotations? -*/
        options.add_bool("INTERNAL_ROTATIONS", true);
        /*- Do attempt to force a two configuration solution by starting with CI coefficents of $\pm \sqrt{\frac{1}{2}}$
           ? -*/
        options.add_bool("FORCE_TWOCON", false);
        /*- The number of singly occupied orbitals, per irrep -*/
        options.add("SOCC", new ArrayType());
        /*- The number of doubly occupied orbitals, per irrep -*/
        options.add("DOCC", new ArrayType());
        /*- The symmetry of the SCF wavefunction.-*/
        options.add_str("WFN_SYM", "1",
                        "A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");
    }
    if (name == "CCENERGY" || options.read_globals()) {
        /*- MODULEDESCRIPTION Computes coupled cluster energies. Called as part of any coupled cluster computation. -*/

        /*- Wavefunction type !expert -*/
        options.add_str("WFN", "NONE",
                        "CCSD CCSD_T CCSD_AT EOM_CCSD BCCD BCCD_T CC2 CC3 EOM_CC2 EOM_CC3 CCSD_MVD");
        /*- Reference wavefunction type -*/
        options.add_str("REFERENCE", "RHF", "RHF ROHF UHF");
        /*- Do use new triples? -*/
        options.add_bool("NEW_TRIPLES", 1);
        /*- Do analyze T2 amplitudes -*/
        options.add_bool("ANALYZE", 0);
        /*- Maximum number of iterations to solve the CC equations -*/
        options.add_int("MAXITER", 50);
        /*- Convergence criterion for energy. See Table :ref:`Post-SCF
        Convergence <table:conv_corl>` for default convergence criteria for
        different calculation types. -*/
        options.add_double("E_CONVERGENCE", 1e-6);
        /*- Convergence criterion for wavefunction (change) in CC amplitude equations. -*/
        options.add_double("R_CONVERGENCE", 1e-7);
        /*- Do restart the coupled-cluster iterations from old $t@@1$ and $t@@2$
        amplitudes?  For geometry optimizations, Brueckner
        calculations, etc. the iterative solution of the CC amplitude
        equations may benefit considerably by reusing old vectors as initial
        guesses.  Assuming that the MO phases remain the same between
        updates, the CC codes will, by default, re-use old vectors, unless
        the user sets RESTART = false. -*/
        options.add_bool("RESTART", 1);
        /*- Do restart the coupled-cluster iterations even if MO phases are screwed up? !expert -*/
        options.add_bool("FORCE_RESTART", 0);
        //#warning CCEnergy ao_basis keyword type was changed.
        /*- The algorithm to use for the $\left\langle VV||VV\right\rangle$ terms
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
        /*- Caching level for libdpd governing the storage of amplitudes,
        integrals, and intermediates in the CC procedure. A value of 0 retains
        no quantities in cache, while a level of 6 attempts to store all
        quantities in cache.  For particularly large calculations, a value of
        0 may help with certain types of memory problems.  The default is 2,
        which means that all four-index quantities with up to two virtual-orbital
        indices (e.g., $\langle ij | ab \rangle$ integrals) may be held in the cache. -*/
        options.add_int("CACHELEVEL", 2);
        /*- Selects the priority type for maintaining the automatic memory
        cache used by the libdpd codes. A value of ``LOW`` selects a "low priority"
        scheme in which the deletion of items from the cache is based on
        pre-programmed priorities. A value of LRU selects a "least recently used"
        scheme in which the oldest item in the cache will be the first one deleted. -*/
        options.add_str("CACHETYPE", "LOW", "LOW LRU");
        /*- Number of threads -*/
        options.add_int("CC_NUM_THREADS", 1);
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
        // options.add_int("LOCAL_FILTER_SINGLES", 1);
        /*- Cutoff value for local-coupled-perturbed-Hartree-Fock -*/
        options.add_double("LOCAL_CPHF_CUTOFF", 0.10);
        /*- Definition of local pair domains, default is BP, Boughton-Pulay. -*/
        options.add_str("LOCAL_PAIRDEF", "BP", "BP RESPONSE");
        /*- Number of important $t@@1$ and $t@@2$ amplitudes to print -*/
        options.add_int("NUM_AMPS_PRINT", 10);
        /*- Convergence criterion for Brueckner orbitals. The convergence is
        determined based on the largest $T_1$ amplitude.  Default adjusts
        depending on |ccenergy__e_convergence|. -*/
        options.add_double("BRUECKNER_ORBS_R_CONVERGENCE", 1e-5);
        /*- Do print the MP2 amplitudes which are the starting guesses for RHF and UHF reference functions? -*/
        options.add_bool("MP2_AMPS_PRINT", 0);
        /*- Do print MP2 and CCSD pair energies for RHF references? -*/
        options.add_bool("PAIR_ENERGIES_PRINT", 0);
        /*- Do build W intermediates required for cc3 in core memory? -*/
        options.add_bool("T3_WS_INCORE", 0);
        /*- Do SCS-MP2 with parameters optimized for nucleic acids? -*/
        options.add_bool("SCSN_MP2", 0);
        /*- Do spin-component-scaled MP2 (SCS-MP2)? -*/
        options.add_bool("SCS_MP2", 0);
        /*- Do spin-component-scaled CCSD -*/
        options.add_bool("SCS_CCSD", 0);
        /*- MP2 opposite-spin scaling value -*/
        options.add_double("MP2_OS_SCALE", 1.20);
        /*- MP2 same-spin scaling value -*/
        options.add_double("MP2_SS_SCALE", 1.0 / 3.0);
        /*- Coupled-cluster opposite-spin scaling value -*/
        options.add_double("CC_OS_SCALE", 1.27);
        /*- Coupled-cluster same-spin scaling value -*/
        options.add_double("CC_SS_SCALE", 1.13);
        /*- Convert ROHF MOs to semicanonical MOs -*/
        options.add_bool("SEMICANONICAL", true);
        /*- Maximum number of iterations for Brueckner CCD. -*/
        options.add_int("BCCD_MAXITER", 50);
    }
    if (name == "DFMP2" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs density-fitted MP2 computations for RHF/UHF/ROHF reference wavefunctions. -*/

        /*- A helpful option, used only in debugging the MADNESS version !expert-*/
        options.add_int("MADMP2_SLEEP", 0);
        /*- Primary basis set -*/
        options.add_str("BASIS", "NONE");
        /*- Auxiliary basis set for MP2 density fitting computations.
        :ref:`Defaults <apdx:basisFamily>` to a RI basis. -*/
        options.add_str("DF_BASIS_MP2", "");
        /*- OS Scale -*/
        options.add_double("MP2_OS_SCALE", 6.0 / 5.0);
        /*- SS Scale  -*/
        options.add_double("MP2_SS_SCALE", 1.0 / 3.0);
        /*- \% of memory for DF-MP2 three-index buffers -*/
        options.add_double("DFMP2_MEM_FACTOR", 0.9);
        /*- Schwarz screening threshold. Mininum absolute value below which TEI are neglected. -*/
        options.add_double("INTS_TOLERANCE", 0.0);
        /*- Minimum error in the 2-norm of the P(2) matrix for corrections to Lia and P. -*/
        options.add_double("DFMP2_P2_TOLERANCE", 0.0);
        /*- Minimum error in the 2-norm of the P matrix for skeleton-core Fock matrix derivatives. -*/
        options.add_double("DFMP2_P_TOLERANCE", 0.0);
        /*- Number of threads to compute integrals with. 0 is wild card -*/
        options.add_int("DF_INTS_NUM_THREADS", 0);
        /*- IO caching for CP corrections, etc !expert -*/
        options.add_str("DF_INTS_IO", "NONE", "NONE SAVE LOAD");
        /*- Do relax the one-particle density matrix? -*/
        options.add_bool("OPDM_RELAX", true);
        /*- Do compute one-particle density matrix? -*/
        options.add_bool("ONEPDM", false);
    }
    if (name == "DFEP2" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs density-fitted EP2 computations for RHF reference wavefunctions. -*/

        /*- Auxiliary basis set for EP2 density fitting computations.
        :ref:`Defaults <apdx:basisFamily>` to a RI basis. -*/
        options.add_str("DF_BASIS_EP2", "");
        /*- Number of Ionization Potentials to compute, starting with the HOMO. -*/
        options.add_int("EP2_NUM_IP", 3);
        /*- Number of Electron Affinities to compute, starting with the LUMO. -*/
        options.add_int("EP2_NUM_EA", 0);
        /*- Explicitly pick orbitals to use in the EP2 method, overrides |dfep2__EP2_NUM_IP| and |dfep2__ep2_num_ea|
           options. Input array should be [[orb1, orb2], [], ...] for each irrep. -*/
        options.add("EP2_ORBITALS", new ArrayType());
        /*- What is the maximum number of iterations? -*/
        options.add_double("EP2_CONVERGENCE", 5.e-5);
        /*- What is the maximum number of iterations? -*/
        options.add_int("EP2_MAXITER", 20);
    }
    if (name == "DLPNO" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs DLPNO-MP2 computations for RHF reference wavefunctions. -*/

        /*- SUBSECTION General Options -*/

        /*- Auxiliary basis set for MP2 density fitting computations.
        :ref:`Defaults <apdx:basisFamily>` to a RI basis. -*/
        options.add_str("DF_BASIS_MP2", "");
        /*- General convergence criteria for DLPNO methods -*/
        options.add_str("PNO_CONVERGENCE", "NORMAL", "LOOSE NORMAL TIGHT");
        /*- Convergence criteria for the Foster-Boys orbital localization -*/
        options.add_double("LOCAL_CONVERGENCE", 1.0E-12);
        /*- Maximum iterations in Foster-Boys localization -*/
        options.add_int("LOCAL_MAXITER", 1000);
        /*- Energy convergence criteria for local MP2 iterations -*/
        options.add_double("E_CONVERGENCE", 1e-6);
        /*- Residual convergence criteria for local MP2 iterations -*/
        options.add_double("R_CONVERGENCE", 1e-6);
        /*- Orbital localizer -*/
        options.add_str("DLPNO_LOCAL_ORBITALS", "BOYS", "BOYS PIPEK_MEZEY");
        /*- Maximum number of iterations to determine the MP2 amplitudes. -*/
        options.add_int("DLPNO_MAXITER", 50);

        /*- SUBSECTION Expert Options -*/

        /*- Which DLPNO Algorithm to run (not set by user) !expert -*/
        options.add_str("DLPNO_ALGORITHM", "MP2", "MP2");
        /*- Occupation number threshold for removing PNOs !expert -*/
        options.add_double("T_CUT_PNO", 1e-8);
        /*- DOI threshold for including PAO (u) in domain of LMO (i) !expert -*/
        options.add_double("T_CUT_DO", 1e-2);
        /*- DOI threshold for treating LMOs (i,j) as interacting !expert -*/
        options.add_double("T_CUT_DO_ij", 1e-5);
        /*- Pair energy threshold (dipole approximation) for treating LMOs (i, j) as interacting !expert -*/
        options.add_double("T_CUT_PRE", 1e-6); 
        /*- DOI threshold for including PAO (u) in domain of LMO (i) during pre-screening !expert -*/
        options.add_double("T_CUT_DO_PRE", 3e-2);
        /*- Mulliken charge threshold for including aux BFs on atom (a) in domain of LMO (i) !expert -*/
        options.add_double("T_CUT_MKN", 1e-3);
        /*- Basis set coefficient threshold for including basis function (m) in domain of LMO (i) !expert -*/
        options.add_double("T_CUT_CLMO", 1e-2);
        /*- Basis set coefficient threshold for including basis function (n) in domain of PAO (u) !expert -*/
        options.add_double("T_CUT_CPAO", 1e-3);
        /*- Overlap matrix threshold for removing linear dependencies !expert -*/
        options.add_double("S_CUT", 1e-8);
        /*- Fock matrix threshold for treating ampltudes as coupled during local MP2 iterations !expert -*/
        options.add_double("F_CUT", 1e-5);

        /*- SUBSECTION DOI Grid Options -*/

        /*- Number of spherical points in DOI grid !expert -*/
        options.add_int("DOI_SPHERICAL_POINTS", 50);
        /*- Number of radial points in DOI grid !expert -*/
        options.add_int("DOI_RADIAL_POINTS", 25);
        /*- Screening criteria for basis function values on DOI grids !expert -*/
        options.add_double("DOI_BASIS_TOLERANCE", 1.0E-10);
        /*- Pruning scheme for DOI grids !expert -*/
        options.add_str("DOI_PRUNING_SCHEME", "ROBUST", "ROBUST TREUTLER NONE FLAT P_GAUSSIAN D_GAUSSIAN P_SLATER D_SLATER LOG_GAUSSIAN LOG_SLATER NONE");
    }
    if (name == "PSIMRCC" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs multireference coupled cluster computations.  This theory
           should be used only by advanced users with a good working knowledge of multireference
           techniques. -*/

        /*- The multiplicity, $M@@S(M@@S+1)$, of the target state.  Must be specified if different from the reference
           $M@@s$. -*/
        options.add_int("CORR_MULTP", 1);
        /*- The molecular charge of the target state -*/
        options.add_int("CORR_CHARGE", 0);
        /*- The amount (percentage) of damping to apply to the amplitude updates.
            0 will result in a full update, 100 will completely stall the update. A
            value around 20 (which corresponds to 20\% of the amplitudes from the
            previous iteration being mixed into the current iteration)
            can help in cases where oscillatory convergence is observed. -*/
        options.add_double("DAMPING_PERCENTAGE", 0.0);
        /*- Maximum number of error vectors stored for DIIS extrapolation -*/
        options.add_int("DIIS_MAX_VECS", 7);
        /*- Number of threads -*/
        options.add_int("CC_NUM_THREADS", 1);
        /*- Which root of the effective hamiltonian is the target state? -*/
        options.add_int("FOLLOW_ROOT", 1);
        /*- Convergence criterion for energy. See Table :ref:`Post-SCF
        Convergence <table:conv_corl>` for default convergence criteria for
        different calculation types. -*/
        options.add_double("E_CONVERGENCE", 1e-6);
        /*- Convergence criterion for amplitudes (residuals). -*/
        options.add_double("R_CONVERGENCE", 1e-9);
        /*- Maximum number of iterations to determine the amplitudes -*/
        options.add_int("MAXITER", 100);
        /*- The number of DIIS vectors needed before extrapolation is performed -*/
        options.add_int("DIIS_START", 2);
        /*- The shift to apply to the denominators, {\it c.f.} Taube and Bartlett, JCP, 130, 144112 (2009) -*/
        options.add_double("TIKHONOW_OMEGA", 0.0);  // Omega = TIKHONOW_OMEGA / 1000
        /*- The cycle after which Tikhonow regularization is stopped.
        Set to zero to allow regularization in all iterations -*/
        options.add_int("TIKHONOW_MAX", 5);
        /*- Do use DIIS extrapolation to accelerate convergence for iterative triples excitations? -*/
        options.add_bool("TRIPLES_DIIS", false);
        /*- Do lock onto a singlet root? -*/
        options.add_bool("LOCK_SINGLET", false);
        /*- Do start from a MP2 guess? -*/
        options.add_bool("MP2_GUESS", true);
        /*- Do use the averaged Fock matrix over all references in (T) computations? -*/
        options.add_bool("FAVG_CCSD_T", false);
        /*- Do include the fourth-order contributions to the effective Hamiltonian? -*/
        options.add_bool("HEFF4", true);
        /*- Do include the off-diagonal corrections in (T) computations? -*/
        options.add_bool("OFFDIAGONAL_CCSD_T", true);
        /*- Do include the diagonal corrections in (T) computations? -*/
        options.add_bool("DIAGONAL_CCSD_T", true);
        /*- Do diagonalize the effective Hamiltonian? -*/
        options.add_bool("DIAGONALIZE_HEFF", false);
        /*- Do use symmetry to map equivalent determinants onto each other, for efficiency? -*/
        options.add_bool("USE_SPIN_SYM", true);
        /*- Do zero the internal amplitudes, i.e., those that map reference determinants onto each other? -*/
        options.add_bool("ZERO_INTERNAL_AMPS", true);
        /*- Do include the terms that couple the reference determinants? -*/
        options.add_bool("COUPLING_TERMS", true);
        /*- Do print the effective Hamiltonian? -*/
        options.add_bool("HEFF_PRINT", false);
        /*- Removed in 1.4. Will raise an error in 1.5. -*/
        options.add_bool("PERTURB_CBS", false);
        /*- Removed in 1.4. Will raise an error in 1.5. -*/
        options.add_bool("PERTURB_CBS_COUPLING", true);
        /*- Do use Tikhonow regularization in (T) computations? !expert -*/
        options.add_bool("TIKHONOW_TRIPLES", false);
        /*- The type of perturbation theory computation to perform -*/
        options.add_str("PT_ENERGY", "SECOND_ORDER",
                        "SECOND_ORDER SCS_SECOND_ORDER PSEUDO_SECOND_ORDER SCS_PSEUDO_SECOND_ORDER");
        /*- The type of correlated wavefunction -*/
        options.add_str("CORR_WFN", "CCSD", "PT2 CCSD MP2-CCSD CCSD_T");
        /*- The type of CCSD(T) computation to perform -*/
        options.add_str("CORR_CCSD_T", "STANDARD", "STANDARD PITTNER");
        /*- The ansatz to use for MRCC computations -*/
        options.add_str("CORR_ANSATZ", "MK", "SR MK BW APBW");
        /*- The order of coupling terms to include in MRCCSDT computations -*/
        options.add_str("COUPLING", "CUBIC", "NONE LINEAR QUADRATIC CUBIC");
        /*- The symmetry of the target wavefunction, specified either by Sch\ |o_dots|\ nflies symbol,
            or irrep number (in Cotton ordering) -*/
        options.add_str("WFN_SYM", "1",
                        "A AG AU AP APP A1 A2 B BG BU B1 B2 B3 B1G B2G B3G B1U B2U B3U 0 1 2 3 4 5 6 7 8");
        /*- The type of algorithm to use for (T) computations -*/
        options.add_str("TRIPLES_ALGORITHM", "RESTRICTED", "SPIN_ADAPTED RESTRICTED UNRESTRICTED");
        /*- How to perform MP2_CCSD computations -*/
        options.add_str("MP2_CCSD_METHOD", "II", "I IA II");
        /*- Whether to use spin symmetry to map equivalent configurations onto each other, for efficiency !expert -*/
        options.add_bool("USE_SPIN_SYMMETRY", true);
        // /*- The number of frozen occupied orbitals per irrep -*/
        // options.add("FROZEN_DOCC", new ArrayType());
        // /*- The number of doubly occupied orbitals per irrep -*/
        // options.add("RESTRICTED_DOCC", new ArrayType());
        // /*- The number of active orbitals per irrep -*/
        // options.add("ACTIVE", new ArrayType());
        // /*- The number of frozen virtual orbitals per irrep -*/
        // options.add("FROZEN_UOCC", new ArrayType());
        /*- -*/
        options.add_int("SMALL_CUTOFF", 0);
        /*- Do disregard updating single excitation amplitudes? -*/
        options.add_bool("NO_SINGLES", false);
    }
    if (name == "OPTKING" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs geometry optimizations and vibrational frequency analyses. -*/

        /*- SUBSECTION Optimization Algorithm -*/

        /*- Maximum number of geometry optimization steps -*/
        options.add_int("GEOM_MAXITER", 50);
        /*- Print all optking parameters. -*/
        options.add_bool("PRINT_OPT_PARAMS", false);
        /*- Specifies minimum search, transition-state search, or IRC following -*/
        options.add_str("OPT_TYPE", "MIN", "MIN TS IRC");
        /*- Geometry optimization step type, either Newton-Raphson or Rational Function Optimization -*/
        options.add_str("STEP_TYPE", "RFO", "RFO RS_I_RFO P_RFO NR SD LINESEARCH");
        /*- Geometry optimization coordinates to use.
            REDUNDANT and INTERNAL are synonyms and the default.
            CARTESIAN uses only cartesian coordinates.
            BOTH uses both redundant and cartesian coordinates.
            CUSTOM is not fully implemented yet - expected optking 0.3.1  -*/
        options.add_str("OPT_COORDINATES", "INTERNAL", "REDUNDANT INTERNAL CARTESIAN BOTH CUSTOM");
        /*- A string formatted as a dicitonary containing a set of coordinates. Coordinates can be
            appended to Optking's coordinate set or used on their own - expected optking 0.3.1. -*/
        options.add_str("CUSTOM_COORDS", "");
        /*- Do follow the initial RFO vector after the first step? -*/
        options.add_bool("RFO_FOLLOW_ROOT", false);
        /*- Root for RFO to follow, 0 being lowest (for a minimum) -*/
        options.add_int("RFO_ROOT", 0);
        /*- Starting level for dynamic optimization (0=nondynamic, higher=>more conservative) -*/
        options.add_int("DYNAMIC_LEVEL", 0);
        /*- IRC step size in bohr(amu)\ $^{1/2}$. -*/
        options.add_double("IRC_STEP_SIZE", 0.2);
        /*- IRC mapping direction -*/
        options.add_str("IRC_DIRECTION", "FORWARD", "FORWARD BACKWARD");
        /*- Maximum number of IRC points to collect before stopping. -*/
        options.add_int("IRC_POINTS", 20);
        /*- Initial maximum step size in bohr or radian along an internal coordinate -*/
        options.add_double("INTRAFRAG_STEP_LIMIT", 0.5);
        /*- Lower bound for dynamic trust radius [au] -*/
        options.add_double("INTRAFRAG_STEP_LIMIT_MIN", 0.001);
        /*- Upper bound for dynamic trust radius [au] -*/
        options.add_double("INTRAFRAG_STEP_LIMIT_MAX", 1.0);
        /*- Maximum step size in bohr or radian along an interfragment coordinate -*/
        options.add_double("INTERFRAG_STEP_LIMIT", 0.5);
        /*- Reduce step size as necessary to ensure back-transformation of internal
            coordinate step to cartesian coordinates. -*/
        options.add_bool("ENSURE_BT_CONVERGENCE", false);
        /*= Do stupid, linear scaling of internal coordinates to step limit (not RS-RFO) -*/
        options.add_bool("SIMPLE_STEP_SCALING", false);
        /*- Set number of consecutive backward steps allowed in optimization -*/
        options.add_int("CONSECUTIVE_BACKSTEPS", 0);
        /*- Eigenvectors of RFO matrix whose final column is smaller than this are ignored. -*/
        options.add_double("RFO_NORMALIZATION_MAX", 100);
        /*- Denominator check for hessian update. -*/
        options.add_double("H_UPDATE_DEN_TOL", 1e-7);
        /*- Absolute maximum value of RS-RFO. -*/
        options.add_double("RSRFO_ALPHA_MAX", 1e8);
        /*- Specify distances between atoms to be frozen (unchanged) -*/
        options.add_str("FROZEN_DISTANCE", "");
        /*- Specify angles between atoms to be frozen (unchanged) -*/
        options.add_str("FROZEN_BEND", "");
        /*- Specify dihedral angles between atoms to be frozen (unchanged) -*/
        options.add_str("FROZEN_DIHEDRAL", "");
        /*- Specify out-of-plane angles between atoms to be frozen (unchanged) -*/
        options.add_str("FROZEN_OOFP", "");
        /*- Specify atom and X, XY, XYZ, ... to be frozen (unchanged) -*/
        options.add_str("FROZEN_CARTESIAN", "");
        /*- Specify range for distances between atoms to be constrained to (eq. value specified)
            analogous to the previous FIXED_DISTANCE -*/
        options.add_str("RANGED_DISTANCE", "");
        /*- Specify range for angles between atoms to be constrained to (eq. value specified)
            analogous to the previous FIXED_BEND -*/
        options.add_str("RANGED_BEND", "");
        /*- Specify range for the dihedral angles between atoms to be constrained to (eq. value specified)
            analogous to the previous FIXED_DIHEDRAL -*/
        options.add_str("RANGED_DIHEDRAL", "");
        /*- Specify range for the out-of-plane angles between atoms to be constrained to (eq. value specified)
            analogous to the old FIXED_<COORD> keyword-*/
        options.add_str("RANGED_OOFP", "");
        /*- Freeze ALL dihedral angles -*/
        options.add_bool("FREEZE_ALL_DIHEDRALS", false);
        /*- Unfreeze a subset of dihedrals - meant for use with freeze_all_dihedrals -*/
        options.add_str("UNFREEZE_DIHEDRALS", "");

        /*- Specify formula for external forces for the distance between atoms -*/
        options.add_str("EXT_FORCE_DISTANCE", "");
        /*- Specify formula for external forces for angles between atoms -*/
        options.add_str("EXT_FORCE_BEND", "");
        /*- Specify formula for external forces for dihedral angles between atoms -*/
        options.add_str("EXT_FORCE_DIHEDRAL", "");
        /*- Specify formula for external forces for out-of-plane angles between atoms -*/
        options.add_str("EXT_FORCE_OOFP", "");
        /*- Symmetry formula for external forces for cartesian coordinates on atoms . -*/
        options.add_str("EXT_FORCE_CARTESIAN", "");
        /*- Tolerance for symmetrizing cartesian geometry between steps -*/
        options.add_double("CARTESIAN_SYM_TOLERANCE", 1e-7);

        /*- SUBSECTION Convergence Control -*/

        /*- Set of optimization criteria. Specification of any MAX_*_G_CONVERGENCE
        or RMS_*_G_CONVERGENCE options will append to overwrite the criteria set here
        unless |optking__flexible_g_convergence| is also on.      See Table :ref:`Geometry Convergence
        <table:optkingconv>` for details. -*/
        options.add_str("G_CONVERGENCE", "QCHEM", "QCHEM MOLPRO GAU GAU_LOOSE GAU_TIGHT INTERFRAG_TIGHT GAU_VERYTIGHT TURBOMOLE CFOUR NWCHEM_LOOSE");
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
        options.add_int("HESS_UPDATE_USE_LAST", 4);
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
        /* - Specify file for Cartesian Hessian (or JSON Schema-like file) to get Hessians from
        Only for experts - use |optking__full_hess_every| instead. Requires |optking__CART_HESS_READ | - */
        options.add_str("HESSIAN_FILE", "");
        /*- Frequency with which to compute the full Hessian in the course
        of a geometry optimization. 0 means to compute the initial Hessian only, 1
        means recompute every step, and N means recompute every N steps. The
        default (-1) is to never compute the full Hessian. -*/
        options.add_int("FULL_HESS_EVERY", -1);
        /*- Model Hessian to guess intrafragment force constants -*/
        options.add_str("INTRAFRAG_HESS", "SCHLEGEL", "FISCHER SCHLEGEL SIMPLE LINDH LINDH_SIMPLE");

        /*- SUBSECTION Fragment/Internal Coordinate Control -*/

        /*- For multi-fragment molecules, treat as single bonded molecule
        or via interfragment coordinates. A primary difference is that in ``MULTI`` mode,
        the interfragment coordinates are not redundant. -*/
        options.add_str("FRAG_MODE", "SINGLE", "SINGLE MULTI");
        /*- Which atoms define the reference points for interfragment coordinates? -*/
        options.add("FRAG_REF_ATOMS", new ArrayType());
        /*- Do freeze all fragments rigid? -*/
        options.add_bool("FREEZE_INTRAFRAG", false);
        /*- Do freeze all interfragment modes? -*/
        options.add_bool("FREEZE_INTERFRAG", false);
        /*- When interfragment coordinates are present, use as reference points either
        principal axes or fixed linear combinations of atoms. -*/
        options.add_str("INTERFRAG_MODE", "FIXED", "FIXED PRINCIPAL_AXES");
        /*- Use 1/R for the interfragment stretching coordinate instead of R -*/
        options.add_bool("INTERFRAG_DIST_INV", false);
        /*- Dictionary to define a dimer. Contains "Natoms per frag", "A Frag", "A Ref Atoms", "B Frag", and "B Ref Atoms" -*/
        options.add_str("INTERFRAG_COORDS", "");
        /*- Specify atoms to use for reference points in interfragment coordinates -*/
        options.add("FRAG_REF_ATOMS", new ArrayType());

        /*- Do add bond coordinates at nearby atoms for non-bonded systems? -*/
        options.add_bool("ADD_AUXILIARY_BONDS", true);
        /*- Re-estimate the Hessian at every step, i.e., ignore the currently stored Hessian. -*/
        options.add_bool("H_GUESS_EVERY", false);
        /*- This factor times standard covalent distance is used to add extra stretch coordinates. -*/
        options.add_double("AUXILIARY_BOND_FACTOR", 2.5);
        /*- Model Hessian to guess interfragment force constants -*/
        options.add_str("INTERFRAG_HESS", "DEFAULT", "DEFAULT FISCHER_LIKE");
        /*- When determining connectivity, a bond is assigned if interatomic distance
        is less than (this number) * sum of covalent radii. -*/
        options.add_double("COVALENT_CONNECT", 1.3);
        /*- When connecting disparate fragments when frag_mode = SIMPLE, a "bond"
        is assigned if interatomic distance is less than (this number) * sum of covalent radii. The
        value is then increased until all the fragments are connected (directly or indirectly). -*/
        options.add_double("INTERFRAGMENT_CONNECT", 1.8);
        /*- Tolerance for whether to reject a set of generated reference atoms due to collinearity -*/
        options.add_double("INTERFRAG_COLLINEAR_TOL", 0.01);

        /*- For now, this is a general maximum distance for the definition of H-bonds -*/
        options.add_double("H_BOND_CONNECT", 4.3);
        /*- Do only generate the internal coordinates and then stop? -*/
        options.add_bool("INTCOS_GENERATE_EXIT", false);
        /*- SUBSECTION Misc. -*/

        /*- Do test B matrix? -*/
        options.add_bool("TEST_B", false);
        /*- Do test derivative B matrix? -*/
        options.add_bool("TEST_DERIVATIVE_B", false);

        /*- Write the optimization history / state to disc -*/
        options.add_bool("WRITE_OPT_RESULT", false);
        /*- Save OptKing's internal classes for possible restart upon error -*/
        options.add_bool("SAVE_OPTIMIZATION", false);
        /*- Restart the optimization from optking's written history -*/
        options.add_double("OPT_RESTART", 0);
        /*- Write Optimization Trajectory -*/
        options.add_bool("WRITE_TRAJECTORY", false);
        /*- Should an xyz trajectory file be kept (useful for visualization)? -*/
        options.add_bool("PRINT_TRAJECTORY_XYZ_FILE", false);
        /*- Write the full history to disk. Produces a non validated OptimizationResult. -*/
        options.add_bool("WRITE_OPT_HISTORY", false);

    }
    if (name == "FINDIF" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs finite difference computations of energy derivative, with respect to nuclear
           displacements for geometry optimizations and vibrational frequency analyses, where the required analytical
           derivatives are not available. -*/

        /*- Number of points for finite-differences (3 or 5) -*/
        options.add_int("POINTS", 3);  // Can we error check integers?
        /*- Displacement size in au for finite-differences. -*/
        options.add_double("DISP_SIZE", 0.005);
        /*- Do write a gradient output file?  If so, the filename will end in
        .grad, and the prefix is determined by |globals__writer_file_label|
        (if set), or else by the name of the output file plus the name of
        the current molecule. -*/
        options.add_bool("GRADIENT_WRITE", false);
        /*- Do write a hessian output file?  If so, the filename will end in
        .hess, and the prefix is determined by |globals__writer_file_label|
        (if set), or else by the name of the output file plus the name of
        the current molecule. -*/
        options.add_bool("HESSIAN_WRITE", false);
        /*- Do write a file containing the normal modes in Molden format?
       If so, the filename will end in .molden_normal_modes, and the prefix is
       determined by |globals__writer_file_label| (if set), or else by the name
       of the output file plus the name of the current molecule. -*/
        options.add_bool("NORMAL_MODES_WRITE", false);
        /*- Do discount rotational degrees of freedom in a finite difference
        frequency calculation. Turned off at non-stationary geometries and
        in the presence of external perturbations. -*/
        options.add_bool("FD_PROJECT", true);
    }
    if (name == "OCC" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs orbital-optimized MPn and CC computations and conventional MPn computations. -*/

        /*- Maximum number of iterations to determine the amplitudes -*/
        options.add_int("CC_MAXITER", 50);
        /*- Maximum number of iterations to determine the orbitals -*/
        options.add_int("MO_MAXITER", 50);
        /*- Caching level for libdpd governing the storage of amplitudes,
        integrals, and intermediates in the CC procedure. A value of 0 retains
        no quantities in cache, while a level of 6 attempts to store all
        quantities in cache.  For particularly large calculations, a value of
        0 may help with certain types of memory problems.  The default is 2,
        which means that all four-index quantities with up to two virtual-orbital
        indices (e.g., $\langle ij | ab \rangle$ integrals) may be held in the cache. -*/
        options.add_int("CACHELEVEL", 2);
        /*- Minimum number of error vectors stored for DIIS extrapolation -*/
        options.add_int("DIIS_MIN_VECS", 2);
        /*- Maximum number of error vectors stored for DIIS extrapolation -*/
        options.add_int("DIIS_MAX_VECS", 6);
        /*- Cutoff value for numerical procedures -*/
        options.add_int("CUTOFF", 14);
        /*- Maximum number of preconditioned conjugate gradient iterations.  -*/
        options.add_int("PCG_MAXITER", 30);
        /*- Maximum number of electron propagator iterations.  -*/
        options.add_int("EP_MAXITER", 30);
        /*- Convergence criterion for energy. See Table :ref:`Post-SCF
        Convergence <table:conv_corl>` for default convergence criteria for
        different calculation types. -*/
        options.add_double("E_CONVERGENCE", 1e-6);
        /*- Convergence criterion for amplitudes (residuals). -*/
        options.add_double("R_CONVERGENCE", 1e-5);
        /*- Convergence criterion for RMS orbital gradient. If this keyword is not
        set by the user, OCC will estimate and use a value required to achieve the
        desired |occ__e_convergence|. The listed default will be used for the default
        value of |occ__e_convergence|. -*/
        options.add_double("RMS_MOGRAD_CONVERGENCE", 1e-4);
        /*- Convergence criterion for maximum orbital gradient. If this keyword is not
        set by the user, OCC will estimate and use a value required to achieve the
        desired |occ__e_convergence|. The listed default will be used for the default
        value of |occ__e_convergence|. -*/
        options.add_double("MAX_MOGRAD_CONVERGENCE", 1e-4);
        /*- Maximum step size in orbital-optimization procedure -*/
        options.add_double("MO_STEP_MAX", 0.5);
        /*- Removed in 1.4. Will raise an error in 1.5. -*/
        options.add_double("LEVEL_SHIFT", 0.02);
        /*- Removed in 1.4. Will raise an error in 1.5. -*/
        options.add_double("MP2_OS_SCALE", 6.0 / 5.0);
        /*- Removed in 1.4. Will raise an error in 1.5. -*/
        options.add_double("MP2_SS_SCALE", 1.0 / 3.0);
        /*- Removed in 1.4. Will raise an error in 1.5. -*/
        options.add_double("MP2_SOS_SCALE", 1.3);
        /*- Removed in 1.4. Will raise an error in 1.5. -*/
        options.add_double("MP2_SOS_SCALE2", 1.2);
        /*- Removed in 1.4. Will raise an error in 1.5. -*/
        options.add_double("CEPA_OS_SCALE", 1.27);
        /*- Removed in 1.4. Will raise an error in 1.5. -*/
        options.add_double("CEPA_SS_SCALE", 1.13);
        /*- Removed in 1.4. Will raise an error in 1.5. -*/
        options.add_double("CEPA_SOS_SCALE", 1.3);
        /*- A custom scaling parameter for opposite-spin terms in OCC. The result goes to a CUSTOM SCS variable, exact name method-dependent. -*/
        options.add_double("OS_SCALE", 1);
        /*- A custom scaling parameter for same-spin terms in OCC. The result goes to a CUSTOM SCS variable, exact name method-dependent. -*/
        options.add_double("SS_SCALE", 1);
        /*- Scaling value for 3rd order energy correction (S. Grimme, Vol. 24, pp. 1529, J. Comput. Chem.) -*/
        options.add_double("E3_SCALE", 0.25);
        /*- Convergence criterion for residual vector of preconditioned conjugate gradient method. -*/
        options.add_double("PCG_CONVERGENCE", 1e-6);
        /*- Damping factor for the orbital gradient (Rendell et al., JCP, vol. 87, pp. 5976, 1987) -*/
        options.add_double("MOGRAD_DAMPING", 1.0);
        /*- mixing parameter for the REMP hybrid perturbation theory, A specifies the Moller-Plesset fraction -*/
        options.add_double("REMP_A", 0.15E0);

        /*- The algorithm for orthogonalization of MOs -*/
        options.add_str("ORTH_TYPE", "MGS", "GS MGS");
        /*- The optimization algorithm. Modified Steepest-Descent (MSD) takes a Newton-Raphson (NR) step
         with a crude approximation to diagonal elements of the MO Hessian. The ORB_RESP option obtains the orbital
         rotation parameters with a crude approximation to all elements of the MO Hessian. Additionally, for both
         methods a DIIS extrapolation will be performed with the DO_DIIS = TRUE option. -*/
        options.add_str("OPT_METHOD", "MSD", "MSD ORB_RESP");
        /*- The algorithm will be used for solving the orbital-response equations. The LINEQ option create the MO
          Hessian and solve the simultaneous linear equations with method choosen by the LINEQ_SOLVER option. The PCG
          option does not create the MO Hessian explicitly, instead it solves the simultaneous equations iteratively
          with the preconditioned conjugate gradient method. -*/
        options.add_str("ORB_RESP_SOLVER", "PCG", "PCG LINEQ");
        /*- Type of PCG beta parameter (Fletcher-Reeves or Polak-Ribiere). -*/
        options.add_str("PCG_BETA_TYPE", "FLETCHER_REEVES", "FLETCHER_REEVES POLAK_RIBIERE");
        /*- Type of the SCS method -*/
        options.add_str("SCS_TYPE", "SCS", "SCS SCSN SCSVDW SCSMI");
        /*- Type of the SOS method -*/
        options.add_str("SOS_TYPE", "SOS", "SOS SOSPI");
        /*- Type of the wavefunction. -*/
        options.add_str("WFN_TYPE", "OMP2", "OMP2 OMP3 OCEPA OMP2.5 REMP OREMP");
        /*- How to take care of the TPDM VVVV-block. The COMPUTE option means it will be computed via an IC/OOC
        algorithm. The DIRECT option (default) means it will not be computed and stored, instead its contribution will
        be directly added to Generalized-Fock Matrix. -*/
        options.add_str("TPDM_ABCD_TYPE", "DIRECT", "DIRECT COMPUTE");
        /*- CEPA type such as CEPA0, CEPA1 etc. currently we have only CEPA0. -*/
        options.add_str("CEPA_TYPE", "CEPA0", "CEPA0");
        /*- Controls the spin scaling set to current energy. This is set by Psi internally. !expert -*/
        options.add_str("SPIN_SCALE_TYPE", "NONE", "NONE CUSTOM SCS SCSN SCSVDW SOS SOSPI");

        /*- Do compute natural orbitals? -*/
        options.add_bool("NAT_ORBS", false);
        /*- Removed in 1.4. Will raise an error in 1.5. -*/
        options.add_bool("DO_LEVEL_SHIFT", true);
        /*- Do print OCC orbital energies? -*/
        options.add_bool("OCC_ORBS_PRINT", false);
        /*- Removed in 1.4. Will raise an error in 1.5. Pass the method name, like scs-mp2, to energy instead. -*/
        options.add_bool("DO_SCS", false);
        /*- Removed in 1.4. Will raise an error in 1.5. Pass the method name, like scs-mp2, to energy instead. -*/
        options.add_bool("DO_SOS", false);
        /*- Do write coefficient matrices to external files for direct reading MOs in a subsequent job? -*/
        options.add_bool("MO_WRITE", false);
        /*- Do read coefficient matrices from external files of a previous OMP2 or OMP3 computation? -*/
        options.add_bool("MO_READ", false);
        /*- Do apply DIIS extrapolation? -*/
        options.add_bool("DO_DIIS", true);
        /*- Do compute CC Lambda energy? In order to this option to be valid one should use "TPDM_ABCD_TYPE = COMPUTE"
        option. -*/
        options.add_bool("CCL_ENERGY", false);
        /*- Do compute OCC poles for ionization potentials? Only valid OMP2. -*/
        options.add_bool("IP_POLES", false);
        /*- Do compute OCC poles for electron affinities? Only valid for OMP2. -*/
        options.add_bool("EA_POLES", false);
        /*- Do compute EP-OCC poles for ionization potentials? Only valid OMP2. -*/
        options.add_bool("EP_IP_POLES", false);
        /*- Do compute EP-OCC poles for electron affinities? Only valid for OMP2. -*/
        options.add_bool("EP_EA_POLES", false);
        /*- Do compute occupied orbital energies based on extended Koopmans' theorem? -*/
        options.add_bool("EKT_IP", false);
        /*- Do compute virtual orbital energies based on extended Koopmans' theorem?  -*/
        options.add_bool("EKT_EA", false);
        /*- Do optimize the orbitals?  -*/
        options.add_bool("ORB_OPT", true);
        /*- Do consider orbital response contributions for PDMs and GFM?  -*/
        options.add_bool("RELAXED", true);
        /*- Do symmetrize the GFM and OPDM in the EKT computations?  -*/
        options.add_bool("SYMMETRIZE", true);
        /*- Do compute one electron properties?  -*/
        options.add_bool("OEPROP", false);
    }
    if (name == "DFOCC" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs density-fitted orbital-optimized MPn and CC computations and conventional MPn
           computations. -*/

        /*- Maximum number of iterations to determine the amplitudes -*/
        options.add_int("CC_MAXITER", 50);
        /*- Maximum number of iterations to determine the orbitals -*/
        options.add_int("MO_MAXITER", 100);
        /*- Maximum number of preconditioned conjugate gradient iterations.  -*/
        options.add_int("PCG_MAXITER", 50);
        /*- Number of vectors used in orbital DIIS -*/
        options.add_int("MO_DIIS_NUM_VECS", 6);
        /*- Minimum number of vectors used in amplitude DIIS -*/
        options.add_int("CC_DIIS_MIN_VECS", 2);
        /*- Maximum number of vectors used in amplitude DIIS -*/
        options.add_int("CC_DIIS_MAX_VECS", 6);
        /*- Cutoff value for DF integrals -*/
        options.add_int("INTEGRAL_CUTOFF", 9);
        /*- Cutoff value for numerical procedures -*/
        options.add_int("CUTOFF", 8);

        /*- Convergence criterion for energy. See Table :ref:`Post-SCF
        Convergence <table:conv_corl>` for default convergence criteria for
        different calculation types. -*/
        options.add_double("E_CONVERGENCE", 1e-6);
        /*- Convergence criterion for amplitudes (residuals). -*/
        options.add_double("R_CONVERGENCE", 1e-5);
        /*- Convergence criterion for RMS orbital gradient. If this keyword is not
        set by the user, DFOCC will estimate and use a value required to achieve the
        desired |dfocc__e_convergence|. The listed default will be used for the default
        value of |dfocc__e_convergence|. -*/
        options.add_double("RMS_MOGRAD_CONVERGENCE", 1e-4);
        /*- Convergence criterion for maximum orbital gradient. If this keyword is not
        set by the user, DFOCC will estimate and use a value required to achieve the
        desired |dfocc__e_convergence|. The listed default will be used for the default
        value of |dfocc__e_convergence|. -*/
        options.add_double("MAX_MOGRAD_CONVERGENCE", 1e-4);
        /*- Maximum step size in orbital-optimization procedure -*/
        options.add_double("MO_STEP_MAX", 0.5);
        /*- Level shift to aid convergence -*/
        options.add_double("LEVEL_SHIFT", 0.02);
        /*- MP2 opposite-spin scaling value -*/
        options.add_double("MP2_OS_SCALE", 6.0 / 5.0);
        /*- MP2 same-spin scaling value -*/
        options.add_double("MP2_SS_SCALE", 1.0 / 3.0);
        /*- MP2 Spin-opposite scaling (SOS) value -*/
        options.add_double("MP2_SOS_SCALE", 1.3);
        /*- Spin-opposite scaling (SOS) value for optimized-MP2 orbitals -*/
        options.add_double("MP2_SOS_SCALE2", 1.2);
        /*- Scaling value for 3rd order energy correction (S. Grimme, Vol. 24, pp. 1529, J. Comput. Chem.) -*/
        options.add_double("E3_SCALE", 0.25);
        /*- OO scaling factor used in MSD -*/
        options.add_double("OO_SCALE", 0.01);
        /*- Convergence criterion for residual vector of preconditioned conjugate gradient method.
        If this keyword is not set by the user, DFOCC will estimate and use a value required to achieve
        |dfocc__r_convergence| residual convergence. The listed default will be used for the default value
        of |dfocc__r_convergence|. -*/
        options.add_double("PCG_CONVERGENCE", 1e-7);
        /*- Regularization parameter -*/
        options.add_double("REG_PARAM", 0.4);
        /*- tolerance for Cholesky decomposition of the ERI tensor -*/
        options.add_double("CHOLESKY_TOLERANCE", 1.0e-4);

        /*- The solver will be used for simultaneous linear equations. -*/
        options.add_str("LINEQ_SOLVER", "CDGESV", "CDGESV FLIN POPLE");
        /*- The algorithm for orthogonalization of MOs -*/
        options.add_str("ORTH_TYPE", "MGS", "GS MGS");
        /*- The orbital optimization algorithm. Presently quasi-Newton-Raphson algorithm available with several Hessian
         * options. -*/
        options.add_str("OPT_METHOD", "QNR", "QNR");
        /*- The algorithm will be used for solving the orbital-response equations. The LINEQ option create the MO
          Hessian and solve the simultaneous linear equations with method choosen by the LINEQ_SOLVER option. The PCG
          option does not create the MO Hessian explicitly, instead it solves the simultaneous equations iteratively
          with the preconditioned conjugate gradient method. -*/
        options.add_str("ORB_RESP_SOLVER", "PCG", "PCG LINEQ");
        /*- Type of the MO Hessian matrix -*/
        options.add_str("HESS_TYPE", "HF", "APPROX_DIAG APPROX_DIAG_EKT APPROX_DIAG_HF HF");
        /*- Type of the SCS method -*/
        options.add_str("SCS_TYPE", "SCS", "SCS SCSN SCSVDW SCSMI");
        /*- Type of the SOS method -*/
        options.add_str("SOS_TYPE", "SOS", "SOS SOSPI");
        /*- Type of the wavefunction. -*/
        options.add_str("WFN_TYPE", "DF-OMP2",
                        "DF-OMP2 DF-OMP3 DF-OLCCD DF-OREMP DF-OMP2.5 DFGRAD DF-CCSD DF-CCD DF-CCSD(T) DF-CCSD(AT) QCHF");
        /*- mixing parameter for the REMP hybrid perturbation theory, A specifies the Moller-Plesset fraction -*/
        options.add_double("REMP_A", 0.15E0);
        /*- CEPA type such as CEPA0, CEPA1 etc. currently we have only CEPA0. -*/
        // options.add_str("CEPA_TYPE","CEPA(0)","CEPA(0)");
        /*- The algorithm that used for 4 index MO TEIs. -*/
        // options.add_str("CONV_TEI_TYPE","DIRECT","DIRECT DISK");
        /*- Type of PCG beta parameter (Fletcher-Reeves or Polak-Ribiere). -*/
        options.add_str("PCG_BETA_TYPE", "FLETCHER_REEVES", "FLETCHER_REEVES POLAK_RIBIERE");
        /*- The algorithm that used to handle mp2 amplitudes. The DIRECT option means compute amplitudes on the fly
         * whenever they are necessary. -*/
        options.add_str("MP2_AMP_TYPE", "DIRECT", "DIRECT CONV");
        /*- Type of the CCSD PPL term. -*/
        options.add_str("PPL_TYPE", "AUTO", "LOW_MEM HIGH_MEM CD AUTO");
        /*- The algorithm to handle (ia|bc) type integrals that used for (T) correction. -*/
        options.add_str("TRIPLES_IABC_TYPE", "DISK", "INCORE AUTO DIRECT DISK");

        /*- Do compute natural orbitals? -*/
        options.add_bool("NAT_ORBS", false);
        /*- Do apply level shifting? -*/
        options.add_bool("DO_LEVEL_SHIFT", true);
        /*- Do print OCC orbital energies? -*/
        options.add_bool("OCC_ORBS_PRINT", false);
        /*- Do perform spin-component-scaled OMP2 (SCS-OMP2)? In all computation, SCS-OMP2 energy is computed
         automatically. However, in order to perform geometry optimizations and frequency computations with SCS-OMP2,
         one needs to set 'DO_SCS' to true -*/
        options.add_bool("DO_SCS", false);
        /*- Do perform spin-opposite-scaled OMP2 (SOS-OMP2)? In all computation, SOS-OMP2 energy is computed
         automatically. However, in order to perform geometry optimizations and frequency computations with SOS-OMP2,
         one needs to set 'DO_SOS' to true -*/
        options.add_bool("DO_SOS", false);
        /*- Do apply DIIS extrapolation? -*/
        options.add_bool("DO_DIIS", true);
        /*- Do compute ionization potentials based on the extended Koopmans' theorem? -*/
        options.add_bool("EKT_IP", false);
        /*- Do optimize the orbitals?  -*/
        options.add_bool("ORB_OPT", true);
        /*- Do use regularized denominators?  -*/
        options.add_bool("REGULARIZATION", false);
        /*- Do read 3-index integrals from SCF files?  -*/
        options.add_bool("READ_SCF_3INDEX", true);
        /*- Do compute one electron properties?  -*/
        options.add_bool("OEPROP", false);
        /*- Do compute $\langle \hat{S}^2 \rangle$ for DF-OMP2/DF-MP2?  -*/
        options.add_bool("COMPUT_S2", false);
        /*- Do perform a QCHF computation?  -*/
        options.add_bool("QCHF", false);
        /*- Do solve lambda amplitude equations?  -*/
        options.add_bool("CC_LAMBDA", false);
        /*- Do write a MOLDEN output file?  If so, the filename will end in
        .molden, and the prefix is determined by |globals__writer_file_label|
        (if set), or else by the name of the output file plus the name of
        the current molecule. -*/
        options.add_bool("MOLDEN_WRITE", false);
        /*- Do Cholesky decomposition of the ERI tensor -*/
        options.add_bool("CHOLESKY", false);
    }
    if (name == "MRCC" || options.read_globals()) {
        /*- MODULEDESCRIPTION Interface to MRCC program written by Mih\ |a_acute|\ ly K\ |a_acute|\ llay. -*/

        /*- Sets the OMP_NUM_THREADS environment variable before calling MRCC.
            If the environment variable :envvar:`OMP_NUM_THREADS` is set prior to calling Psi4 then
            that value is used. When set, this option overrides everything. Be aware
            the ``-n`` command-line option described in section :ref:`sec:threading`
            does not affect MRCC.
            !expert -*/
        options.add_int("MRCC_OMP_NUM_THREADS", 1);

        /*- Convergence criterion for energy. See Table :ref:`Post-SCF
        Convergence <table:conv_corl>` for default convergence criteria for
        different calculation types. This becomes ``tol`` (option \#16) in
        fort.56. -*/
        options.add_double("E_CONVERGENCE", 1e-6);

        /*- Schwarz screening threshold. Mininum absolute value below which TEI are neglected. -*/
        options.add_double("INTS_TOLERANCE", 1.0E-12);

        /*- Maximum excitation level. This is used ONLY if it is explicitly set by the user.
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

        /*- Number of root in case of open shell system. This becomes ``ndoub`` (option \#13)
           int fort.56. -*/
        options.add_int("MRCC_NUM_DOUBLET_ROOTS", 0);

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
            and will only work with :py:func:`~psi4.driver.energy`.
            This becomes ``CC/CI`` (option \#5) in fort.56.
            See Table :ref:`MRCC_METHOD <table:mrcc__mrcc_method>` for details.
            !expert -*/
        options.add_int("MRCC_METHOD", 1);
    }
    if (name == "FNOCC" || options.read_globals()) {
        /*- Do time each cc diagram? -*/
        options.add_bool("CC_TIMINGS", false);
        /*- Convergence criterion for CC energy. See Table :ref:`Post-SCF
        Convergence <table:conv_corl>` for default convergence criteria for
        different calculation types.  Note that convergence is
        met only when |fnocc__e_convergence| and |fnocc__r_convergence|
        are satisfied. -*/
        options.add_double("E_CONVERGENCE", 1.0e-6);
        /*- Maximum number of iterations for Brueckner orbitals optimization -*/
        options.add_int("BRUECKNER_MAXITER", 20);
        /*- Convergence for the CC amplitudes.  Note that convergence is
            met only when |fnocc__e_convergence| and |fnocc__r_convergence|
            are satisfied. -*/
        options.add_double("R_CONVERGENCE", 1.0e-7);
        /*- Maximum number of CC iterations -*/
        options.add_int("MAXITER", 100);
        /*- Desired number of DIIS vectors -*/
        options.add_int("DIIS_MAX_VECS", 8);
        /*- Do use low memory option for triples contribution? Note that this
            option is enabled automatically if the memory requirements of the
            conventional algorithm would exceed the available resources.
            The low memory algorithm is faster in general and has been turned
            on by default starting September 2020. -*/
        options.add_bool("TRIPLES_LOW_MEMORY", true);
        /*- Do compute triples contribution? !expert -*/
        options.add_bool("COMPUTE_TRIPLES", true);
        /*- Do compute MP4 triples contribution? !expert -*/
        options.add_bool("COMPUTE_MP4_TRIPLES", false);
        /*- Do use MP2 NOs to truncate virtual space for QCISD/CCSD and (T)? -*/
        options.add_bool("NAT_ORBS", false);
        /*- Cutoff for occupation of MP2 virtual NOs in FNO-QCISD/CCSD(T).
            Virtual NOs with occupations less than |fnocc__occ_tolerance|
            will be discarded. This option is only used if |fnocc__nat_orbs| =
            true. -*/
        options.add_double("OCC_TOLERANCE", 1.0e-6);
        /*- Cutoff for occupation of MP2 virtual NOs in FNO-QCISD/CCSD(T).
            The number of virtual NOs is chosen so the occupation of the
            truncated virtual space is |fnocc__occ_percentage| percent of
            occupation of the original MP2 virtual space. This option is only
            used if |fnocc__nat_orbs| = true. This keyword overrides
            |fnocc__occ_tolerance|. -*/
        options.add_double("OCC_PERCENTAGE", 99.0);
        /*- An array containing the number of virtual natural orbitals per irrep
        (in Cotton order) so a user can specify the number of retained
        natural orbitals rather than determining them with |fnocc__occ_tolerance|.
        This keyword overrides |fnocc__occ_tolerance| and
        |fnocc__occ_percentage|. -*/
        options.add("ACTIVE_NAT_ORBS", new ArrayType());
        /*- Do SCS-MP2? -*/
        options.add_bool("SCS_MP2", false);
        /*- Do SCS-CCSD? -*/
        options.add_bool("SCS_CCSD", false);
        /*- Do SCS-CEPA? Note that the scaling factors will be identical
        to those for SCS-CCSD. -*/
        options.add_bool("SCS_CEPA", false);
        /*- Opposite-spin scaling factor for SCS-MP2 -*/
        options.add_double("MP2_SCALE_OS", 1.20);
        /*- Same-spin scaling factor for SCS-MP2 -*/
        options.add_double("MP2_SCALE_SS", 1.0 / 3.0);
        /*- Oppposite-spin scaling factor for SCS-CCSD -*/
        options.add_double("CC_SCALE_OS", 1.27);
        /*- Same-spin scaling factor for SCS-CCSD -*/
        options.add_double("CC_SCALE_SS", 1.13);
        /*- do only evaluate mp2 energy? !expert -*/
        options.add_bool("RUN_MP2", false);
        /*- do only evaluate mp3 energy? !expert -*/
        options.add_bool("RUN_MP3", false);
        /*- do only evaluate mp4 energy? !expert -*/
        options.add_bool("RUN_MP4", false);
        /*- do ccsd rather than qcisd? !expert -*/
        options.add_bool("RUN_CCSD", false);
        /*- Use 3-index integrals to generate 4-index ERI's?
        This keyword is used for testing purposes only.  Changing its
        value will have no effect on the computation.  !expert -*/
        options.add_bool("USE_DF_INTS", false);

        /*- Do use density fitting or cholesky decomposition in CC? This
        keyword is used internally by the driver. Changing its value
        will have no effect on the computation. -*/
        options.add_bool("DFCC", false);
        /*- Auxilliary basis for df-ccsd(t). -*/
        options.add_str("DF_BASIS_CC", "");
        /*- tolerance for Cholesky decomposition of the ERI tensor -*/
        options.add_double("CHOLESKY_TOLERANCE", 1.0e-4);

        /*- Is this a CEPA job? This parameter is used internally
        by the pythond driver.  Changing its value won't have any
        effect on the procedure. !expert -*/
        options.add_bool("RUN_CEPA", false);
        /*- Which coupled-pair method is called?  This parameter is
        used internally by the python driver.  Changing its value
        won't have any effect on the procedure. !expert -*/
        options.add_str("CEPA_LEVEL", "CEPA(0)");
        /*- Compute the dipole moment? Note that dipole moments
        are only available in the FNOCC module for the ACPF,
        AQCC, CISD, and CEPA(0) methods. -*/
        options.add_bool("DIPMOM", false);
        /*- Flag to exclude singly excited configurations from a
        coupled-pair computation.  -*/
        options.add_bool("CEPA_NO_SINGLES", false);
    }
    if (name == "THERMO" || options.read_globals()) {
        /*- Temperature in Kelvin for thermodynamic analysis. Note that 273.15
        is the value for IUPAC STP. -*/
        options.add_double("T", 298.15);
        /*- Pressure in Pascal for thermodynamic analysis. Note that 100000.
        is the value for IUPAC STP. -*/
        options.add_double("P", 101325);
        /*- Rotational symmetry number for thermodynamic analysis. Default is set
        from the full point group (e.g., Td for methane) as opposed to the computational
        point group (e.g., C2v for methane). Default takes into account symmetry
        reduction through asymmetric isotopic substitution and is unaffected by
        user-set symmetry on molecule, so this option is the sole way to influence
        the symmetry-dependent aspects of the thermodynamic analysis.
        Note that this factor is handled differently among quantum chemistry software. -*/
        options.add_int("ROTATIONAL_SYMMETRY_NUMBER", 1);
    }
    if (name == "CFOUR" || options.read_globals()) {
        /*- MODULEDESCRIPTION Interface to CFOUR program written by Stanton and Gauss.
        Keyword descriptions taken from the
        `CFOUR Website <http://slater.chemie.uni-mainz.de/cfour/index.php?n=Main.ListOfKeywordsInAlphabeticalOrder>`_
        and extended by interface comments. -*/

        /*- SUBSECTION Psi4 Control of CFOUR -*/

        /*- Sets the OMP_NUM_THREADS environment variable before calling CFOUR.
            If the environment variable :envvar:`OMP_NUM_THREADS` is set prior to calling Psi4 then
            that value is used. When set, this option overrides everything. Be aware
            the ``-n`` command-line option described in section :ref:`sec:threading`
            does not affect CFOUR.
            !expert -*/
        options.add_int("CFOUR_OMP_NUM_THREADS", 1);

        /*- Do translate set Psi4 options to their cfour counterparts. -*/
        options.add_bool("TRANSLATE_PSI4", true);

        /*- SUBSECTION CFOUR Internal -*/

        /*- Specifies the way the :math:`\langle ab||cd \rangle` molecular orbital
        integrals are handled in post-MP2 calculations. STANDARD (= 0) uses
        directly the corresponding MO integrals and thus results in an
        algorithm which in particular for large-scale calculations results
        in excessive use of disk space (storage of all :math:`\langle ab||cd\rangle`
        integrals. AOBASIS (=2) uses an AO-based algorithm to evaluate all
        terms involving the :math:`\langle ab||cd\rangle` integrals and
        significantly reduces the amount of disk storage.  The use of
        ABCDTYPE=AOBASIS is strongly recommended for all CC calculations up
        to CCSD(T) and has been implemented for energy, gradient,
        second-derivative, and excitation energy calculations. -*/
        options.add_str("CFOUR_ABCDTYPE", "STANDARD", "STANDARD AOBASIS");

        /*- Specifies the active orbitals used in a TCSCF calculation and
        has to be used in combination with the keyword |cfour__cfour_core_orbitals|. The
        active orbitals are specified by either NIRREP or 2*NIRREP integers
        specifying the number of active orbitals of each symmetry type,
        where NIRREP is the number of irreducible representations in the
        computational point group. If there are no orbitals of a particular
        symmetry type a zero must be entered. For more information and an
        example see |cfour__cfour_occupation| . -*/
        options.add("CFOUR_ACTIVE_ORBI", new ArrayType());

        /*- Specifies treatment of anharmonic effects by calculating cubic
        and/or quartic force fields. VIBROT (=3) requests calculation
        of only those cubic constants of the form :math:`\phi_{nij}`, where n is a
        totally symmetric coordinate. These are sufficient to determine the
        vibration-rotation interaction constants needed to calculate
        vibrational corrections to rotational constants, but are *not*
        sufficient to generate the corresponding cubic constants of
        isotopologs that have a lower point-group symmetry (*i.e.* HOD
        isotopolog of water). VPT2 (=1, note that the old value
        CUBIC can be still used and is equivalent to VPT2)
        generates all cubic constants and all quartic constants apart from
        those of the form :math:`\phi_{ijkl}`, which is enough for: 1) generation of
        cubic constants of isotopologs (see manual entries associated with
        anharmonic calculations for an example); 2) calculation of
        vibrational energy levels with VPT2. This keyword also directs the
        program to analyze resonances and calculate intensities of one- and
        two-quantum transitions. FULLQUARTIC (=2) (not part of the
        public release) is largely self-explanatory; it directs the program
        to calculate all quartic constants. This is sufficient (but this has
        not been implemented) to generate the full quartic force field of
        all isotopologs. -*/
        options.add_str("CFOUR_ANHARMONIC", "OFF", "CUBIC VPT2 FULLQUARTIC VIBROT OFF");

        /*- Specifies which algorithm is used for |cfour__cfour_anharmonic| =VIBROT, VPT2,
        and FULLQUARTIC calculations. If STANDARD (=0) is chosen,
        then simply invoking ``xcfour`` will cause a complete job to be run with
        all second-derivative calculations being done in series. If
        PARALLEL (=1), then the job stops after the second-derivative
        calculation at the reference geometry and generates out all input
        geometries for the remaining calculation. These can be then
        processed in "parallel" (currently not recommended).  Note that it
        is recommended to carry out all calculations with
        PARALLEL, even when the actual calculation is carried
        out in a sequential mode. -*/
        options.add_str("CFOUR_ANH_ALGORITHM", "STANDARD", "STANDARD PARALLEL");

        /*- Specifies whether the anharmonic force field is calculated using
        analytic gradients (=FIRST) or analytic Hessians (=SECOND). -*/
        options.add_str("CFOUR_ANH_DERIVATIVES", "SECOND", "FIRST SECOND");

        /*- Controls the stepsize used in anharmonic force field
        calculations. The value is specified in reduced normal coordinates,
        which are dimensionless. The actual stepsize used in the calculation
        is :math:`\times 10^6` the integer value specified. -*/
        options.add_int("CFOUR_ANH_STEPSIZE", 50000);

        /*- Specifies whether non-abelian symmetry is to be exploited in
        determining displacements for |cfour__cfour_anharmonic| =VIBROT or VPT2 calculations. If
        set to NONABELIAN (=0), maximum advantage will be taken of symmetry
        and the full set of cubic force constants will be generated from a
        skeleton set by application of the totally symmetric projection
        operator. If set to ABELIAN (=1), only the operations of the
        abelian subgroup will be exploited. Note: It is important to point out
        that the symmetrization currently works only for cubic constants.
        Therefore, if you require quartic force constants (for frequency
        calculations), you *must* use the ABELIAN option. Moreover, the latter
        work for only asymmetric tops and linear molecules. -*/
        options.add_str("CFOUR_ANH_SYMMETRY", "ABELIAN", "ABELIAN NONABELIAN");

        /*- Can be used to control the algorithm used by CFOUR when terms
        involving :math:`\langle ab||cd\rangle` molecular orbital integrals are
        calculated in the atomic orbital basis (see |cfour__cfour_abcdtype|).
        MULTIPASS (= 0) uses an approach where the AO integral file is read
        a number of times in order to ensure maximal vectorization and is
        usually the optimal strategy on supercomputers; SINGLEPASS (= 1)
        determines the contributions with only a single pass through the AO
        integrals, but at the cost of significantly reduced vectorization.
        In general, however, SINGLEPASS is definitely preferable on
        workstations with RISC architectures. (Default : MULTIPASS on all
        64-bit machines (e.g., CRAY-YMP) ; SINGLEPASS on all 32-bit machines
        (e.g., IBM-RS6000, HP-735, SGI-Indigo, DEC alphastations)).
        SPARSE_AO (=2) uses a sparse matrix algorithm which first rearranges
        the integral matrix in order to get "well-occupied" and "very
        sparse" blocks. "Well-occupied" blocks will be multiplied by matrix
        multiplication while in "very sparse" blocks only the non-zero
        elements are considered. The computational time is further reduced
        using symmetrized and anti-symmetrized integral and amplitude
        matrices in the multiplication. Substantial saving is assumed if
        SPARSE_AO (=2) is used. -*/
        options.add_str("CFOUR_AO_LADDERS", "SINGLEPASS", "MULTIPASS SINGLEPASS");

        /*- Experimental Use!  ON (=1) requests and averaged SCF over two
        states. So far only implemented for degenerate doublet-Pi states and
        used in conjunction with SOPERT. -*/
        options.add_bool("CFOUR_AV_SCF", false);

        /*- Specifies the AO basis used in the calculation. One can either
        specify a basis known to CFOUR or via BASIS=SPECIAL (=0) requests an
        arbitrary basis (see non-standard basis-set input). However, the
        latter must be available in the supplied GENBAS file. As standard
        basis sets, currently the following are available.
        **Psi4 Interface:** Recommended to use instead |mints__basis| for
        larger basis set selection and greater flexibility. When |mints__basis|
        used, |cfour__cfour_spherical| is set appropriately. -*/
        options.add_str("CFOUR_BASIS", "SPECIAL",
                        "STO-3G 3-21G 4-31G 6-31G 6-31G* 6-31G** 6-311G 6-311G* 6-311G** DZ DZP TZ TZP TZ2P PVDZ PVTZ "
                        "PVQZ PV5Z PV6Z PCVDZ PCVTZ PCVQZ PCV5Z PCV6Z AUG-PVDZ AUG-PVTZ AUG-PVTZ AUG-PVQZ AUG-PV5Z "
                        "AUG-PV6Z D-AUG-PVDZ D-AUG-PVTZ D-AUG-PVQZ D-AUG-PV5Z D-AUG-PV6Z cc-pVDZ cc-pVTZ cc-pVQZ "
                        "cc-pV5Z cc-pV6Z cc-pCVDZ cc-pCVTZ cc-pCVQZ cc-pCV5Z cc-pCV6Z PWCVDZ PWCVTZ PWCVQZ PWCV5Z "
                        "PWCV6Z PwCVDZ PwCVTZ PwCVQZ PwCV5Z PwCV6Z svp dzp tzp tzp2p qz2p pz3d2f 13s9p4d3f WMR ANO0 "
                        "ANO1 ANO2 EVEN_TEMPERED SPECIAL");

        /*- experimental use -*/
        // BREIT

        /*- Specifies the convergence criterion in Brueckner based CC
        calculations. The calculation is considered to be converged when the
        absolute value of largest single excitation amplitudes falls below
        $10^N$, where NNN is the value associated with the keyword. -*/
        options.add_int("CFOUR_BRUCK_CONV", 4);

        /*- Specifies whether Brueckner orbitals are to be determined for
        the specified CC method. OFF(=0) Brueckner orbitals are not to be
        determined, ON (=1) they are to be determined. -*/
        options.add_bool("CFOUR_BRUECKNER", false);

        // experimental use
        // BUFFERSIZE

        /*- Defines the level of calculation to be performed.
        **Psi4 Interface:** Keyword set from argument of computation
        command: CCSD if ``energy('c4-ccsd')``, *etc.* See :ref:`Energy
        (CFOUR) <table:energy_cfour>` and :ref:`Gradient (CFOUR)
        <table:grad_cfour>`. for all available. -*/
        options.add_str("CFOUR_CALC_LEVEL", "SCF",
                        "SCF HF MBPT(2) MP2 MBPT(3) MP3 SDQ-MBPT(4) SDQ-MP4 MBPT(4) MP4 CCD CCSD CCSD(T) CCSDT-1 "
                        "CCSDT-1b CCSDT-2 CCSDT-3 CCSDT-4 CCSDT CC2 CC3 QCISD QCISD(T) CID CISD UCC(4) B-CCD");

        /*- The number of records held in the i/o cache used by the post-SCF
        programs. The maximum number of records which can be held is 100. -*/
        options.add_int("CFOUR_CACHE_RECS", 10);

        // CCORBOPT
        // experimental use

        /*- Specifies the convergence criterion for the CC amplitude
        equations. The amplitudes are considered to be converged when the
        maximum of all (absolute) changes in the amplitudes is less than
        $10^N$, where $N$ is the value associated with the keyword. -*/
        options.add_int("CFOUR_CC_CONV", 7);

        /*- Specifies the maximum number of expansion vectors used in the
        iterative subspace to enhance convergence in the solution of the CC
        equations. -*/
        options.add_int("CFOUR_CC_EXPORDER", 5);

        /*- Specifies the type of convergence acceleration used to solve the
        CC equations. RLE (=0) uses the RLE methods of Purvis and Bartlett,
        DIIS (=1) uses the DIIS approach by Pulay, NOJACOBI (=2) uses RLE
        with continuous extrapolation, OFF (=3) uses no convergence
        acceleration. In general, DIIS provides the best results and is
        recommended, while OFF often results in poor convergence and thus
        cannot be recommended. -*/
        options.add_str("CFOUR_CC_EXTRAPOLATION", "DIIS", "RLE DIIS NOJACOBI OFF");

        /*- Specifies the maximum number of iterations in solving the CC
        amplitude equations. -*/
        options.add_int("CFOUR_CC_MAXCYC", 50);

        /*- Specifies which CC program is used. The available options are
        VCC (=0), ECC (=1), MRCC (=2), and EXTERNAL (=3). The default for
        all calculations is currently VCC which requests usage of ``xvcc``, but
        in many cases (e.g., for CCSD and CCSD(T)) ECC should be preferred
        due to the better performance of ``xecc`` (available currently for CCSD,
        CCSD+T, CCSD(T), and closed-shell CCSDT-n, CC3, and CCSDT). MRCC and
        External are intended for CC programs outside the CFOUR suite, e.g.,
        the general CC module mrcc written by M. Kallay (Budapest, Hungary).
        Default: VCC Note: Using the option ECC is not recommended for ROHF
        gradients. That is, if you are doing a geometry optimization with
        ROHF as your reference wave function then it is safe to use the
        option VCC.
        **Psi4 Interface:** Keyword set according to best practice for the
        computational method |cfour__cfour_calc_level|, reference
        |cfour__cfour_reference| (NYI) and derivative level
        |cfour__cfour_deriv_level| according to Table :ref:`Best Practices
        <table:cfour_cc_program>` when method specified by argument to
        computation command (*e.g.*, when ``energy('c4-ccsd')`` requested
        but not when ``energy('cfour')`` requested). Value can always be set
        explicitly. -*/
        options.add_str("CFOUR_CC_PROGRAM", "VCC", "VCC ECC NCC MRCC EXTERNAL");

        /*- Specifies the molecular charge.
        **Psi4 Interface:** Keyword set from active molecule. -*/
        options.add_int("CFOUR_CHARGE", 0);

        /*- Specifies the convergence threshold as :math:`10^{-N}` for CIS
        calculations. -*/
        options.add_int("CFOUR_CIS_CONV", 5);

        // COMM_SIZE
        // experimental use

        /*- Signifies that one or more "continuum" orbitals should be added
        to the calculation. VIRTUAL and DVIRTUAL specify one or two orbital
        which should be initially unoccupied (in the SCF calculation), while
        OCCUPIED and DOCCUPIED specify one or two orbitals which should be
        initially occupied. -*/
        options.add_str("CFOUR_CONTINUUM", "NONE", "NONE VIRTUAL DVIRTUAL OCCUPIED DOCCUPIED");

        /*- Specifies the contraction scheme used by the integral and
        integral derivative program. SEGMENTED (=0) uses a segmented
        contraction scheme; GENERAL (=1) uses a general contraction scheme,
        and UNCONTRACTED (=2) uses the corresponding uncontracted sets. Note
        that even for truly segmented basis sets, the integral programs run
        significantly faster in the GENERAL mode. -*/
        options.add_str("CFOUR_CONTRACTION", "GENERAL", "SEGMENTED GENERAL UNCONTRACTED");

        /*- Identical to |cfour__cfour_geo_conv|. -*/
        options.add_int("CFOUR_CONVERGENCE", 4);

        /*- Specifies the type of coordinates used in the input file ZMAT.
        Value INTERNAL (=0) means that the geometry is supplied in the
        usual Z-matrix format, while CARTESIAN (=1) means that the geometry
        is given in Cartesian coordinates. A third option is XYZINT (=2) for
        which a Z-matrix connectivity is defined, but with values of the
        internal coordinates defined implicitly by supplying Cartesian
        coordinates. Note that geometry optimizations are currently only
        possible for INTERNAL and XYZ2INT.
        **Psi4 Interface:** Keyword set from active molecule, always CARTESIAN.
        Above restrictions on geometry optimizations no longer apply. -*/
        options.add_str("CFOUR_COORDINATES", "INTERNAL", "INTERNAL CARTESIAN XYZINT");

        /*- Specifies the core orbitals used in a TCSCF calculation and has
        to be used in combination with the keyword |cfour__cfour_active_orbi|. The core
        orbitals are specified by either NIRREP or 2*NIRREP integers
        specifying the number of core orbitals of each symmetry type, where
        NIRREP is the number of irreducible representations in the
        computational point group. If there are no orbitals of a particular
        symmetry type a zero must be entered. For more information and an
        example see |cfour__cfour_occupation|. -*/
        options.add("CFOUR_CORE_ORBITALS", new ArrayType());

        /*- Specifies the convergence criterion for the iterative solution
        of the CPHF and Z-vector equations. The solutions are considered to
        be converged when the residual norm of the error vector falls below
        $10^N$. -*/
        options.add_int("CFOUR_CPHF_CONVER", 12);

        /*- Specifies the maximum number of cycles allowed for the solution
        of the CPHF- and/or Z-vector equations. -*/
        options.add_int("CFOUR_CPHF_MAXCYC", 64);

        /*- Specifies whether or not Hessian matrix is transformed
        (nonlinearly) to curvilinear internal coordinates. A value of 0 (or
        OFF) turns the transformation off if the analytic force constants
        are not available, while it is always performed if CURVILINEAR=1 (or
        ON). Values higher than 1 (or NO) unconditionally turn the
        transformation off.(Default: ON if analytic Hessian is available,
        OFF otherwise). -*/
        options.add_bool("CFOUR_CURVILINEAR", true);

        /*- Specifies whether the diagonal Born-Oppenheimer correction
        (DBOC) to the energy is evaluated (ON =1) or not (OFF =0). DBOC
        calculations are currently only available for HF-SCF and CCSD using
        RHF or UHF reference functions. -*/
        options.add_bool("CFOUR_DBOC", false);

        /*- Specifies whether the Dipole Coupling Tensor (DCT) is calculated
        (ON =1) or not (OFF =0). -*/
        options.add_bool("CFOUR_DCT", false);

        /*- Specifies whether or not energy derivatives are to be calculated
        and if so whether first or second derivatives are computed. ZERO (=
        0) derivatives are not calculated, FIRST (=1) first derivatives are
        calculated, SECOND (=2) second derivatives are calculated.  Note
        that this keyword usually needs not be set in any calculation since
        it is automatically set if the appropriate other options in the
        CFOUR namelist are turned on.
        **Psi4 Interface:** Keyword set from type of computation command:
        ZERO if :py:func:`~psi4.driver.energy`, FIRST if
        :py:func:`~psi4.driver.gradient` or :py:func:`~psi4.driver.optimize`,
        *etc.* -*/
        options.add_str("CFOUR_DERIV_LEVEL", "ZERO", "ZERO FIRST SECOND");

        /*- Specifies whether orbital-relaxed (RELAXED =0) or
        orbital-unrelaxed (UNRELAXED =1) derivatives are computed in the CC
        calculation. -*/
        options.add_str("CFOUR_DIFF_TYPE", "RELAXED", "RELAXED UNRELAXED");

        // DIRECT
        // experimental use

        // DIAG_MRCC
        // experimental use

        /*- Specifies which molecular orbitals will be dropped from the
        post-SCF calculation. The orbitals are numbered in ascending order
        from the most stable (negative energy) to the most unstable (largest
        positive energy). Individual orbitals must be separated with a dash,
        while x>y means orbitals x through y inclusive. For example, the
        string ``1>10-55-58>64``, would result in orbitals
        1,2,3,4,5,6,7,8,9,10,55,58,59,60,61,62,63 and 64 being dropped. For
        UHF calculations, the appropriate orbitals are deleted for both spin
        cases. No dropped virtual MOs are currently allowed for gradient or
        property calculations.
        **Psi4 Interface:** The array above is specified in PSI as
        (white space tolerant) [1,2,3,4,5,6,7,8,9,10,55,58,59,60,61,62,63,64]. -*/
        options.add("CFOUR_DROPMO", new ArrayType());

        // EA_CALC
        // experimental use

        // EA_SYM
        // experimental use

        /*- Specifies whether effective core potentials (pseudopotentials)
        are used (ON, =1) or not (OFF, =0). -*/
        options.add_bool("CFOUR_ECP", false);

        /*- Specifies which eigenvector of the totally symmetric part of the
        block-factored Hessian is to be followed uphill in a transition
        state search. Eigenvectors are indexed by their eigenvalues -- the
        lowest eigenvalue is 1, the next lowest is 2, etc. The default is 1,
        which should always be used if you are not looking for a specific
        transition state which you know corresponds to motion along a
        different mode. In the future, relatively sophisticated generation
        of a guessed eigenvector will be implemented, but this is the way
        things are for now. Of course, this keyword has no
        meaning if |cfour__cfour_method| is not set to TS. -*/
        options.add_int("CFOUR_EIGENVECTOR", 1);

        /*- Experimental use, ON = 1 requests the evaluation of electrical
        anharmonicities -*/
        options.add_bool("CFOUR_EL_ANHARM", false);

        /*- Specifies the threshold used in converging CC-LR/EOM-CC
        calculations. The iterative diagonalization is continued until the
        RMS residual falls below $10^{-N}$ with $N$ as the value specified with
        this keyword. -*/
        options.add_int("CFOUR_ESTATE_CONV", 5);

        // EOM_NSING
        // experimental use

        // EOM_NTRIP
        // experimental use

        /*- Controls whether non-iterative triples corrections are applied
        after various types of EOM-CCSD calculation. Works with |cfour__cfour_excite| set to EOMIP, might
        work with EOMEE, certainly doesn't work with EOMEA. Use with great
        caution, preferably after having a few drinks. -*/
        options.add_bool("CFOUR_EOM_NONIT", false);

        ///*- For experimental use only. Selects the iterative diagonalization
        // algorithm for the EOMEE calculations. If set to DAVIDSON, the
        // general modified Davidson technique is used. If set to MULTIROOT, a
        // multi-root Davidson approach is invoked that evaluates all roots of
        // a symmetry block simultaneously. This approach is much more stable
        // if the roots are energetically close to each other. -*/
        // options.add_str("CFOUR_EOM_NSTATES", "", "DAVIDSON MULTIROOT");

        ///*- Selects the excited state the EOMEE properties are calculated
        // for. Only valid if EOM_NSTATES = MULTIROOT is set. It always refers
        // to the corresponding state of the last symmetry block considered. -*/
        // options.add_int("CFOUR_EOM_PROPSTA");

        // EOMFOLLOW
        // experimental use

        // ESTATE_DIAG
        // experimental use

        // ESTATE_LOCK
        // experimental use

        /*- The maximum number of expansion vectors used in the solution of
        EOMCC equations (Default: 20, hard-coded to 4 in triples
        calculations) -*/
        options.add_int("CFOUR_ESTATE_MAXCYC", 20);

        /*- This keyword applies only to EOM-CC calculations and specifies
        whether any excited or ionized state one-electron properties are to
        be calculated. Proper use of this keyword requires a relatively
        advanced knowledge of quantum chemistry and the available options
        are discussed here. The options are: OFF (=0) [no properties or
        transition moments are calculated]; EXPECTATION (=1) [transition
        moments and dipole strengths are calculated along with selected
        one-electron properties which are evaluated as expectation values];
        UNRELAXED (=2) [selected one-electron properties are calculated in
        an approximation that neglects relaxation of molecular orbitals];
        RESPONSE (=3) [selected one-electron properties are calculated as
        analytic first derivatives of the energy]. Except for EOMCC
        calculations on two-electron systems (which are exact), properties
        obtained by the three approaches will not be equivalent. The default
        value for this keyword is slightly complicated. For TDA
        calculations, the default is EXPECTATION since the evaluation of
        transition moments involves only a negligible amount of additional
        computation relative to the evaluation of the excitation energies.
        For EOMCC, the default is OFF since evaluation of any transition
        moments or properties requires approximately twice the computational
        time. Transition moments and dipole strengths are evaluated by
        default for all values of ESTATE_PROP other than OFF. -*/
        options.add_str("CFOUR_ESTATE_PROP", "", "OFF EXPECTATION UNRELAXED RESPONSE");

        /*- Specifies the number of excited states which are to be
        determined in each irreducible representation of the computational
        subgroup. The program attempts to find all of the lowest roots, but
        this is not guaranteed because the eigenvalue problem is not solved
        by direct matrix diagonalization, but rather by an iterative
        (modified Davidson) algorithm. For excited state gradient
        calculations, only one root (clearly) is used. In such a case, one
        and only one non-zero entry in the string can be used, and this
        value is usually set to one (*i.e.* 0/1/0/0). (However
        sometimes one wants to calculate the gradient for, say, the second
        root of a given symmetry, and in such a case, one could use
        0/2/0/0. What happens is that both roots are calculated,
        but only the second one is used in the subsequent density matrix and
        gradient calculation.) The format used for this keyword is identical
        to that used in |cfour__cfour_occupation|. For example, for a
        computational subgroup having four symmetry species, the string
        3/1/0/2 specifies that 6 total roots should be searched
        for, three in the first block, one in the second block, and two in
        the fourth block. It is also important to note that the ``%excite*``
        input, if present, takes precedence over this keyword.
        Default: All zeros.
        **Psi4 Interface:** The array above is specified in PSI as
        (white space tolerant) [3,1,0,2]. -*/
        options.add("CFOUR_ESTATE_SYM", new ArrayType());

        /*- Specifies whether just the excitation energies (OFF, =0) or in
        addition transition moments (EXPECTATION, =1) are calculated. Note
        that this keyword should not be used in excited-state calculations
        involving analytic gradients and that transition moments are
        essentially only available for EOM-CCSD/CCSD-LR. -*/
        options.add_str("CFOUR_ESTATE_TRANS", "OFF", "OFF EXPECTATION");

        /*- Tells the program, in the course of a geometry optimization, to
        calculate the Hessian explicitly every N cycles. 0 means never
        calculated explicitly.
        **Psi4 Interface:** Geometry optimizations run through PSI (except in
        sandwich mode) use PSI's optimizer and so this keyword has no effect.
        Use :ref:`optking <apdx:optking>` keywords instead,
        particularly |optking__full_hess_every|. -*/
        options.add_int("CFOUR_EVAL_HESS", 0);

        /*- Specifies in CC calculations using mrcc the excitation level if
        the calculation level has been chosen as CC(n), CI(n), or CCn(n). -*/
        options.add_int("CFOUR_EXCITATION", 0);

        /*- Specifies the type of EOM-CC/LR-CC treatment to be performed.
        Available options are NONE (=0), EOMEE (=3, the EOM-CC/CC-LR
        approach for the treatment of excited states), EOMIP (=4, the
        EOM-CC/CC-LR approach for the treatment of ionized states), EOMEA
        (=7, the EOM-CC/CC-LR approach for the treatment of
        electron-attached states). -*/
        options.add_str("CFOUR_EXCITE", "NONE", "NONE EOMEE EOMIP EOMEA");

        /*- Specifies the strength of a Fermi-Contact perturbation as
        required for finite-field calculations of spin densities and the FC
        contributions to indirect spin-spin coupling constants. The value
        must be specified as an integer and the FC strength used by the
        program will be the value of the keyword $\times 10^{-6}$. The atom for
        which the FC perturbation is switched on is specified in the ZMAT
        file after the CFOUR command line and potential basis set input, as
        follows

        %spin density
         N

        with N as the number of atom (in (X5,I3) format) in the order they
        are written by JODA to the MOL file. Be aware that for some atoms,
        the calculation has to be run in lower symmetry or even without
        symmetry. (Default : 0) -*/
        options.add_int("CFOUR_FC_FIELD", 0);

        /*- Specifies the algorithm used to compute the harmonic force
        constants in finite-difference calculations.GRADONLY (=0) evaluates
        the force constants and dipole moment derivatives by numerical
        differentiation of analytic gradients; ENERONLY (=1) evaluates the
        force constants by second differences of energies (dipole moment
        derivatives are not evaluated); while MIXED (=2) evaluates 1x1
        blocks of symmetry-blocked force constants by second differences pf
        energies and all other elements by first differences of gradients.
        the GRADONLY and MIXED approaches may, of course, only be used hwen
        using computational methods for which analytic gradients are
        available. -*/
        options.add_str("CFOUR_FD_CALCTYPE", "GRADONLY", "GRADONLY ENERONLY MIXED");

        /*- Requests that only vibrational frequencies of certain symmetry
        types are evaluated in a VIBRATION=FINDIF calculation. The numbers
        of the irreducible representations for which vibrational analysis is
        to be performed are separated by slashes. For example,
        FD_IRREP=1/3/4 means compute the frequencies of modes transforming
        as the first, third, and fourth irreducible representations. If a
        symmetry is specified for which there are no vibrational modes, the
        program will terminate. The labels of the irreducible
        representations for this keyword are not usually the same as those
        used in the rest of the calculation. Moreover, for some point
        groups, for example, those of linear molecules, the two sets of
        labels refer to different subgroups. There is as yet no
        straightforward way to determine what they will be without starting
        a calculation. If one runs the ``xjoda`` and then the ``xsymcor``
        executables, the relevant irreducible representations will be
        listed. If all vibrational frequencies are desired, this keyword
        need not be included.  Default : compute vibrational frequencies for
        all irreducible representations -*/
        options.add("CFOUR_FD_IRREPS", new ArrayType());

        /*- Specifies whether or not rotational degrees of freedoms are
        projected out from the symmetry-adapted coordinates in a finite
        difference calculations. ON (=0) uses rotationally projected
        coordinates, while OFF (=1) retains the rotational degrees of
        freedom. At a stationary point on the potential energy surface, both
        options will give equivalent harmonic force fields, but OFF should
        be used at non-stationary points. -*/
        options.add_str("CFOUR_FD_PROJECT", "ON", "ON OFF");

        /*- Specifies the step length in mass-weighted coordinates (in
        :math:`10^{-4} amu^{1/2} bohr` ) used in generating the force constant matrix
        by finite difference of Cartesian gradients. -*/
        options.add_int("CFOUR_FD_STEPSIZE", 5);

        /*- In finite difference calculations using the FINDIF option, this
        keyword specifies the point group to be used in generating the
        symmetry-adapted vibrational coordinates. FULL (= 0) specifies the
        full molecular point group, COMP (= 1) specifies the Abelian
        subgroup used in the electronic structure calculation. -*/
        options.add_str("CFOUR_FD_USEGROUP", "FULL", "FULL COMP");

        /*- This specifies the physical length (in integer words) of the
        records used in the word-addressable direct access files used by
        CFOUR. This value should always be chosen as a multiple of 512
        bytes, as your local system manager certainly understands. -*/
        options.add_int("CFOUR_FILE_RECSIZ", 2048);

        /*- This option allows the splitting of files. Input is required
        in the form N1/N2/N3/N4/N5, where N1, N2, N3, N4, and N5 specify
        the number of files in which ``MOINTS``, ``GAMLAM``, ``MOABCD``,
        ``DERINT``, and ``DERGAM`` are split, respectively. -*/
        options.add_str("CFOUR_FILE_STRIPE", "0/0/0/0/0");

        /*- Specifies the field strength for a perturbation (defined within
        a ``%perturbation`` section). The value must be given as an integer, and
        the field strength used by the program will be then the value of the
        keyword $\times 10^{-6}$. -*/
        options.add_int("CFOUR_FINITE_PERTURBATION", 0);

        /*- This option is used to control the algorithm used for
        construction of the Fock matrix in SCF calculations. PK (=0) uses
        the PK-supermatrix approach while AO (=1) constructs the matrix
        directly from the basis function integrals. In general, PK is
        somewhat faster, but results in considerable use of disk space when
        out-of-core algorithms are required. (Default: FOCK). -*/
        options.add_str("CFOUR_FOCK", "", "PK AO");

        /*- FREQ_ALGORIT experimental use -*/
        options.add_str("CFOUR_FREQ_ALGORITHM", "STANDARD", "STANDARD PARALLEL");

        /*- Specifies whether in the correlation treatment all electron (OFF
        =0) or only the valence electrons (ON =1) are considered. This
        keyword provides an alternative to the |cfour__cfour_dropmo| keyword, as it allows
        frozen-core calculation without explicitly specifying the
        corresponding inner-shell orbitals. -*/
        options.add_bool("CFOUR_FROZEN_CORE", false);

        /*- Specifies whether in the correlation treatment all virtual
        orbitals (OFF =0) or only a subset of virtual orbitals (ON =1) are
        used. In the latter case, the threshold for deleting virtual
        orbitals based on the orbital energy needs to be specified in a
        ``%frozen_virt`` section. -*/
        options.add_bool("CFOUR_FROZEN_VIRT", false);

        /*- Used to control the handling and storage of two-particle density
        matrix elements with four virtual indices $\Gamma(abcd)$. DISK (=0)
        directs the program to calculate and store all elements of
        $\Gamma(abcd)$, while DIRECT (=1) tells the program to use
        alternative algorithms in which $\Gamma(abcd)$ is calculated and
        used "on the fly". Note that this option might be not available
        for all type of calculations. -*/
        options.add_str("CFOUR_GAMMA_ABCD", "DISK", "DISK DIRECT");

        // GAMMA_ABCI
        // see GAMMA_ABCD

        /*- This keyword applies only to Hydrogen and Helium atoms and
        specifies the number of contracted Gaussian functions per shell.
        There is usually no need to use this keyword, but it can be useful
        for using a subset of the functions in a particular entry in the
        ``GENBAS`` file, particularly for generally contracted WMR basis sets.
        For example, if entry H:BASIS in the ``GENBAS`` file contains 7
        contracted s functions, 4 p functions and a single d function, then
        setting GENBAS_1=730 would eliminate the last p function and the d
        function. Default: use the unaltered ``GENBAS`` entry. -*/
        options.add_str("CFOUR_GENBAS_1", "");

        /*- This keyword performs the same function as |cfour__cfour_genbas_1|
        above, but applies to second-row atoms. -*/
        options.add_str("CFOUR_GENBAS_2", "");

        /*- This keyword performs the same function as |cfour__cfour_genbas_1| and
        |cfour__cfour_genbas_2| , but applies to third-row atoms. -*/
        options.add_str("CFOUR_GENBAS_3", "");

        /*- This keyword performs the same function as |cfour__cfour_genbas_1| ,
        |cfour__cfour_genbas_2| , and |cfour__cfour_genbas_3| , but applies
        to fourth-row atoms. -*/
        options.add_str("CFOUR_GENBAS_4", "");

        /*- Specifies the convergence criterion for geometry optimization.
        The optimization terminates when the RMS gradient is below $10^{-N}$
        Hartree/bohr, where $N$ is the specified value.
        **Psi4 Interface:** Geometry optimizations run through PSI (except in
        sandwich mode) use PSI's optimizer and so this keyword has no effect.
        Use :ref:`optking <apdx:optking>` keywords instead,
        particularly |optking__g_convergence| =CFOUR, which should be equivalent
        except for different internal coordinate definitions. -*/
        options.add_int("CFOUR_GEO_CONV", 5);

        /*- Specifies largest step (in millibohr) which is allowed in
        geometry optimizations.
        **Psi4 Interface:** Geometry optimizations run through PSI (except in
        sandwich mode) use PSI's optimizer and so this keyword has no effect.
        Use :ref:`optking <apdx:optking>` keywords instead,
        particularly |optking__intrafrag_step_limit|. -*/
        options.add_int("CFOUR_GEO_MAXSTEP", 300);

        /*- Specifies the used geometry optimization methods. The following
        values are permitted: NR (=0) --- straightforward Newton-Raphson
        search for minimum; RFA (=1) --- Rational Function Approximation
        search for minimum (this method can be used to find minima when the
        initial structure is in a region where the Hessian index is
        nonzero); TS (=2) Cerjan-Miller eigenvector following search for a
        transition state (can be started in a region where the Hessian index
        is not equal to unity); MANR (=3) --- Morse-adjusted Newton-Raphson
        search for minimum (very efficient minimization scheme, particularly
        if the Hessian is available); SINGLE_POINT (=5) for a single-point
        energy calculation. ENERONLY (=6) requests a geometry optimization
        based on single-point energy calculations.  Default: SINGLE-POINT
        (NR as soon as variables are marked to be optimized). -*/
        options.add_str("CFOUR_GEO_METHOD", "SINGLE_POINT", "NR RFA TS MANR SINGLE_POINT ENERONLY");

        /*- Specifies the maximum allowed number of geometry optimization cycles.
        **Psi4 Interface:** Geometry optimizations run through PSI (except in
        sandwich mode) use PSI's optimizer and so this keyword has no effect.
        Use :ref:`optking <apdx:optking>` keywords instead,
        particularly |optking__geom_maxiter|. -*/
        options.add_int("CFOUR_GEO_MAXCYC", 50);

        /*- Specifies whether gauge-including atomic orbitals are used (ON)
        or not (OFF). Default: ON for |cfour__cfour_props| =NMR  and =MAGNETIC,
        otherwise OFF -*/
        options.add_str("CFOUR_GIAO", "", "ON OFF");

        // GIMIC
        // experimental use

        /*- Keyword used to control type of grid calculation (see later
        section in this manual). Options are OFF (=0), no grid calculation;
        CARTESIAN (=1), steps are in Cartesian coordinates (which must be
        run with |cfour__cfour_coordinates| =CARTESIAN); INTERNAL (=2), steps are in Z-matrix
        internal coordinates; QUADRATURE (=3) steps are chosen for an
        integration based on Gauss-Hermite quadrature. (Default: OFF) -*/
        options.add_str("CFOUR_GRID", "OFF", "OFF CARTESIAN INTERNAL QUADRATURE");

        // GRID_ALGO
        // experimental use

        /*- Where the initial SCF eigenvectors are read from. MOREAD means
        to read from the disk (the ``JOBARC`` file) and CORE means to use a
        core Hamiltonian initial guess. If MOREAD is chosen but no disk file
        is present, the core Hamiltonian is used. (Default: MOREAD) -*/
        options.add_str("CFOUR_GUESS", "MOREAD", "MOREAD CORE");

        /*- This keyword determines which action is taken by the linear
        response program. ON (=1) the full effective Hamiltonian is
        calculated and written to disk; OFF (=0) the "lambda" linear
        response equations are solved. -*/
        options.add_bool("CFOUR_HBAR", false);

        // HESS_TYPE
        // experimental use

        /*- Control analysis of the stability of RHF, ROHF and UHF
        wavefunctions, as well as a possible search for a lower SCF
        solution. There are three possible options for this keyword. OFF
        (=0) does nothing, while ON (=1) performs a stability analysis and
        returns the number of negative eigenvalues in the orbital rotation
        Hessian. A third option, FOLLOW (=2) performs the stability
        analysis and then proceeds to rotate the SCF orbitals in the
        direction of a particular negative eigenvalue of the orbital
        rotation Hessian (see the explanation of keyword |cfour__cfour_rot_evec|), after
        which the SCF is rerun. -*/
        options.add_str("CFOUR_HFSTABILITY", "OFF", "OFF ON FOLLOW");

        // HF2_FILE
        // experimental use

        /*- This keyword can be used to significantly reduce disk i/o, and
        should be implemented very soon. The following options are
        available: OFF (= 0), no special algorithms are used (the default
        case); ALL (=1) all quantities except the $\langle ab\vert\vert
        cd\rangle$ molecular integral lists are held in core; PARTIAL (= 2),
        the T2 and T1 vectors are held in core throughout the calculation;
        (=4) all quantities except the $\langle ab\vert\vert cd\rangle$ and
        $\langle ab\vert\vert ci\rangle$ integrals are held in core; (=5)
        $\langle ij\vert\vert kl\rangle$ and $\langle ij\vert\vert
        ka\rangle$ and two-index quantities are held in core; (=6) all
        direct access files (``MOINTS``, ``GAMLAM``, etc.) are held in core. At
        present, these options have been implemented only in the energy code
        ``xvcc`` and the excitation energy code ``xvee``. (Default: 0) -*/
        options.add_str("CFOUR_INCORE", "OFF", "OFF ALL PARTIAL");

        /*- Specifies whether an input for mrcc is written (ON, =0) or not
        (OFF, =1) if |cfour__cfour_cc_program| =EXTERNAL has been specified. -*/
        options.add_bool("CFOUR_INPUT_MRCC", true);

        /*- This keyword defines what type of integral input will be written
        by ``xjoda``. VMOL (=1) has to be used with the programs of CFOUR. Using
        ARGOS (=0), input for Pitzer's ARGOS integral program will be
        written. (Default: VMOL). -*/
        options.add_str("CFOUR_INTEGRALS", "VMOL", "VMOL ARGOS");

        /*- Controls amount of debug printing performed by ``xjoda``. The higher
        the number, the more information is printed. Values of 25 or higher
        generally do not produce anything of interest to the general user.
        Do not set JODA_PRINT to 999 as this will cause the core vector to
        be dumped to disk. -*/
        options.add_int("CFOUR_JODA_PRINT", 0);

        /*- Convergence threshold for linear equations controlled by
        LINEQ_TYPE. Equations are iterated until smallest residual falls
        below $10^{-N}$, where $N$ is the value associated with this keyword. -*/
        options.add_int("CFOUR_LINEQ_CONV", 7);

        ///*- Maximum subspace dimension for linear equation solutions. -*/
        // options.add_int("CFOUR_LINEQ_EXPOR");

        /*- Determines the algorithm used to solve linear equations
        ( $\Lambda$ and derivative $T$ and $\Lambda$ ). POPLE (=0) uses
        Pople's method of successively orthogonalized basis vectors, while
        DIIS (=1) uses Pulay's DIIS method. The latter offers the practical
        advantage of requiring much less disk space, although it is not
        guaranteed to converge. Moreover, POPLE has not been tested for some
        time and should definitely be checked! (Default : DIIS) -*/
        options.add_str("CFOUR_LINEQ_TYPE", "DIIS", "POPLE DIIS");

        /*- The maximum number of iterations in all linear CC equations. -*/
        options.add_int("CFOUR_LINEQ_MAXCY", 50);

        // LINDEP_TOL

        /*- This keyword is used by the SCF program to determine if the
        orbital occupancy (by symmetry block) is allowed to change in the
        course of the calculation. ON (=1) locks the occupation to that set
        by the keyword |cfour__cfour_occupation| (or the initial guess if
        omitted); OFF (= 0) permits the occupation to change. (Default : 1
        if the occupation is specified with |cfour__cfour_occupation| and for
        second and later steps of optimizations; 0 if |cfour__cfour_occupation|
        omitted.) -*/
        options.add_bool("CFOUR_LOCK_ORBOCC", false);

        /*- Identical to |cfour__cfour_geo_maxstep|. -*/
        options.add_int("CFOUR_MAXSTEP", 300);

        /*- Specifies the amount of core memory used in integer words
        (default) or in the units specified via the keyword |cfour__cfour_mem_unit|.
        Default: 100 000 000 (approximately 381 or 762 MB for 32 or 64 bit
        machines, respectively).
        **Psi4 Interface:** Keyword set in MB from memory input command when
        given. -*/
        options.add_int("CFOUR_MEMORY_SIZE", 100000000);

        /*- Specifies the units in which the amount of requested core memory
        is given. Possible choices are INTEGERWORDS (default), kB, MB, GB,
        and TB.
        **Psi4 Interface:** Keyword set from memory input command when
        given, always MB. -*/
        options.add_str("CFOUR_MEM_UNIT", "INTEGERWORDS", "INTEGERWORDS KB MB GB TB");

        /*- Specifies the geometry optimization strategy. Four values are
        permitted: NR (=0) -- Straightforward Newton-Raphson search for
        minimum; RFA (=1) -- Rational Function Approximation search for
        minimum (this method can be used to find minima when the initial
        structure is in a region where the Hessian index is nonzero); TS (=2)
        Cerjan-Miller eigenvector following search for a transition
        state (can be started in a region where the Hessian index is not
        equal to unity); MANR (=3) -- Morse-adjusted Newton-Raphson search
        for minimum (very efficient minimization scheme, particularly if the
        Hessian is available); 4 is currently unavailable;
        SINGLE_POINT (=5) is a single point calculation.
        **Psi4 Interface:** Geometry optimizations run through PSI (except in
        sandwich mode) use PSI's optimizer and so this keyword has no effect.
        Use :ref:`optking <apdx:optking>` keywords instead,
        particularly |optking__opt_type| and |optking__step_type|. -*/
        options.add_str("CFOUR_METHOD", "SINGLE_POINT", "NR RFA TS MANR SINGLE_POINT");

        /*- Specifies the type of MRCC calculation. MK performs a MR-CC
        calculation based on Mukherjee's ansatz. -*/
        options.add_bool("CFOUR_MRCC", false);

        /*- Specifies the spin multiplicity.
        **Psi4 Interface:** Keyword set from active molecule. -*/
        options.add_int("CFOUR_MULTIPLICITY", 1);

        /*- Calculation of non-adiabatic coupling. In case of ON (=1) the
        method by Ichino, Gauss, Stanton is used to obtain the lambda
        coupling, while in case of LVC (=3) the lambda coupling is computed
        by means of the algorithm by Tajti and Szalay. Furthermore, NACV (=2)
        requests the computation of the full non-adiabatic coupling. Note
        that for calculations using LVC or NACV
        options the multiroot diagonalization has to be used, as requested
        via the keyword CFOUR_EOM_NSTATES (dne?) =MULTIROOT. -*/
        options.add_str("CFOUR_NACOUPLING", "OFF", "ON NACV LVC");

        /*- Specifies what to do if negative eigenvalues are encountered in
        the totally symmetric Hessian during an NR or MANR
        geometry-optimization search. If ABORT (=0), the job will
        terminate with an error message; if SWITCH (=1) the program
        will just switch the eigenvalue to its absolute value and keep
        plugging away (this is strongly discouraged!); and if RFA
        (=2), the keyword |cfour__cfour_geo_method| is switched to RFA internally and the
        optimization is continued.
        **Psi4 Interface:** Geometry optimizations run through PSI (except in
        sandwich mode) use PSI's optimizer and so this keyword has no effect.
        Use :ref:`optking <apdx:optking>` keywords instead. -*/
        options.add_str("CFOUR_NEGEVAL", "ABORT", "ABORT SWITCH RFA");

        /*- All components of spherical AO's are normalized to 1. This
        feature can help with numerical convergence issues if AO integrals
        are involved. Currently only working for single-point energy
        calculations. -*/
        options.add_bool("CFOUR_NEWNORM", false);

        /*- Specifies whether the reference function used in the correlation
        energy calculation satisfies the (spin-orbital) HF equations or not.
        Usually there is no need to set this parameter (OFF = 0 and ON =1),
        since standard non-HF reference functions (QRHF and ROHF) set this
        flag automatically. -*/
        options.add_bool("CFOUR_NONHF", false);

        /*- Specifies how many t amplitudes will be printed for each spin
        case and excitation level. For =N, The largest N amplitudes for each spin
        case and excitation level will be printed. -*/
        options.add_int("CFOUR_NTOP_TAMP", 15);

        /*- Specifies the orbital occupancy of the reference function in
        terms of the occupation numbers of the orbitals and their
        irreducible representations. The occupancy is specified by either
        NIRREP or 2*NIRREP integers specifying the number of occupied
        orbitals of each symmetry type, where NIRREP is the number of
        irreducible representations in the computational point group. If
        there are no orbitals of a particular symmetry type a zero must be
        entered. If the reference function is for an open-shell system, two
        strings of NIRREP occupation numbers separated by a slash are input
        for the $\alpha$ and $\beta$ sets of orbitals.  An example of the use of
        the OCCUPATION keyword for the water molecule would be
        OCCUPATION=3-1-1-0. For the :math:`^2A_1` water cation, an open-shell
        system, the keyword would be specified by
        OCCUPATION=3-1-1-0/2-1-1-0. It should be noted that the ``xvmol``
        integral program orders the irreducible representations in a strange
        way, which most users do not perceive to be a logical order. Hence,
        it is usually advisable initially to run just a single point
        integral and HF-SCF calculation in order to determine the number and
        ordering of the irreducible representations.  The occupation keyword
        may be omitted, in which case an initial orbital occupancy is
        determined by diagonalization of the core Hamiltonian. In many
        cases, HF-SCF calculations run with the core Hamiltonian guess will
        usually converge to the lowest energy HF-SCF solution, but this
        should not be blindly assumed.  (Default: The occupation is given
        by the core Hamiltonian initial guess).
        **Psi4 Interface:** The arrays above are specified in PSI as
        (white space tolerant) [3,1,1,0] and [[3,1,1,0],[3,0,1,0]]. -*/
        options.add("CFOUR_OCCUPATION", new ArrayType());

        /*- Specifies which kind of open-shell CC treatment is employed. The
        default is a spin-orbital CC treatment (SPIN-ORBITAL =1) which is
        the only possible choice for UHF-CC schemes anyways. For ROHF-CC
        treatments, the possible options are beside the standard
        spin-orbital scheme a spin-restricted CC approach (SR-CC=3), as well
        as a corresponding linear approximation (which in the literature
        usually is referred to as partially-spin-adapted CC scheme)
        (PSA-CC=1). SR-CC and PSA-CC are within the CCSD approximation
        restricted to excitations defined by the first-order interacting
        space arguments. With the keywords PSA-CC_FULL (=2) or SR-CC_FULL (=6)
        inclusion of the so called "pseudo-triples" beyond the first-order
        interacting space is also possible.  The two-determinant CC method
        for open-shell singlet states can be activated by TD-CC (=8). -*/
        options.add_str("CFOUR_OPEN-SHELL", "SPIN-ORBITAL", "SPIN-ORBITAL SR-CC PSA-CC_FULL SR-CC_FULL TD-CC");

        /*- Identical to |cfour__cfour_geo_maxcyc|. -*/
        options.add_int("CFOUR_OPT_MAXCYC", 50);

        /*- Specifies the type of molecular orbitals used in post-HF
        calculations. STANDARD (=0) requests usage of the orbitals (from a
        corresponding HF-SCF calculation) without any modification. These
        are in the case of RHF/UHF the usual canonical HF orbitals and in
        the case of ROHF calculations the standard ROHF-orbitals with equal
        spatial parts for both the $\alpha$ and the $\beta$ spin orbitals.
        SEMICANONICAL (=1) forces in ROHF type calculations a transformation
        to so-called semicanonical orbitals which diagonalize the
        occupied-occupied and virtual-virtual blocks of the usual
        Fock-matrices. The use of semicanonical orbitals is, for example,
        required for ROHF-CCSD(T) calculations and for those calculations
        also automatically set. LOCAL requests a localization of the HF
        orbitals and this is currently done according to the Pipek-Mezey
        localization criterion.  Note that it is strongly recommended not to
        use this keyword unless you know what are you doing.  Default:
        STANDARD except for ROHF-CCSD(T) and ROHF-MP4 calculations for which
        SEMICANONICAL is the default. -*/
        options.add_str("CFOUR_ORBITALS", "STANDARD", "STANDARD SEMICANONICAL");

        // PARALLEL
        // experimental use

        // PARA_PRINT
        // experimental use

        // PARA_INT
        // experimental use

        /*- Specifies the type of perturbed orbitals used in energy
        derivative calculations. STANDARD means that the gradient
        formulation assumes that the perturbed orbitals are not those in
        which the (perturbed) Fock matrix is diagonal. CANONICAL means that
        the perturbed orbitals are assumed to be canonical. This keyword is
        set automatically to CANONICAL in derivative calculations with
        methods which include triple excitations (MBPT[4]/MP4, CCSD+T[CCSD],
        CCSD[T], QCISD[T] and all iterative schemes like CCSDT-n and CC3)
        apart from CCSDT. IJ_CANONICAL requests a canonical
        perturbed-orbital treatment only for the occupied-occupied block of
        the unperturbed density matrix in analytic derivative calculations.
        For testing purposes, it is possible to force the use standard
        perturbed orbitals even in case of iterative triple excitations via
        the option FORCE_STANDA (dne?).  Note also that in case of unrelaxed
        derivatives standard orbitals must be used.  Default : STANDARD for
        all methods without triples (except CCSDT), CANONICAL for all
        methods with triples in case of relaxed derivatives. -*/
        options.add_str("CFOUR_PERT_ORB", "", "STANDARD CANONICAL IJ_CANONICAL");

        /*- Specifies either single (=1, or SINGLE) or double (=2, DOUBLE)
        sided numerical differentiation in the finite difference evaluation
        of the Hessian. Two-sided numerical differentiation is considerably
        more accurate than the single-sided method, and its use is strongly
        recommended for production work. -*/
        options.add_str("CFOUR_POINTS", "DOUBLE", "SINGLE DOUBLE");

        /*- Controls the amount of printing in the energy and energy
        derivative calculation programs. Using a value of 1 will produce a
        modest amount of additional output over the default value of 0,
        which includes some useful information such as SCF eigenvectors,
        Fock matrix elements, etc. -*/
        options.add_int("CFOUR_PRINT", 0);

        /*- Specifies whether and which molecular property is calculated.
        OFF (=0) means that no property is calculated, FIRST_ORDER (=1)
        requests computation of various one-electron first-order properties
        (e.g., dipole moment, quadrupole moment, electric field gradient,
        spin densities,etc.), SECOND_ORDER (=2, in the next release replaced
        by STAT_POL) computes static electric polarizabilities, DYNAMICAL
        (=7, in the next release replaced by DYN_POL) requests the
        calculation of frequency-dependent polarizabilities (note that here
        an additional input of the frequency is required), NMR (=5) requests
        the calculation of NMR chemical shifts/chemical shielding tensors
        (by default using GIAOs), J_FC requests the calculation of the
        Fermi-Contact contribution to indirect spin-spin coupling constants,
        J_SD the calculation of the corresponding spin-dipole contribution,
        and J_SO the calculation of the corresponding spin-orbit
        contribution to J; HYPERPOL (=22) invokes a calculation of static
        hyperpolarizabilities, DYN_HYP (=23) requests the calculation of
        frequency-dependent hyperpolarizabilities, SHG (=24) the calculation
        of hyperpolarizabilities related to the second-harmonic
        generation, OPT_REC (=25) the computation of hyperpolarizabilities
        related to optical rectification, VERDET (=26) the calculation of
        Verdet constants. -*/
        options.add_str("CFOUR_PROPS", "OFF", "OFF FIRST_ORDER SECOND_ORDER NMR HYPERPOL DYN_HYP SHG OPT_REC VERDET");

        /*- Allows storage of property integrals computed in ``xvdint`` on
        internal files (e.g., ``MOINTS`` and ``GAMLAM``, default choice INTERNAL,
        =0) or on external files (EXTERNAL, =1). -*/
        options.add_str("CFOUR_PROP_INTEGRAL", "INTERNAL", "INTERNAL EXTERNAL");

        /*- The presence of this keyword specifies that a QRHF based CC
        calculation, or alternatively, an SCF calculation that uses the
        |cfour__cfour_qrhfgues| option, is to be performed. -*/
        options.add("CFOUR_QRHF_GENERAL", new ArrayType());

        /*- If this keyword is set to ON (=1), then the QRHF orbitals
        specified by the |cfour__cfour_qrhf_general|, |cfour__cfour_qrhf_orbital|
        and CFOUR_QRHF_SPIN (nyi?) keywords
        are used as a starting guess for a restarted SCF procedure. This can
        be an extremely useful way to converge "difficult" SCF solutions,
        such as those that correspond to states that are not the lowest
        states of a given symmetry. Note that when this option is used, the
        calculation that is performed is not a QRHF-CC calculation; it is
        instead a UHF-based or ROHF-based calculation, depending on what
        type of reference is specified by the |cfour__cfour_reference| keyword. The QRHF
        aspect of the calculation is used simply as a device to converge the
        orbitals. -*/
        options.add_bool("CFOUR_QRHFGUES", false);

        /*- By default, in QRHF calculations, electrons are removed from the
        highest occupied orbital in a symmetry block (symmetry block HOMO),
        while electrons are added to the lowest unoccupied orbital within a
        symmetry block (symmetry block LUMO). The purpose of the
        QRHF_ORBITAL keyword is to allow additional flexibility in choosing
        which orbitals will have their occupation numbers altered. The value
        of this keyword gives the offset with respect to the default orbital
        for the orbital which will be depopulated (or populated) in QRHF-CC
        calculations. For calculations involving more than one removal or
        addition of electrons, values are separated by commas and correspond
        to the |cfour__cfour_qrhf_general| input on a one-to-one basis. For example,
        specifying |cfour__cfour_qrhf_general| =2/-4, QRHF_ORBITAL=3/2 means that an electron
        will be added to the third lowest virtual in symmetry block 2 and
        another will be removed from the second highest occupied orbital in
        symmetry block 4. Examples given later in this manual further
        illustrate the QRHF input options and may help to clarify any
        confusion resulting from this documentation. (Default : 1) -*/
        options.add("CFOUR_QRHF_ORBITAL", new ArrayType());

        ///*- Specifies the spin of the electrons modified by the QRHF_GENERAL
        // and QRHF_ORBITAL keywords, where a value of 1 means $\alpha$ spin,
        // while 2 corresponds to a $\beta$ electron. By default, electrons that
        // are removed are assigned to $\beta$ spin, while added electrons are
        //$\alpha$. Note that this option allows one to construct low-spin
        // determinants, which generally are unsuitable for single-reference
        // coupled-cluster calculations. An important exception is the
        // open-shell singlet coupled-cluster method (see keyword
        // OPEN-SHELL=TD-CC above). -*/
        // options.add_int("CFOUR_QRHF_SPIN");

        /*- ON (=1) requests a calculation of Raman intensities based on the
        geometrical derivatives of the static polarizability tensor, while
        DYN (=2) requests a calculation of Raman intensities based on the
        derivatives of the dynamical polarizability tensor. -*/
        options.add_str("CFOUR_RAMAN_INT", "OFF", "ON DYN OFF");

        /*- Specifies whether Raman intensities are calculated with orbital
        relaxation with respect to the electric field perturbation (RELAXED,
        = 1) or without orbital relaxation (UNRELAXED, = 0). -*/
        options.add_str("CFOUR_RAMAN_ORB", "UNRELAXED", "RELAXED UNRELAXED");

        /*- Specifies whether or not relaxed density natural orbitals are to
        be computed. This option only has meaning for a correlated
        calculation. For =0, Do not compute. For =1, compute. -*/
        options.add_bool("CFOUR_RDO", true);

        /*- Specifies the type of SCF calculation to be performed. RHF (= 0)
        requests a restricted Hartree-Fock reference; UHF (= 1) an
        unrestricted Hartree-Fock reference; ROHF (= 2) a restricted
        open-shell Hartree-Fock calculation; TCSCF (=3) a
        two-configurational SCF calculation, and CASSCF (=4) a
        complete-active space SCF calculations (currently not implemented).
        **Psi4 Interface:** Keyword subject to translation from value of
        |scf__reference| unless set explicitly. -*/
        options.add_str("CFOUR_REFERENCE", "RHF", "RHF UHF ROHF TCSCF CASSCF");

        /*- Specifies the treatment of relativistic effects. The default is
        a non-relativistic treatment (OFF), while perturbational treatments
        are invoked via MVD1 (mass-velocity and 1-electron Darwin
        contribution), MVD2 (mass-velocity and 1- and 2-electron Darwin
        contribution), DPT2 (second-order direct perturbation theory
        approach), SF-DPT4 (scalar-relativistic part of fourth-order direct
        perturbation theory, DPT4 (full fourth-order DPT including
        spin-orbit corrections), SF-DPT6 (scalar-relativistic part of
        sixth-order direct perturbation theory), SFREE (spin-free
        treatment), X2C1E (spin-free X2C-1e treatment), or DPT (synonym with
        DPT2). -*/
        options.add_str("CFOUR_RELATIVISTIC", "OFF", "OFF MVD1 MVd2 DPT2 SF-DPT4 DPT4 SF-DPT6 SFREE X2C1E DPT");

        /*- Specifies whether the relaxed density matrix is computed for
        correlated wave functions. OFF (= 0) The relaxed density will not be
        computed, ON (= 1) it will be computed. -*/
        options.add_bool("CFOUR_RELAX_DENS", false);

        /*- This option can be used to convert an analytically calculated
        gradient vector to a particular normal coordinate representation. A
        useful application is to calculate the gradient of an electronically
        excited state in the normal coordinate representation of the ground
        electronic state, as this provides a first approximation to
        resonance Raman intensities (hence the name of the keyword).
        Calculations that use the this option require the externally
        supplied force constant matrix ``FCMFINAL``, which is written to disk
        during the course of both analytic and finite-difference vibrational
        frequency calculations. No such transformation is performed if OFF
        (=0); while ON (=1) directs the program to evaluate the gradient and
        transform it to the chosen set of normal coordinates. A warning
        message is printed if the force constant matrix is unavailable. -*/
        options.add_bool("CFOUR_RES_RAMAN", false);

        // RESET_FLAGS
        // experimental use

        /*- Offers the possibility to restart a CC calculation which stopped
        for various reasons, e.g. time limit, in the correlation part.
        However, note that a restart which is specified by ON (= 1) needs
        the following files of the previous unfinished calculation: ``JOBARC``,
        ``JAINDX``, ``MOINTS``, and ``MOABCD``. -*/
        options.add_bool("CFOUR_RESTART_CC", false);

        /*- Specifies which eigenvector of the orbital rotation Hessian is
        to be used to rotate the original SCF orbitals. By default, it will
        use that associated with the lowest eigenvalue of the totally
        symmetric part of the block-factored Hessian, as this choice often
        leads to the lowest energy SCF solution. For RHF stability checks,
        only those instabilities which correspond to RHF solutions will be
        considered. It is important to understand that following
        non-symmetric eigenvectors lowers the symmetry of the wavefunction
        and that following RHF --> UHF stabilities leads to a UHF solution.
        To converge the SCF roots associated with such instabilities, one
        must run the calculation in reduced symmetry and as a closed-shell
        UHF case, respectively. Value *n* directs the program to follow the
        vector associated with the *n*\ th lowest eigenvalue having the proper
        symmetry (totally symmetric) and spin (RHF-->RHF or UHF-->UHF)
        properties. 0 means use the lowest eigenvalue. -*/
        options.add_int("CFOUR_ROT_EVEC", 0);

        /*- Tells CFOUR whether to delete large files (AO integrals and
        ``MOINTS`` file for now) when they are no longer needed. OFF (=0) They
        will not be saved, ON (=1) they will be saved. -*/
        options.add_bool("CFOUR_SAVE_INTS", false);

        /*- Controls whether step scaling is based on the absolute step
        length (1-norm) (=0 or MAG(S)) or the largest individual step in the
        internal coordinate space (=1 or MAX(S)). -*/
        options.add_str("CFOUR_SCALE_ON", "MAG(S)", "MAG(S) MAX(S)");

        /*- Specifies the convergence criterion for the HF-SCF equations.
        Equations are considered converged when the maximum change in
        density matrix elements is less than $10^{-N}$.
        **Psi4 Interface:** Keyword subject to translation from value of
        |scf__d_convergence| unless set explicitly. -*/
        options.add_int("CFOUR_SCF_CONV", 7);

        /*- Controls the damping (in the first iterations (specified by
        |cfour__cfour_scf_expstart| via
        :math:`D_{new} = D_{old} + X/1000 * (D_{new} - D_{old})` with $X$
        as the value specified by the keyword. The default value is
        currently 1000 (no damping), but a value of 500 is recommended in
        particular for transition metal compounds where the SCF convergence
        is often troublesome.
        **Psi4 Interface:** Keyword subject to translation from value of
        |scf__damping_percentage| unless set explicitly. -*/
        options.add_int("CFOUR_SCF_DAMPING", 1000);

        /*- Specifies the number of density matrices to be used in the
        DIIS convergence acceleration procedure. -*/
        options.add_int("CFOUR_SCF_EXPORDER", 6);

        /*- Specifies the first iteration in which the DIIS convergence
        acceleration procedure is applied. -*/
        options.add_int("CFOUR_SCF_EXPSTART", 8);

        /*- Specifies whether or not the DIIS extrapolation is used to
        accelerate convergence of the SCF procedure. OFF (=0) means do not use
        DIIS, ON (=1) means use DIIS. -*/
        options.add_bool("CFOUR_SCF_EXTRAPOLATION", true);

        /*- Specifies the maximum number of SCF iterations.
        **Psi4 Interface:** Keyword subject to translation from value of
        |scf__maxiter| unless set explicitly.-*/
        options.add_int("CFOUR_SCF_MAXCYC", 150);

        /*- Specifies the strength of a spin-dipole perturbation as required
        for finite-field calculations of the SD contributions to indirect
        spin-spin coupling constants. The value must be specified as an
        integer and the SD strength used by the program will be the value of
        the keyword $\times 10^{-6}$. (Default : 0, currently not implemented) -*/
        options.add_int("CFOUR_SD_FIELD", 0);

        // SOPERT
        // Experimental Use!
        // Default : OFF.
        ///*- Perturbative treatment of spin-orbit splittings in doublet-pi
        // states via multireference coupled-cluster theory. MKMRCC (=1)
        // requests a treatment based on Mukherjee's multireference
        // coupled-cluster theory. EMRCCSO (=2) requests the expectation value
        // of a similarity transformed spin-orbit operator. Please note that
        // symmetric orbitals are needed, e.g., using AV_SCF. For more
        // information on the theory see J. Chem. Phys. 136, 111103 (2012). -*/

        /*- Specifies whether spherical harmonic (5d, 7f, 9g, etc.) or
        Cartesian (6d, 10f, 15g, etc.) basis functions are to be used. ON (=
        1) uses spherical harmonics, OFF (= 0) uses Cartesians.
        **Psi4 Interface:** Keyword set according to basis design when
        |mints__basis| is used instead of |cfour__cfour_basis|. Keyword
        subject to translation from value of |globals__puream| unless set
        explicitly. -*/
        options.add_bool("CFOUR_SPHERICAL", true);

        /*- Controls whether excitation energy calculations allow for a
        "spin flip" which changes the $M_s$ quantum number. Such
        calculations have some advantages for biradicals and are currently
        implemented (together with gradients) for CIS and CIS(D)
        calculations. Options are OFF and ON. -*/
        options.add_bool("CFOUR_SPIN_FLIP", false);

        /*- Experimental Use!  ON (=1) requests calculation of one-electron
        spin-orbit integrals. MEANSO additionally gives a mean-field
        treatment of the two-electron terms (spin-orbit mean field treatment
        as described Mol. Phys. 98, 1823-1833 (2000)). -*/
        options.add_str("CFOUR_SPIN_ORBIT", "OFF", "ON MEANSO OFF");

        // SPINORBIT
        // experimental use

        /*- ON (=1) requests the spin-component scaled variant of the MP2
        approach. This keyword has only an effect when |cfour__cfour_calc_level| =MP2 is
        specified and must be used together with |cfour__cfour_reference| =UHF. -*/
        options.add_bool("CFOUR_SPIN_SCAL", false);

        /*- Specifies whether nuclear spin-rotation tensors are computed
        within a NMR chemical shift calculation (ON, =1) or not (OFF, =9).
        In the case of electronic g-tensor calculations for open-shell
        molecules this keyword controls the calculation of the electronic
        spin-rotation tensor. -*/
        options.add_bool("CFOUR_SPINROTATION", false);

        /*- Specifies an Abelian subgroup to be used in a calculation.
        Acceptable arguments are DEFAULT (=0); C1 (= 1); C2 (= 2); CS (= 3);
        CI (= 4); C2V (= 5); C2H (= 6); D2 (= 7) and D2H (= 8). Use of C1 is
        of course equivalent to setting |cfour__cfour_symmetry| =OFF in the input. The
        DEFAULT option (which is the default) uses the highest order Abelian
        subgroup. -*/
        options.add_str("CFOUR_SUBGROUP", "DEFAULT", "DEFAULT C1 C2 CS CI C2V C2H D2 D2H OFF");

        ///*- Is a somewhat complicated keyword to use. Allowed values are 0,
        // 1, and 2, which specify the $x$, $y$, and $z$ axes, respectively. The
        // meaning of the keyword is best described by example: Suppose one is
        // running a calculation on water, and wishes to run it in the $C_s$
        // point group with the "special" plane being the one which bisects the
        // H-O-H bond angle. Now, what SUBGRPAXIS does is to specify which
        // Cartesian direction in the $C_{2v}$ frame becomes the special
        // direction in the $C_s$ frame. CFOUR will orient water in the $yz$
        // plane, so one wants the $y$ axis in the $C_{2v}$ frame to be the $z$
        // axis in the $C_s$ frame. Hence, for this case, one would specify
        // SUBGRPAXIS=2. Use of this keyword may be facilitated by studying
        // section D1 of this chapter, entitled "Molecular Orientation".
        // However, when the true Abelian subgroup is either $C_{2v}$ or
        //$D_{2h}$, the CFOUR orientation is not well defined, and it may be
        // necessary to run the ``xjoda`` executable directly two times. If
        // SUBGROUP=0 in the first pass, then the reference orientation for the
        // true Abelian subgroup can be determined and the appropriate value of
        // SUBGRPAXIS selected. -*/
        // options.add_int("CFOUR_SUBGRPAXIS");

        /*- In principle can be used to force the SCF to converge a solution
        for which the density matrix transforms as the totally symmetric
        representation of the point group (i.e. no broken symmetry
        solutions). The code seems to work in most cases, but has currently
        been implemented for point groups with E type representation and not
        for those with triply-, quadruply- or pentuply-degenerate
        representations. Extending the code to those cases is probably
        straightforward, and the reader is encouraged to do so if (s)he is
        so inclined. SYM_CHECK=0 "forces" the high-symmetry solution.
        SYM_CHECK=OVERRIDE (=1) doesn't. The latter is the default. -*/
        options.add_bool("CFOUR_SYM_CHECK", true);

        /*- Specifies what subgroup of the full point group is to be used in
        the energy and/or gradient calculation (the computational point
        group). OFF (=1) forces a no symmetry run (in $C_1$ ) and ON (=0) runs
        the calculation in the largest self-adjoint subgroup ( $D_{2h}$ and its
        subgroups). -*/
        options.add_bool("CFOUR_SYMMETRY", true);

        /*- Specifies how often the largest $t$ amplitudes are to be printed.
        For =0, amplitudes are printed at the beginning and end of the run.
        For =1, amplitudes are printed every iteration. For =2, amplitudes are
        printed every other iteration, etc. -*/
        options.add_int("CFOUR_TAMP_SUM", 5);

        // TDHF
        // experimental use

        // TESTSUITE
        //(currently not available)

        /*- Specifies whether to calculate finite-temperature thermodynamic
        corrections after a frequency calculation. OFF (=0) skips this; ON
        (=1) gives abbreviated output; and VERBOSE (=2) gives elaborate
        output that is separated by translation, rotation and vibration.
        Default: ON (currently not available in public version) -*/
        options.add_str("CFOUR_THERMOCHEMISTRY", "ON", "OFF ON VERBOSE");

        // TRANGRAD
        // experimental use

        /*- Specifies whether or not translational invariance is exploited
        in geometrical derivative calculations. USE(=0) specifies that
        translational invariance is exploited, while IGNORE (=1) turns it
        off. -*/
        options.add_str("CFOUR_TRANS_INV", "USE", "USE IGNORE");

        /*- Specifies whether in a correlated NMR chemical shift
        calculations all perturbations are treated at once or sequentially.
        Available option are SIMULTANEOUS (=0) and SEQUENTIAL (=1). The
        latter is at least preferred for large-scale calculations, as it has
        less demands on the available disk space. -*/
        options.add_str("CFOUR_TREAT_PERT", "SIMULTANEOUS", "SIMULTANEOUS SEQUENTIAL");

        /*- Specifies whether the T3 amplitudes are included ON (=1) or not
        included OFF (=0) in the DIIS convergence acceleration during CCSDT
        calculations. Inclusion of T3 speeds up convergence and allows tight
        convergence, but on the other hand it increases disk space
        requirements. Note that this keyword is only available with module
        ``xecc``. -*/
        options.add_bool("CFOUR_T3_EXTRAPOL", false);

        /*- Specifies the threshold value (given as an integer) for the
        treatment of CPHF coefficients in second derivative calculations
        using perturbed canonical orbitals. If a CPHF coefficient is above
        the threshold, the corresponding orbital rotation is treated (at the
        expense of additional CPU cost) using the standard non-canonical
        procedures, while orbital pairs corresponding to CPHF coefficients
        below the threshold are treated using perturbed canonical
        representation.  Default: 25 (Default: 1 in the developer version) -*/
        options.add_int("CFOUR_UIJ_THRESHOLD", 25);

        /*- Specifies the units used for molecular geometry input. ANGSTROM
        (= 0) uses Angstrom units, BOHR (= 1) specifies atomic units.
        **Psi4 Interface:** Keyword set from active molecule, always ANGSTROM. -*/
        options.add_str("CFOUR_UNITS", "ANGSTROM", "ANGSTROM BOHR");

        // UNOS
        // experimental use

        /*- Specifies whether or not the Hessian update is carried out. OFF
        (= 0) uses the initial Hessian (however supplied, either the default
        guess or a ``FCMINT`` file), ON (= 1) updates it during subsequent
        optimization cycles. (not in current public version). -*/
        options.add_bool("CFOUR_UPDATE_HESSIAN", true);

        // VIB_ALGORIT
        // experimental use

        // VIBPHASE
        // experimental use

        /*- Specifies whether (harmonic) vibrational frequencies are
        calculated or not. If the default NO (=0) is specified then no
        frequencies are calculated. For ANALYTIC, vibrational frequencies
        are determined from analytically computed second derivatives, and
        for FINDIF (=2) vibrational frequencies are calculated from a force
        field obtained by numerical differentiation of analytically
        evaluated gradients (or even single-point energies) using
        symmetry-adapted mass-weighted Cartesian coordinates. If vibrational
        frequencies are calculated, a normal mode analysis using the
        computed force-constant matrix is performed, rotationally projected
        frequencies are computed, infrared intensities are determined, and
        zero-point energies (ZPE) are evaluated. -*/
        options.add_str("CFOUR_VIBRATION", "NO", "NO ANALYTIC FINDIF EXACT");

        /*- This keyword defines what type of integral transformation is to
        be performed in the program ``xvtran``. FULL/PARTIAL (=0) allows the
        transformation program to choose the appropriate type of
        transformation, while FULL (=1) requires a full integral
        transformation and PARTIAL (=2) means a MBPT(2)-specific
        transformation where the :math:`(ab \vert cd)` integrals are not formed. -*/
        options.add_str("CFOUR_VTRAN", "FULL/PARTIAL", "FULL/PARTIAL FULL PARTIAL");

        /*- Specifies the X-component of an external electric field. The
        value must be specified as an integer and the field used by the
        program will be the value of the keyword :math:`\times 10^{-6}`. This allows
        field strengths :math:`|\varepsilon| > 10^{-6}` to be used. -*/
        options.add_int("CFOUR_XFIELD", 0);

        /*- The tolerance for storing transformed integrals. Integrals less
        than $10^{-N}$ are neglected and not stored on disk. -*/
        options.add_int("CFOUR_XFORM_TOL", 11);

        /*- Specifies the Y-component of an external electric field. The
        value must be specified as an integer and the field used by the
        program will be the value of the keyword :math:`\times 10^{-6}`. This allows
        field strengths :math:`|\varepsilon| > 10^{-6}` to be used. -*/
        options.add_int("CFOUR_YFIELD", 0);

        /*- Specifies the Z-component of an external electric field. The
        value must be specified as an integer and the field used by the
        program will be the value of the keyword :math:`\times 10^{-6}`. This allows
        field strengths :math:`|\varepsilon| > 10^{-6}` to be used. -*/
        options.add_int("CFOUR_ZFIELD", 0);
    }
    if (name == "EFP" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs effective fragment potential
        computations through calls to Kaliman's libefp library. -*/

        /*- The amount of information printed to the output file. -*/
        options.add_int("PRINT", 1);
        /*- Do include electrostatics energy term in EFP computation? -*/
        options.add_bool("EFP_ELST", true);
        /*- Do include exchange repulsion energy term in EFP computation? -*/
        options.add_bool("EFP_EXCH", true);
        /*- Do include polarization energy term in EFP computation? (EFP_POL c. v1.1) -*/
        options.add_bool("EFP_IND", true);
        /*- Do include dispersion energy term in EFP computation? -*/
        options.add_bool("EFP_DISP", true);
        /*- Fragment-fragment electrostatic damping type. ``SCREEN``
        is a damping formula based on screen group in the EFP potential.
        ``OVERLAP`` is damping that computes charge penetration energy. -*/
        options.add_str("EFP_ELST_DAMPING", "SCREEN", "SCREEN OVERLAP OFF");
        /*- Fragment-fragment polarization damping type. ``TT`` is a
        damping formula like Tang and Toennies. (EFP_POL_DAMPING c. v1.1) -*/
        options.add_str("EFP_IND_DAMPING", "TT", "TT OFF");
        /*- Fragment-fragment dispersion damping type. ``TT`` is a damping
        formula by Tang and Toennies. ``OVERLAP`` is overlap-based
        dispersion damping. -*/
        options.add_str("EFP_DISP_DAMPING", "OVERLAP", "TT OVERLAP OFF");
        /*- Do include electrostatics energy term in QM/EFP computation? (QMEFP_ELST c. v1.1) -*/
        options.add_bool("EFP_QM_ELST", true);
        /*- Do include polarization energy term in QM/EFP computation? (QMEFP_POL c. v1.1) -*/
        options.add_bool("EFP_QM_IND", true);
        /*- Do EFP gradient? !expert -*/
        options.add_str("DERTYPE", "NONE", "NONE FIRST");
        /*- Do turn on QM/EFP terms? !expert -*/
        options.add_bool("QMEFP", false);
    }
    if (name == "DMRG" || options.read_globals()) {
        /*- MODULEDESCRIPTION Performs a DMRG computation
         through calls to Wouters's CheMPS2 library. -*/

        /*- The DMRG wavefunction multiplicity in the form (2S+1) -*/
        options.add_int("DMRG_MULTIPLICITY", -1);

        /*- The DMRG wavefunction irrep uses the same conventions as PSI4. How convenient :-).
            Just to avoid confusion, it's copied here. It can also be found on
            http://sebwouters.github.io/CheMPS2/doxygen/classCheMPS2_1_1Irreps.html .
            Symmetry Conventions        Irrep Number & Name
            Group Number & Name         0     1     2     3     4     5     6     7
            0: c1                       A
            1: ci                       Ag     Au
            2: c2                       A     B
            3: cs                       A'     A''
            4: d2                       A     B1     B2     B3
            5: c2v                      A1     A2     B1     B2
            6: c2h                      Ag     Bg     Au     Bu
            7: d2h                      Ag     B1g     B2g     B3g     Au     B1u     B2u     B3u
        -*/
        options.add_int("DMRG_IRREP", -1);

        /*- The number of reduced renormalized basis states to be
            retained during successive DMRG instructions -*/
        options.add("DMRG_SWEEP_STATES", new ArrayType());

        /*- The energy convergence to stop an instruction
            during successive DMRG instructions -*/
        options.add("DMRG_SWEEP_ENERGY_CONV", new ArrayType());

        /*- The density RMS convergence to stop an instruction
            during successive DMRG instructions -*/
        options.add_double("DMRG_SCF_GRAD_THR", 1.e-6);

        /*- The maximum number of sweeps to stop an instruction
            during successive DMRG instructions -*/
        options.add("DMRG_SWEEP_MAX_SWEEPS", new ArrayType());

        /*- The noise prefactors for successive DMRG instructions -*/
        options.add("DMRG_SWEEP_NOISE_PREFAC", new ArrayType());

        /*- The residual tolerances for the Davidson diagonalization during DMRG instructions -*/
        options.add("DMRG_SWEEP_DVDSON_RTOL", new ArrayType());

        /*- Whether or not to print the correlation functions after the DMRG calculation -*/
        options.add_bool("DMRG_PRINT_CORR", false);

        /*- Whether or not to create intermediary MPS checkpoints -*/
        options.add_bool("DMRG_MPS_WRITE", false);

        /*- Whether or not to store the unitary on disk (convenient for restarting). -*/
        options.add_bool("DMRG_UNITARY_WRITE", true);

        /*- Whether or not to use DIIS for DMRG. -*/
        options.add_bool("DMRG_DIIS", false);

        /*- When the update norm is smaller than this value DIIS starts. -*/
        options.add_double("DMRG_SCF_DIIS_THR", 1e-2);

        /*- Whether or not to store the DIIS checkpoint on disk (convenient for restarting). -*/
        options.add_bool("DMRG_DIIS_WRITE", true);

        /*- Maximum number of DMRG iterations -*/
        options.add_int("DMRG_SCF_MAX_ITER", 100);

        /*- Which root is targeted: 0 means ground state, 1 first excited state, etc. -*/
        options.add_int("DMRG_EXCITATION", 0);

        /*- Whether or not to use state-averaging for roots >=2 with DMRG-SCF. -*/
        options.add_bool("DMRG_SCF_STATE_AVG", true);

        /*- Which active space to use for DMRG calculations:
               --> input with SCF rotations (INPUT);
               --> natural orbitals (NO);
               --> localized and ordered orbitals (LOC) -*/
        options.add_str("DMRG_SCF_ACTIVE_SPACE", "INPUT", "INPUT NO LOC");

        /*- Whether to start the active space localization process from a random unitary matrix instead of a unit matrix. -*/
        options.add_bool("DMRG_LOCAL_INIT", true);

        /*- Do calculate the DMRG-CASPT2 energy after the DMRGSCF calculations are done? -*/
        options.add_bool("DMRG_CASPT2_CALC", false);

        /*- Whether to calculate the DMRG-CASPT2 energy after the DMRGSCF calculations are done. -*/
        options.add_str("DMRG_CASPT2_ORBS", "PSEUDOCANONICAL", "PSEUDOCANONICAL ACTIVE");

        /*- CASPT2 IPEA shift -*/
        options.add_double("DMRG_CASPT2_IPEA", 0.0);

        /*- CASPT2 Imaginary shift -*/
        options.add_double("DMRG_CASPT2_IMAG", 0.0);

        /*- DMRG-CI or converged DMRG-SCF orbitals in molden format -*/
        options.add_bool("DMRG_MOLDEN_WRITE", false);

        /*- Print out the density matrix in the AO basis -*/
        options.add_bool("DMRG_OPDM_AO_PRINT", false);
    }

    return true;
}

}  // namespace psi

// clang-format on
//  LocalWords:  Psi4
