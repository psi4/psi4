.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2023 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This file is part of Psi4.
.. #
.. # Psi4 is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU Lesser General Public License as published by
.. # the Free Software Foundation, version 3.
.. #
.. # Psi4 is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU Lesser General Public License for more details.
.. #
.. # You should have received a copy of the GNU Lesser General Public License along
.. # with Psi4; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #

.. include:: autodoc_abbr_options_c.rst

.. index:: SCF, HF, Hartree--Fock
.. _`sec:scf`:

HF: Hartree--Fock Theory
========================

.. codeauthor:: Justin M. Turney, Robert M. Parrish, and Andrew C. Simmonett
.. sectionauthor:: Robert M. Parrish and Jerome F. Gonthier

*Module:* :ref:`Keywords <apdx:scf>`, :ref:`PSI Variables <apdx:scf_psivar>`, :source:`LIBSCF_SOLVER <psi4/src/psi4/libscf_solver>`, :source:`LIBMINTS <psi4/src/psi4/libmints>`, :source:`LIBFOCK <psi4/src/psi4/libfock>`, :source:`LIBDIIS <psi4/src/psi4/libdiis>`

.. _`sec:scfintro`:

Introduction
~~~~~~~~~~~~

Self-Consistent-Field (SCF) theory forms the cornerstone of *ab initio* quantum
chemistry. Here SCF refers both to conventional Hartree--Fock (HF) molecular
orbital theory and also to generalized Kohn--Sham Density Functional Theory
(KS-DFT). |PSIfour| contains a wholly rewritten SCF code, including many of the
most popular spin specializations, several efficient numerical methods for
treating Fock Matrix construction, and a brand new KS-DFT code featuring many of
the most popular DFT functional technologies.

An illustrative example of using the SCF module is as follows::

    molecule {
    0 3
    O
    O 1 1.21
    }

    set {
    basis cc-pvdz
    guess sad
    reference uhf
    scf_type direct
    }

    energy('scf')

This will run a UHF computation for triplet molecular oxygen (the ground state)
using a Direct algorithm for the Electron Repulsion Integrals (ERI) and starting
from a Superposition of Atomic Densities (SAD) guess. DF integrals are
automatically used to converge the DF-SCF solution before the Direct algorithm is
activated.  After printing all manner of titles, geometries, sizings, and
algorithm choices, the SCF finally reaches the iterations::

                           Total Energy        Delta E     RMS |[F,P]|

   @DF-UHF iter   0:  -149.80032977420572   -1.49800e+02   1.48808e-01
   @DF-UHF iter   1:  -149.59496320631871    2.05367e-01   2.58009e-02
   @DF-UHF iter   2:  -149.62349901753706   -2.85358e-02   6.68980e-03 DIIS
   @DF-UHF iter   3:  -149.62639942687878   -2.90041e-03   2.19285e-03 DIIS
   @DF-UHF iter   4:  -149.62689561367233   -4.96187e-04   5.99497e-04 DIIS
   @DF-UHF iter   5:  -149.62694151275420   -4.58991e-05   1.27338e-04 DIIS
   @DF-UHF iter   6:  -149.62694337910040   -1.86635e-06   1.65616e-05 DIIS
   @DF-UHF iter   7:  -149.62694340915198   -3.00516e-08   2.68990e-06 DIIS
   @DF-UHF iter   8:  -149.62694340999315   -8.41169e-10   2.61249e-07 DIIS

  DF guess converged.
  ...

   @UHF iter   9:  -149.62730705472407   -3.63645e-04   8.63697e-05 DIIS
   @UHF iter  10:  -149.62730737348096   -3.18757e-07   1.50223e-05 DIIS
   @UHF iter  11:  -149.62730738537113   -1.18902e-08   3.80466e-06 DIIS
   @UHF iter  12:  -149.62730738624032   -8.69193e-10   7.06634e-07 DIIS

The first set of iterations are from the DF portion of the computation, the
second set uses the exact (but much slower) Direct algorithm. Within the DF portion
of the computation, the zeroth-iteration uses a non-idempotent density matrix
obtained from the SAD guess, so the energy is unphysically low. However, the
first true iteration is quite close to the final DF energy, highlighting the
efficiency of the SAD guess. Pulay's DIIS procedure is then used to accelerate
SCF convergence, with the DF phase reaching convergence in eight true
iterations. When used together, SAD and DIIS are usually sufficient to converge
the SCF for all but the most difficult systems. Additional convergence
techniques are available for more difficult cases, and are detailed below. At
this point, the code switches on the requested Direct integrals technology, which
requires only four full iterations to reach convergence, starting from the DF
guess. This hybrid DF/Direct procedure can significantly accelerate SCF
computations requiring exact integrals.

After the iterations are completed, a number of one-electron properties are
printed, and some bookkeeping is performed to set up possible correlated
computations. Additional one-electron properties are available by increasing the
|globals__print| option. Also printed are the occupied and virtual orbital energies,
which are useful in elucidating the stability and reactivity of the system.

.. index::
   pair: SCF; theory

.. _`sec:scftheory`:

Theory
~~~~~~

The objective of Hartree--Fock (HF) Theory is to produce the optimized Molecular
Orbitals (MOs) :math:`\{\psi_i\}`,

.. math:: \psi_i(\vec x_1) = C_{\mu i} \phi_{\mu} (\vec x_1).

Here, :math:`\{\phi_{\mu}\}` are the basis functions, which, in |PSIfour| are
contracted Cartesian Gaussian functions often referred to as Atomic Orbitals
(AOs). The matrix :math:`C_{\mu i}` contains the MO coefficients, which are the
constrained variational parameters in Hartree--Fock. The molecular orbitals are
used to build the simplest possible antisymmetric wavefunction, a single Slater
determinant,

.. math:: | \Psi_0 \rangle =
    \frac{1}{\sqrt{N!}} \left | \begin{array}{cccc}
    \psi_1 (\vec x_1) & \psi_2(\vec x_1) & \ldots & \psi_N (\vec x_1) \\
    \psi_1 (\vec x_2) & \psi_2(\vec x_2) & \ldots & \psi_N (\vec x_2) \\
    \vdots & \vdots & \ddots & \vdots \\
    \psi_1 (\vec x_N) & \psi_2(\vec x_N) & \ldots & \psi_N (\vec x_N) \\
    \end{array}\right |

This form for the Hartree--Fock wavefunction is actually entirely equivalent to
treating the electron correlation as a mean field repulsion in
:math:`\mathbb{R}^6` instead of a more complicated effect in
:math:`\mathbb{R}^N`\ .

Considering the electronic Hamiltonian,

.. math:: \hat H = \sum_{i} -\frac{1}{2} \nabla_i^2 + \sum_{i} \sum_{A} -
    \frac{Z_A}{r_{iA}} + \sum_{i>j} \frac{1}{r_{ij}},

the Hartree--Fock energy is, by Slater's rules,

.. math:: E_{\mathrm{HF}} =
    \langle \Psi_0 | \hat H | \Psi_0 \rangle
    = \sum_{i} \langle i | \hat h | i \rangle
    + \frac 1 2 \sum_{i,j} [ii|jj] - [ij|ji]

.. math:: =
    D_{\mu\nu}^\alpha \left(H_{\mu\nu} + F_{\mu\nu}^{\alpha} \right)
    + D_{\mu\nu}^\beta \left(H_{\mu\nu} + F_{\mu\nu}^{\beta} \right)

Here :math:`H` is the AO-basis one-electron potential, encapsulating both
electron-nuclear attraction and kinetic energy,

.. math:: H_{\mu\nu} =
    \left(\mu \left| -\frac{1}{2} \nabla^2 + \sum_{A} -\frac{Z_A}{r_{1A}} \right
    | \nu \right),

:math:`D` is the AO-basis density matrix, build from the occupied orbital
coefficients,

.. math:: D_{\mu\nu}^{\alpha} =
    C_{\mu i}^{\alpha} C_{\nu i}^{\alpha},

and :math:`F` is the Fock matrix, which is the effective one-body potential at
the current value of the density,

.. math:: F_{\mu\nu}^{\alpha} = H_{\mu\nu}
    + \underbrace{\left(D_{\lambda\sigma}^{\alpha} + D_{\lambda\sigma}^{\beta}\right)
    (\mu\nu|\lambda\sigma)}_{J}
    + \underbrace{D_{\lambda\sigma}^{\alpha} (\mu\lambda|\sigma\nu)}_{K^{\alpha}}

Here the tensor :math:`(\mu\nu|\lambda\sigma)` is an AO Electron-Repulsion
Integral (ERI) in chemists' notation,

.. math:: (\mu\nu|\lambda\sigma) = \iint_{\mathbb{R}^6}
    \phi_{\mu} (\vec r_1)
    \phi_{\nu} (\vec r_1)
    \frac{1}{r_{12}}
    \phi_{\lambda} (\vec r_2)
    \phi_{\sigma} (\vec r_2)
    \ \mathrm{d}^3 r_1
    \ \mathrm{d}^3 r_2.

The MO coefficients are found as the generalized eigenvectors of the Fock Matrix,

.. math:: F^\alpha C^\alpha = S C^\alpha \epsilon^\alpha

The eigenvalues :math:`\epsilon` are the orbital energies, and the metric matrix
:math:`S` is the AO-basis overlap matrix

.. math:: S_{\mu\nu} = (\mu | \nu )

Note that the Fock Matrix depends on the density (both alpha and beta), and
therefore the orbitals. Because of this, SCF is a nonlinear procedure, which
terminates when the generating orbitals are self-consistent with the Fock matrix
they generate.

The formation of the Coulomb matrix :math:`J` and the exchange matrix
:math:`K^{\alpha}` dominate the computational effort of the SCF procedure. For
very large systems, diagonalization of the Fock matrix can also present a
significant hurdle.

.. _`sec:scfinput`:

Minimal Input
~~~~~~~~~~~~~

Minimal input for a Hartree--Fock computation is a molecule block, basis set
option, and a call to ``energy('scf')``::

    molecule {
    He
    }

    set basis sto-3g

    energy('scf')

This will run a Restricted Hartree--Fock (RHF) on neutral singlet Helium in
:math:`D_{2h}` spatial symmetry with a minimal ``STO-3G`` basis, 1.0E-6
energy and density convergence criteria (since single-point, see
:ref:`SCF Convergence & Algorithm <table:conv_scf>`), a DF ERI algorithm, symmetric
orthogonalization, DIIS, and a core Hamiltonian guess. For more
information on any of these options, see the relevant section below.

.. _`sec:scfsymm`:

Spin/Symmetry Treatment
~~~~~~~~~~~~~~~~~~~~~~~

|PSIfour| implements the most popular spin specializations of Hartree--Fock
theory, including:

Restricted Hartree--Fock (RHF) [Default]
  Appropriate only for closed-shell singlet systems, but twice as efficient
  as the other flavors, as the alpha and beta densities are constrained to be
  identical.
Unrestricted Hartree--Fock (UHF)
  Appropriate for most open-shell systems and fairly easy to converge.
  The spatial parts of the alpha and beta orbitals are fully independent of each
  other, which allows a considerable amount of flexibility in the wavefunction.
  However, this flexibility comes at the cost of spin symmetry; UHF wavefunctions
  need not be eigenfunctions of the :math:`\hat S^2` operator. The deviation of
  this operator from its expectation value is printed on the output file. If the
  deviation is greater than a few hundredths, it is advisable to switch to a
  ROHF to avoid this "spin-contamination" problem.
Restricted Open-Shell Hartree--Fock (ROHF)
  Appropriate for open-shell systems where spin-contamination is problem.
  Sometimes more difficult to converge, and assumes uniformly positive spin
  polarization (the alpha and beta doubly-occupied orbitals are identical).
Constrained Unrestricted Hartree--Fock (CUHF)
  A variant of ROHF that starts from a UHF ansatz and is therefore often
  easier to converge.

These can be invoked by the |scf__reference| keyword, which defaults to ``RHF``.
The charge and multiplicity may either be specified in the molecule definition::

    molecule h {
    0 2  # Neutral doublet
    H
    }

or, dynamically, by setting the relevant attributes in the Python molecule
object::

    h.set_molecular_charge(0)
    h.set_multiplicity(2)

Abelian spatial symmetry is fully supported in |PSIfour| and can be used to
obtain physical interpretation of the molecular orbitals, to assist in difficult
convergence cases, and, in some methods, to obtain significant performance
gains. The point group of the molecule is inferred when reading the molecule
section, and may be overridden by the :ref:`symmetry <sec:moleculeKeywords>` flag, as in::

    molecule h {
    0 2
    H
    symmetry c1
    }

or by the ``reset_point_group`` Python molecule attribute::

    h.reset_point_group('c2v')

During the SCF procedure, the occupation of orbitals is typically determined by
the Aufbau principal across all spatial symmetries. This may result in the
occupation shifting between iterations. If the occupations are known *a priori*,
they may be clamped throughout the procedure by using the |globals__docc| and
|globals__socc| options. For instance, all good quantum chemists know that
:math:`C_{2v}` water is
actually,::

    molecule h2o {
    0 1
    O
    H 1 1.0
    H 1 1.0 2 104.5
    }

    set {
    docc [3, 0, 1, 1]  # 1A1 2A1 1B1 3A1 1B2
    basis cc-pvdz
    }

    energy('scf')

.. _`sec:scfbrokensymm`:

Broken Symmetry
~~~~~~~~~~~~~~~

For certain problems, such diradicals, allowing the spin-up and spin-down
orbitals to differ in closed-shell computations can be advantageous;
this is known as symmetry breaking.  The resulting unrestricted wavefunction
will often provide superior energetics, due to the increased flexibility,
but it will suffer non-physical spin contamination from higher multiplicity states.
A convenient approach to break symmetry is to perform a UHF or UKS calculation
with the guess HOMO and LUMO orbitals mixed.
Mixing of the guess orbitals can be requested by setting the |scf__guess_mix|
keyword to true::

    set reference uhf
    set guess_mix true
    energy('scf')

.. _`sec:scflindep`:

Orthogonalization
~~~~~~~~~~~~~~~~~

One of the first steps in the SCF procedure is the determination of an
orthogonal basis (known as the OSO basis) from the atomic orbital
basis (known as the AO basis). The Molecular Orbital basis (MO basis)
is then built as a particular unitary transformation of the OSO
basis. In |PSIfour|, the determination of the OSO basis is
accomplished via either symmetric, canonical, or partial Cholesky
orthogonalization.

Symmetric orthogonalization uses the symmetric inverse square root of
the overlap matrix for the orthogonalization matrix. Use of symmetric
orthogonalization always yields the same number of OSO functions (and
thereby MOs) as AO functions. However, this may lead to numerical
problems if the overlap matrix has small eigenvalues, which may occur
for large systems or for systems where diffuse basis sets are used.

This problem may be avoided by using canonical orthogonalization, in
which an asymmetric inverse square root of the overlap matrix is
formed, with numerical stability enhanced by the elimination of
eigenvectors corresponding to very small eigenvalues. As a few
combinations of AO basis functions may be discarded, the number of
canonical-orthogonalized OSOs and MOs may be slightly smaller than the
number of AOs.

When the basis set is too overcomplete, the eigendecomposition of the
overlap matrix is no longer numerically stable. In this case the
partial Cholesky decomposition can be used to pick a subset of basis
functions that span a sufficiently complete set, see
[Lehtola:2019:241102]_ and [Lehtola:2020:032504]_. This subset can then
be orthonormalized as usual; the rest of the basis functions are
hidden from the calculation. The Cholesky approach allows reaching
accurate energies even in the presence of significant linear
dependencies [Lehtola:2020:04224]_.

In |PSIfour|, symmetric orthogonalization is used by default, unless
the smallest overlap eigenvalue falls below the user-supplied double
option |scf__s_tolerance|, which defaults to 1E-7. If the smallest
eigenvalue is below this cutoff, canonical orthogonalization is
forced, and all eigenvectors corresponding to eigenvalues below the
cutoff are eliminated.

If the eigendecomposition is detected to be numerically unstable - the
reciprocal condition number of the overlap matrix to be smaller than
the machine epsilon - the partial Cholesky decomposition is undertaken
until |scf__s_cholesky_tolerance|, which defaults to 1E-8.

Use of symmetric, canonical, and partial Cholesky orthogonalization
can be forced by setting the |scf__s_orthogonalization| option to
``SYMMETRIC``, ``CANONICAL``, or ``PARTIALCHOLESKY``,
respectively.

Note that in practice, the MOs and OSOs are built separately within
each irrep from the symmetry-adapted combinations of AOs known as
Unique Symmetry Orbitals (USOs).  For canonical orthogonalization,
this implies that the number of MOs and OSOs per irrep may be slightly
smaller than the number of USOs per irrep.

A contrived example demonstrating OSOs/MOs vs. AOs with symmetry is shown
below::

    molecule h2o {
    0 1
    O
    H 1 1.0
    H 1 1.0 2 104.5
    symmetry c2  # Two irreps is easier to comprehend
    }

    set {
    s_tolerance 0.0001      # Set an unreasonably tight
                            # tolerance to force canonical
    basis aug-cc-pv5z       # This diffuse basis will have
                            # small-ish eigenvalues for even H2O
    print 3
    }

    energy('scf')

Output::

  ==> Pre-Iterations <==

   -------------------------------------------------------
    Irrep   Nso     Nmo     Nalpha   Nbeta   Ndocc  Nsocc
   -------------------------------------------------------
     A        145     145       0       0       0       0
     B        142     142       0       0       0       0
   -------------------------------------------------------
    Total     287     287       5       5       5       0
   -------------------------------------------------------

  ...

  Minimum eigenvalue in the overlap matrix is 1.6888063568E-05.
  Using Canonical Orthogonalization with cutoff of 1.0000000000E-04.
  Irrep 0, 1 of 145 possible MOs eliminated.
  Irrep 1, 2 of 142 possible MOs eliminated.
  Overall, 3 of 287 possible MOs eliminated.

In this example, there are 287 AO basis functions after spherical harmonics are
applied. These are used to produce 287 symmetry adapted USOs, 145 of which are
assigned to irrep A, and 142 of which are assigned to irrep B. Within irrep A,
144 OSOs fall above the eigenvalue cutoff, and within irrep B 140 OSOs fall
above the eigenvalue cutoff. In total, 284 molecular orbitals are chosen from
287 AOs/USOs.

.. _`sec:scfguess`:

Initial Guess
~~~~~~~~~~~~~

In each step of the SCF procedure, a new Fock or Kohn--Sham potential is built
according to the previous density, following which the potential is diagonalized
to produce new molecular orbitals, from which a new density is computed. This
procedure is continued until either convergence is reached or a preset maximum
number of iterations is exceeded. Convergence is determined by both change in
energy and root-mean-square change in density matrix values, which must be below
the user-specified |scf__e_convergence| and |scf__d_convergence|, respectively.
The maximum number of iterations is specified by the |scf__maxiter| option. It
should be noted that SCF is a chaotic process, and, as such, often requires
careful selection of initial orbitals and damping during iterations to ensure
convergence. This is particularly likely for large systems, metallic systems,
multireference systems, open-shell systems, anions, and systems with diffuse
basis sets.

For initial orbital selection, several options are available. These include:

CORE
    Diagonalization of the core Hamiltonian, removing even mean-field electron
    repulsion. Simple, but often too far from the final solution for larger
    systems. This is the default for single atoms.
SAD [:term:`Default <GUESS (SCF)>`]
    Superposition of Atomic Densities. Builds the initial density as the
    spin-averaged sum of atomic UHF computations in the current basis. If an
    open-shell system, uniform scaling of the spin-averaged density matrices is
    performed. If orbitals are needed (*e.g.*, in density fitting), a partial
    Cholesky factorization of the density matrices is used. Often extremely
    accurate, particularly for closed-shell systems. This is the default for
    systems of more than one atom.
SADNO
    Natural orbitals from Superposition of Atomic Densities. Similar
    to the above, but it forms natural orbitals from the SAD density
    matrix to get proper orbitals which are used to start the
    calculation, see [Lehtola:2019:1593]_.
GWH
    A generalized Wolfsberg-Helmholtz modification of the core
    Hamiltonian matrix. Usually less accurate than the core guess: the
    latter is exact for one-electron systems, GWH is not; see
    [Lehtola:2019:1593]_).
HUCKEL
    An extended H\ |u_dots|\ ckel guess based on on-the-fly atomic UHF
    calculations alike SAD, see [Lehtola:2019:1593]_.
READ
    Read the previous orbitals from a ``wfn`` file, casting from
    one basis to another if needed. Useful for starting anion
    computations from neutral orbitals, or after small geometry
    changes. At present, casting from a different molecular point
    group is not supported.  This becomes the default for the second
    and later iterations of geometry optimizations.
SAP
    Superposition of Atomic Potentials. This is essentially a
    modification of the core Hamiltonian, which includes screening
    effects by using a radially screened effective atomic charge. The
    screening effects have been calculated at the complete basis set
    limit with finite-element calculations, see [Lehtola:2019:25945]_
    and [Lehtola:2020:012516]_. The guess and its implementation have
    been described in [Lehtola:2019:1593]_. The guess is evaluated on a
    DFT quadrature grid, so the guess energy depends slightly on the
    used DFT quadrature. The current implementation is based on
    exchange-only local density calculations that are but nanohartree
    away from the complete basis set limit [Lehtola:2020:012516]_.

These are all set by the |scf__guess| keyword. Also, an automatic Python
procedure has been developed for converging the SCF in a small basis, and then
casting up to the true basis. This can be done by adding
|scf__basis_guess| = SMALL_BASIS to the options list. We recommend the
3-21G or pcseg-0 basis for the small basis due to its efficient mix of flexibility and
compactness. An example of performing an RHF solution of water by SAD guessing
in a 3-21G basis and then casting up to cc-pVTZ is shown below::

    molecule h2o {
    0 1
    O
    H 1 1.0
    H 1 1.0 2 104.5
    }

    set {
    basis cc-pvtz
    basis_guess 3-21G
    guess sad
    }

    energy('scf')

.. _`sec:scfrestart`:

Restarting the SCF
~~~~~~~~~~~~~~~~~~

Reading orbital data from a previous calculations is done via the ``restart_file`` option,
where the actual file is a serialized ``wfn`` object (see :ref:`saving the wfn <sec:save_wfn>`)
By default, the orbital data file of the converged SCF(``psi.PID.name.180.npy``) is deleted
after |PSIfour| exits or the ``clean()`` function is called. The orbital guess is automatically
set to ``READ`` when ``restart_file`` is given a ``wfn`` file.
To write the orbitals after every iteration and keep the orbitals from the last iteration, the ``write_orbitals`` options is available: ::

  energy('scf', write_orbitals='my_mos'),

which writes a ``Wavefunction`` object converted (serialized) to a numpy file called ``my_mos.npy``.
The restart can then be done as follows: ::

  energy('scf', restart_file='my_mos')

Specifying the ``.npy`` suffix when writing and reading restart files is optional.

Alternatively, the restart can also be done from any previously saved ``wfn`` object. ::

  energy, scf_wfn = energy('scf',return_wfn=True)
  scf_wfn.to_file('my_wfn')
  energy('scf', restart_file='my_wfn')


For advanced users manipulating or writing custom wavefunction files, note
that |PSIfour| expects the numpy file on disk to have the ``.npy`` extension, not, e.g., `.npz`.


.. index:: DIIS, MOM, damping

.. _`sec:scfconv`:

Convergence Stabilization
~~~~~~~~~~~~~~~~~~~~~~~~~

A summary of Psi's supported convergence stabilization techniques is presented below:

DIIS [On by Default]
    DIIS uses previous iterates of the Fock matrix together
    with an error criterion based on the orbital gradient to produce an informed
    estimate of the next Fock Matrix. DIIS is almost always necessary to converge
    the SCF procedure and is therefore turned on by default. In rare cases, the
    DIIS algorithm may need to be modified or turned off altogether, which may be
    accomplished via :term:`options <DIIS (SCF)>`.
ADIIS [On by Default]
    ADIIS uses previous iterates of the Fock and density matrices to produce an
    informed estimate of the next Fock matrix. ADIIS estimates are based on minimizing
    an energy estimate rather than zeroing the residual, so this performs best in the early
    iterations. By default, Psi will start using ADIIS before blending the ADIIS step with
    the DIIS step, eventually using the pure DIIS step. The closely-related EDIIS procedure
    may be used instead by setting |scf__scf_initial_accelerator|. This is formally identical
    to ADIIS for HF, but the methods will differ for more general DFT.
MOM [Off by Default]
    MOM was developed to combat a particular class of convergence failure:
    occupation flipping. In some cases, midway though the SCF procedure, a partially
    converged orbital which should be occupied in the fully-optimized SCF solution
    has a slightly higher orbital eigenvalue than some other orbital which should be
    destined to be a virtual orbital. This results in the virtual orbital being
    spuriously occupied for one or more iterations. Sometimes this resolves itself
    without help, other times the occupation flips back and forth between two, four,
    or more orbitals. This is typically visible in the output as a non-converging
    SCF which eventually settles down to steady oscillation between two (or more)
    different total energies. This behavior can be ameliorated by choosing occupied
    orbitals by "shape" instead of by orbital eigenvalue, *i.e.*, by choosing the set
    of new orbitals which looks most like some previously known "good" set.  The
    "good" set is typically the occupied orbitals from one of the oscillating
    iterations with the lowest total energy. For an oscillating system where the
    lowest total energy occurs on iterations :math:`N,N+2,\ldots`, invoking
    |scf__mom_start| ``N`` can often rescue the convergence of the SCF. MOM can
    be used in concert with DIIS, though care should be taken to not turn MOM on
    until the oscillatory behavior begins.
Damping [Off by Default]
    In some cases, a static mixing of Fock Matrices from adjacent
    iterations can quench oscillations. This mixing, known as
    "damping" can be activated by setting the
    |scf__damping_percentage| keyword to a nonzero percent. Damping is
    turned off when the DIIS error is smaller than
    |scf__damping_convergence|.
Level shifting [Off by default]
    A commonly used alternative to damping is to use level shifting,
    which decreases the mixing of occupied and unoccupied orbitals in
    the SCF update by moving the unoccupied orbitals up in energy. It
    can be shown that the SCF procedure always converges with a
    suitably large level shift; however, the larger the shift is, the
    slower the convergence becomes, and the calculation may end up
    converging onto a higher lying SCF solution. Because of this, in
    practice level shifting is most useful in the initial phase of the
    calculation to reduce the orbital error enough for DIIS to work
    well. The level shift is controlled by the parameter
    |scf__level_shift|, and it is turned off when the DIIS error is
    smaller than |scf__level_shift_cutoff|. Reasonable values for
    the shift and convergence threshold are 5.0 and 1e-2,
    respectively.
SOSCF [Off by Default]
    See :ref:`sec:soscf`

.. _`sec:scferi`:

ERI Algorithms
~~~~~~~~~~~~~~

The key difficulty in the SCF procedure is treatment of the four-index ERI
contributions to the Fock Matrix. A number of algorithms are available in
|PSIfour| for these terms. The algorithm is selected by the |globals__scf_type|
keyword. Some such algorithms consist of a single algorithm applied to
the consutruction of both the Coulomb and Exchange constructions
to the Fock Matrix, such as follows:

PK [:ref:`Default <table:conv_scf>`]
    An out-of-core, presorted algorithm using exact ERIs. Quite fast for a
    zero-error algorithm if enough memory is available. Integrals are
    generated only once, and symmetry is utilized to reduce number of
    integrals.
OUT_OF_CORE
    An out-of-core, unsorted algorithm using exact ERIs. Overcomes the
    memory bottleneck of the current PK algorithm. Integrals are generated
    only once, and symmetry is utilized to reduce number of integrals.
DIRECT
    A threaded, sieved, integral-direct algorithm, with full permutational
    symmetry. This algorithm is brand new, but seems to be reasonably fast
    up to 1500 basis functions, uses zero disk (if DF pre-iterations are
    turned off), and can obtain significant
    speedups with negligible error loss if |scf__ints_tolerance|
    is set to 1.0E-8 or so.
DF [:ref:`Default <table:conv_scf>`]
    A density-fitted algorithm designed for computations with thousands of
    basis functions. This algorithm is highly optimized, and is threaded
    with a mixture of parallel BLAS and OpenMP. Note that this algorithm
    should use the -JKFIT series of auxiliary bases, *not* the -RI or
    -MP2FIT bases. The default guess for auxiliary basis set should work
    for most bases, otherwise the |scf__df_basis_scf| keyword can
    be used to manually specify the auxiliary basis.  This algorithm is
    preferred unless either absolute accuracy is required
    [:math:`\gtrsim`\ CCSD(T)] or a -JKFIT auxiliary basis is unavailable
    for the orbital basis/atoms involved.
CD
    A threaded algorithm using approximate ERIs obtained by Cholesky
    decomposition of the ERI tensor.  The accuracy of the Cholesky
    decomposition is controlled by the keyword |scf__cholesky_tolerance|.
    This algorithm is similar to the DF algorithm, but it is not suitable
    for gradient computations.  The algorithm to obtain the Cholesky
    vectors is not designed for computations with thousands of basis
    functions.

|PSIfour| also features the capability to use arbitrary combinations 
of specialized algorithms that construct either the Coulomb
or the Exchange matrix separately. These algorithms can be mixed and
matched arbitrarily to optimize performance based on the calculation
in question. In general, such combinations of algorithms display lower
scaling then their combined counterparts, and thus perform better 
with larger systems. Such combinations of algorithms can be utilized 
by setting the |globals__scf_type| keyword to ``J_alg+K_alg``, 
where *J_alg* is the name of the separate Coulomb construction algoritm to use, 
and *K_alg* is the name of the separate Exchange construction algorithm to use.
More information about the machinery of separate Coulomb and Exchange matrix
construction algorithms in Psi4 can be found in the :ref:`sec:compositejk` section.

Specialized algorithms available to construct the Coulomb term
separately are as follows:

DIRECTDFJ
    An integral-direct algorithm constructing the Coulomb term based on [Weigend:2002:4285]_
    The DIRECTDFJ algorithm combines the benefits of integral-direct SCF approaches 
    with that of density-fitting. Specifically, DFJ utilizes no I/O and displays 
    strong performance with large system size through a combination of 
    effective parallelization and utilization of density-fitting to minimize 
    ERI computational cost. See the :ref:`sec:scfddfj` section for more information.

Specialized algorithms available to construct the Exchange term
separately are as follows:

COSX
    An algorithm based on the semi-numerical "chain of spheres exchange" (COSX)
    approach described in [Neese:2009:98]_. The coulomb term is computed with a
    direct density-fitting algorithm. The COSX algorithm uses no I/O, scales
    well with system size, and requires minimal memory, making it ideal for
    large systems and multi-core CPUs. See :ref:`sec:scfcosx` for more information.
LINK
    An implementation of the linear-scaling "Linear Exchange" (LinK)
    algorithm described in [Ochsenfeld:1998:1663]_. The LINK algorithm provides 
    many of the benefits of integral-direct SCF algorithms, including no disk I/O, 
    low memory usage, and effective parallelization. Additionally, the
    LINK implementation scales well with system size 
    while simultaneously providing a formally-exact computation of the 
    Exchange term. See :ref:`sec:scflink` for more information.

In some cases the above algorithms have multiple implementations that return
the same result, but are optimal under different molecules sizes and hardware
configurations. Psi4 will automatically detect the correct algorithm to run and
only expert users should manually select the below implementations. The DF
algorithm has the following two implementations

MEM_DF
    A DF algorithm optimized around memory layout and is optimal as long as
    there is sufficient memory to hold the three-index DF tensors in memory. This
    algorithm may be faster for builds that require disk if SSDs are used.
DISK_DF
    A DF algorithm (the default DF algorithm before Psi4 1.2) optimized to
    minimize Disk IO by sacrificing some performance due to memory layout.

Note that these algorithms have both in-memory and on-disk options, but
performance penalties up to a factor of 2.5 can be found if the incorrect
algorithm is chosen. It is therefore highly recommended that the keyword "DF"
be selected in all cases so that the correct implementation can be selected by
|PSIfours| internal routines. Expert users can manually switch between MEM_DF and
DISK_DF; however, they may find documented exceptions during use as several
post SCF algorithms require a specific implementation. Additionally, expert users 
can manually switch between the in-memory and on-disk options *within* MEM_DF or DISK_DF using 
the |scf__scf_subtype| option. Using ``SCF_SUBTYPE = AUTO``, where |PSIfour| 
automatically selects the in-memory or on-disk option for MEM_DF/DISK_DF based on memory and molecule, is the default 
and recommended option. However, the in-memory or on-disk algorithms for MEM_DF and DISK_DF can be forced by using
``SCF_SUBTYPE = INCORE`` or ``SCF_SUBTYPE = OUT_OF_CORE``, respectively.
Note that an exception will be thrown if 
``SCF_SUBTYPE = INCORE`` is used without allocating sufficient memory to 
|PSIfour|.

For some of these algorithms, Schwarz and/or density sieving can be used to
identify negligible integral contributions in extended systems. To activate
sieving, set the |scf__ints_tolerance| keyword to your desired cutoff
(1.0E-12 is recommended for most applications). To choose the type of sieving, set 
the |globals__screening| keyword to your desired option. For Schwarz screening, set it
to ``SCHWARZ``, for CSAM, ``CSAM``, and for density matrix-based screening, ``DENSITY``.

SCHWARZ
    Uses the Cauchy-Schwarz inequality to calculate an upper bounded value of a shell quartet,

.. math:: (PQ|RS) <= \sqrt{(PQ|PQ)(RS|RS)}

CSAM
    An extension of the Schwarz estimate that also screens over the long range 1/r operator, described in [Thompson:2017:144101]_.

DENSITY
    An extension of the Schwarz estimate that also screens over elements of the density matrix.
    For the RHF case, described in [Haser:1989:104]_

.. math:: CON(PQ|RS) <= \sqrt{(PQ|PQ)(RS|RS)} \cdot DCON(PQ, RS)

.. math:: DCON(PQ, RS) = max(4D_{PQ}, 4D_{RS}, D_{PR}, D_{PS}, D_{QR}, D_{QS})

When using density-matrix based integral screening, it is useful to build the J and K matrices
incrementally, also described in [Haser:1989:104]_, using the difference in the density matrix between iterations, rather than the
full density matrix. To turn on this option, set |scf__incfock| to ``true``.

We have added the automatic capability to use the extremely fast DF
code for intermediate convergence of the orbitals, for |globals__scf_type|
``DIRECT``. At the moment, the code defaults to cc-pVDZ-JKFIT as the
auxiliary basis, unless the user specifies |scf__df_basis_scf| manually. For
some atoms, cc-pVDZ-JKFIT is not defined, so a very large fitting basis of last
resort will be used.
To avoid this, either set |scf__df_basis_scf| to an auxiliary
basis set defined for all atoms in the system, or set |scf__df_scf_guess|
to false, which disables this acceleration entirely.

.. _`sec:compositejk`:

"Composite" Coulomb and Exchange Matrix Construction Algorithms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Alongside the more traditional SCF algorithms in Psi4, in which both the Coulomb (J) and the Exchange (K) contributions to the Fock matrix are computed
simultaneously using a single algorithm, "composite" algorithms also exist for the construction of the Fock matrix. These "composite" algorithms
consist of a collection of SCF algorithms that construct one of J or K, but not both. Two such algorithms, one to construct J and one 
to construct K, can then be combined to form the full Fock matrix during the SCF procedure. 

Generally, algorithms that are specifically designed to construct one of J or K exploit properties specific to J or K to improve performance.
Exploitation of these properties can lead to either a reduced prefactor or a reduced algorithmic scaling for computing that specific
matrix. Since composite J+K algorithm combinations often involve algorithms which reduce the algorithmic scaling of either J or K (or both), 
composite algorithm combinations tend to display better scaling properties than single algorithms which compute J and K simultaneously. 
However, building J and K separately via a composite algorithm combination requires the redundant recomputation of ERIs, a facet which is not 
necessary for combined Fock build algorithms. The overall consequence, then, is that composite combinations of separate J and K construction 
algorithms tend to perform better than combined Fock build algorithms for sufficiently large system sizes, while performing worse for smaller systems. 

Multiple separate J and K construction algorithms to construct the Fock matrix in a composite fashion are available for use in Psi4, and they can be 
mixed and matched with each other freely. The algorithms availble for use in a composite fashion are as follows:

Separate J Construction: DIRECTDFJ 

Separate K Construction: COSX, LINK

Two composite algorithms can be combined with each other by setting the |globals__scf_type| keyword to ``J_alg+K_alg``, where *J_alg* is an available
algorithm for the separate construction of J, and *K_alg* is an available algorithm for the separate construction of K. For example, if one wanted
to run an integral-direct density-fitted algorithm for Coulomb construction (i.e., DIRECTDFJ) in conjunction with the Linear Exchange algorithm
for Exchange construction (i.e., LINK), one would use the |globals__scf_type| keyword ``DIRECTDFJ+LINK``.  

.. _`sec:scfddfj`:

Integral-Direct Density-Fitted Coulomb Construction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Resolution of the Identity (RI) can be used to decompose the normally 4-center ERI tensor into a combination of 3-center and 2-center components.
By reducing the dimensionality of the ERI tensor, application of the RI (often referred to as density-fitting, or DF) can be used to greatly speed up
SCF calculations. The reduction in ERI tensor rank also makes DF an appealing option for conventional SCF calculations, where the ERIs are stored 
in core or on disk. However, even when using DF, I/O becomes a significant bottleneck for systems of a sufficient size when performing conventional SCF
calculations. In principle, though, DF approaches can be utilized in an integral-direct context, gaining the benefits of DF methods without suffering the
I/O bottlenecks that conventional DF methods will eventually run into. One such approach, outlined by Weigend in [Weigend:2002:4285]_,
is available for use in Psi4 for the separate construction of the Coulomb contribution to the Fock matrix.  This implementation can be used alongside 
Psi4's separate Exchange construction algorithms for composite Fock matrix construction by using the keyword DIRECTDFJ as the Coulomb construction 
algorithm when specifying |globals__scf_type| to use a composite algorithm combination (``DIRECTDFJ+K_alg``). 

DIRECTDFJ supports multiple capabilities to improve performance. Specifically, DIRECTDFJ allows for a combination of density-matrix based ERI 
screening (set |globals__screening| to ``DENSITY``) and incremental Fock matrix construction (set |scf__incfock| to ``TRUE``). These two, when combined,
enable more aggressive screening of ERI contributions to the Coulomb matrix and thus greatly improve performance.

.. _`sec:scfcosx`:

COSX Exchange
~~~~~~~~~~~~~

The semi-numerical COSX algorithm described in [Neese:2009:98]_ evaluates
two-electron ERIs analytically over one electron coordinate and numerically
over the other electron coordinate, and belongs to the family of pseudospectral 
methods originally suggested by Friesner. In COSX, numerical integration is performed on standard
DFT quadrature grids, which are described in :ref:`sec:dft`.
Both the accuracy of the COSX algorithm and also the computational
cost are directly determined by the size of the integration grid, so selection
of the grid is important. This COSX implementation uses two separate grids.
The SCF algorithm is first converged on a smaller grid, followed by a final SCF
iteration on a larger grid. This results in numerical errors comparable to
performing the entire SCF on the expensive larger grid at a computational cost
much closer to the smaller grid. The size of the initial grid is controlled by the
keywords |scf__cosx_radial_points_initial| and |scf__cosx_spherical_points_initial|.
The final grid is controlled by |scf__cosx_radial_points_final| and
|scf__cosx_spherical_points_final|. The defaults for both grids aim to balance
cost and accuracy.

Screening thresholds over integrals, densities, and basis extents are set
with the |scf__cosx_ints_tolerance|, |scf__cosx_density_tolerance|, and
|scf__cosx_basis_tolerance| keywords, respectively. |scf__cosx_ints_tolerance|
is the most consequential of the three thresholds in both cost and accuracy.
This keyword determines screening of negligible one-electron integrals.
|scf__cosx_density_tolerance| controls the threshold for significant
shell pairs in the density matrix. Lastly, |scf__cosx_basis_tolerance| is
a cutoff for the value of basis functions at grid points. This keyword is
used to determine the radial extent of the each basis shell, and it is the
COSX analogue to |scf__dft_basis_tolerance|.

The |scf__incfock| keyword (defaults to ``false``) increases performance
by constructing the Fock matrix from differences in the density matrix, which
are more amenable to screening. This option is disabled by default because of
potential SCF convergence issues, particularly when using diffuse basis functions.
The |scf__cosx_overlap_fitting| keyword (defaults to ``true``) reduces numerical
integration errors using the method described in [Izsak:2011:144105]_ and is
always recommended.

.. _`sec:scflink`:

Linear Exchange
~~~~~~~~~~~~~~~

Large SCF calculations can benefit from specialized screening procedures that further reduce the scaling of the ERI contribution to the Fock matrix.
LinK, the linear-scaling exchange method described in [Ochsenfeld:1998:1663]_, is available in Psi4 in conjunction with composite algorithms that build J (|globals__scf_type| set to ``J_alg+LINK``).
LinK achieves linear-scaling by exploiting shell pair sparsity in the density matrix and overlap sparsity between shell pairs. Specifically, LinK exploits the fact that the Exchange term
requires only a linear-scaling number of significant elements through reformulating the
shell quartet screening process to scale linearly with system size.
LinK is most competitive when used with non-diffuse orbital basis sets, since orbital and density overlaps decay slower with diffuse functions.
LinK is especially powerful when combined with density-matrix based ERI screening (set |globals__screening| to ``DENSITY``) and incremental Fock builds (set |scf__incfock| to ``TRUE``), which decrease the number of significant two-electron integrals to calculate.

To control the LinK algorithm, here are the list of options provided.
  
  |scf__linK_ints_tolerance|: The integral screening tolerance used for sparsity-prep in the LinK algorithm. Defaults to the |scf__ints_tolerance| option.

.. index::
    single: SOSCF

.. _`sec:soscf`:

Second-order Convergence
~~~~~~~~~~~~~~~~~~~~~~~~

Second-order convergence takes into account both the gradient and Hessian to
take a full Newton step with respect to the orbital parameters. This results in
quadratic convergence with respect to density for SCF methods. For cases where
normal acceleration methods either fail or take many iterations to converge,
second-order can reduce the total time to solution.

Solving second-order (SO) methods exactly would require an inversion of the
orbital Hessian (an expensive :math:`\mathbb{N}^6` operation); however, these
equations are normally solved iteratively where each iteration costs the same
as a normal Fock build (:math:`\mathbb{N}^4`). The overall SOSCF operation is
thus broken down into micro- and macroiterations where the microiterations
refer to solving the SOSCF equations and macroiterations are the construction
of a new Fock matrix based on the orbitals from a SOSCF step.

SOSCF requires that all elements of the gradient to be less than one before the
method is valid. To this end, pre-SOSCF SCF iterations use normal
gradient-based extrapolation procedures (*e.g.*, DIIS) until the gradient
conditions are met. Note that while the total number of macroiterations will be
less for SOSCF than gradient-based convergence acceleration, the cost of solving
the microiterations typically results in the overall cost being greater for
SOSCF than for gradient-based methods. Therefore, SOSCF should only be used if
it is difficult to locate a stable minimum.

SOSCF is available for all HF and DFT references with the exception of meta-
GGA functionals. To enable, set the option |scf__soscf| to ``true``.
Additional options to modify the number of microiterations taken are as
follows:

    |scf__soscf_start_convergence|: when to start SOSCF based on the current density RMS

    |scf__soscf_max_iter|: the maximum number of SOSCF microiterations per macroiteration

    |scf__soscf_conv|: the relative convergence tolerance of the SOSCF microiterations

    |scf__soscf_print|: option to print the microiterations or not


.. _`sec:scfstability_doc`:

Stability Analysis
~~~~~~~~~~~~~~~~~~

SCF algorithms attempt to minimize the gradient of the energy with respect
to orbital variation parameters. At convergence, the gradient should be approximately zero
given a convergence criterion. Although this is enough to make sure the SCF converged to a
stationary point, this is not a sufficient condition for a minimal SCF solution. It may be
a saddle point or a maximum.

To ensure that a minimum has been found, the electronic Hessian, *i.e.* the matrix of second
derivatives of the energy with respect to orbital variation parameters, must be computed.
If one or more eigenvalues of the electronic Hessian are negative, the SCF solution is not a minimum.
In that case, orbital parameters can be varied along the lowest Hessian eigenvector to lower the energy.

Orbital variation parameters are usually constrained. For example, in RHF the
spatial parts of the :math:`\alpha` and :math:`\beta` orbitals are the same. In
UHF, the orbital coefficients are usually constrained to be real. A stability analysis
can check whether a lower SCF solution exists while respecting the constraints of the original
solution; this is an internal instability. If one or more constraints have to be relaxed to reach
a lower-energy solution, there is an external instability. In |PSIfour|, the only external instability
that can be checked at present is the RHF :math:`\rightarrow` UHF one.

Currently, two algorithms exist in |PSIfour| for stability analysis: the original
Direct Inversion and the newly implemented Davidson algorithms. We will first describe
options common to both algorithms. To request a stability analysis at the end of the SCF,
set the keyword |scf__stability_analysis|. Value ``CHECK`` only computes the electronic
Hessian eigenvalue and checks if an actual SCF minimum has been found, while value ``FOLLOW``
rotates the converged orbitals along the lowest eigenvector, then invokes the SCF
procedure again to lower the energy. In case the minimization does not succeed
or ends up on the same unstable solution, you can tune the scale factor for the orbital
rotation through the keyword |scf__follow_step_scale|.
The rotation angle is :math:`\frac{\pi}{2}\mbox{ } \cdot` (|scf__follow_step_scale|). The default value of
0.5 usually provides a good guess, and modification is only recommended in difficult cases.
The default behavior for the stability code is to stop after trying to reoptimize the orbitals once
if the instability still exists. For more attempts, set |scf__max_attempts|;
the default value of 1 is recommended. In case the SCF ends up in the same minimum, modification
of |scf__follow_step_scale| is recommended over increasing |scf__max_attempts|.

.. note:: Setting the option |scf__stability_analysis| to ``FOLLOW`` is only avalible for UHF. When using
   RHF and ROHF instabilities can be checked, but not followed. If you want to attempt to find a lower energy solution
   you should re-run the calculation with |scf__reference| set to ``UHF``.

The main algorithm available in |PSIfour| is the Direct Inversion algorithm. It can *only*
work with |globals__scf_type| ``PK``, and it explicitly builds the full electronic Hessian
matrix before explicitly inverting it. As such, this algorithm is very slow and it should
be avoided whenever possible. Direct Inversion is automatically invoked if the newer algorithm
is not available.

The Davidson algorithm for stability analysis was implemented recently.
Only the lowest eigenvalues of the electronic Hessian are computed, and Hessian-vector
products are computed instead of the full Hessian. This algorithm is thus
much more efficient than the Direct Inversion, but at present, it is only available for UHF :math:`\rightarrow` UHF stability
analysis. The capabilities of both algorithms are summarized below:

.. _`table:stab_methods`:

.. table:: Stability analysis methods available in |PSIfour|

    +------------------+------------------+----------------------------------------------+---------------------------+---------------------+
    |     Algorithm    | |scf__reference| |     Stability checked                        | |scf__stability_analysis| | |globals__scf_type| |
    +==================+==================+==============================================+===========================+=====================+
    |                  |       RHF        | Internal, External (:math:`\rightarrow` UHF) | ``CHECK``                 |   PK only           |
    +                  +------------------+----------------------------------------------+---------------------------+---------------------+
    | Direct Inversion |       ROHF       | Internal                                     | ``CHECK``                 |   PK only           |
    +------------------+------------------+----------------------------------------------+---------------------------+---------------------+
    |   Davidson       |       UHF        | Internal                                     | ``CHECK`` or ``FOLLOW``   |   Anything          |
    +------------------+------------------+----------------------------------------------+---------------------------+---------------------+

The best algorithm is automatically selected, *i.e.* Davidson for UHF :math:`\rightarrow` UHF and Direct Inversion otherwise.

In addition to the options available for Direct Inversion, the Davidson algorithm can automatically
adapt |scf__follow_step_scale| to find a new SCF minimum. If |scf__max_attempts| > 1, additional attempts
will automatically increment |scf__follow_step_scale| by 0.2 every time the SCF falls back to the previously
found unstable minimum. The increment can be adjusted by setting |scf__follow_step_increment|.
The default value is 0.2; adjust if needed to try different values of |scf__follow_step_scale| in a single computation.

The Davidson solver for the eigenvalues is controlled through several keywords.

  |scf__solver_maxiter|: maximum number of iterations

  |scf__solver_convergence|: eigenvector convergence threshold

  |scf__solver_n_root|: Solve for N eigenvectors in each irreducible representation

  |scf__solver_roots_per_irrep|: The number of eigenvectors to solve in each irreducible representation. An array of as many integers as there are irreducible representations.

  |scf__solver_n_guess|: Use N guess vectors, this needs to be larger than the number of roots so that the lowest ones can be captured reliably. Defaults to 4 guess vectors per root.

.. warning:: Prior to Dec 2022, v1.7, Psi4 had a different set of keywords controlling instability analysis, and those were included in the CPHF module rather than the SCF module.
.. warning:: Extending Davidson instability analysis to Kohn-Sham references is under development. As of 1.7, only LDA functionals are currently supported.

In case convergence problems are encountered during the Davidson procedure, file a bug report.
If the solver seems to converge on the wrong eigenvalue, try increasing |scf__solver_n_guess|.
Otherwise, if the solver is almost converged but reaches the maximum number of iterations, try increasing
|scf__solver_maxiter|.


.. _`sec:scf-ecps`:

Effective core potentials (ECPs)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

|PSIfour| supports the use of effective core potentials to describe the
innermost electrons in heavy elements.
ECPs are only available if |PSIfour| is compiled with the :ref:`LibECPInt <cmake:ecpint>` library.
If a basis set is designed to use an
effective core potential, the ECP definition should be simply placed alongside
the orbital basis set definition, *c.f.* :ref:`sec:basissets-ecps`.  All
information related to the definition and number of core electrons will
automatically be detected and no further input is required to use the
ECP-containing basis set.  See :srcsample:`scf-ecp` and :srcsample:`dfmp2-ecp`
for examples of computations with ECP-containing basis sets.

.. warning:: Prior to May 2022, v1.6, Psi4 used a built-in ECP code. Analytic derivatives of ECPs were not available. The HF and DFT derivatives were implemented in a semi-numerical scheme, where numerical ECP gradients were added to analytic SCF gradients. For post-SCF methods, the entire gradient computation needed to be run as finite difference of energies.

.. warning:: As of May 2022, v1.6, Psi4 uses the LibECPInt library, and analytic derivatives and Hessians of ECPs are available. Analytic derivatives of molecular systems including ECPs should be available whenever the method has analytic derivatives, but these have so far only been verified for HF and DFT.

.. warning:: ECPs have not been tested with projected basis set guesses or with FI-SAPT calculations.  If you require this functionality, please contact the developers on GitHub and/or the `forum <http://forum.psicode.org>`_.

.. _`sec:scfqmmm`:

External potentials and QM/MM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the implementation of :ref:`EFP <sec:libefp>` for accurate QM/MM
computations, |PSIfour| can perform more rudimentary QM/MM procedures via the
|scf__extern| keyword.  The following snippet, extracted from the
:srcsample:`extern1` test case, demonstrates its use for a TIP3P external potential::

    import numpy as np
    external_potentials = [
        [-0.834, np.array([1.649232019048,0.0,-2.356023604706]) / psi_bohr2angstroms],
        [ 0.417, np.array([0.544757019107,0.0,-3.799961446760]) / psi_bohr2angstroms],
        [ 0.417, np.array([0.544757019107,0.0,-0.912085762652]) / psi_bohr2angstroms]]

    gradient('scf', external_potentials=external_potentials)

The ``external_potentials`` array has three rows for three separate
particles, and it is passed to the SCF code on the last line. The
rows are composed of the atomic charge, x coordinate, y coordinate,
and z coordinate in that order. The atomic charge and coordinates are
specified in atomic units, [e] and [a0]. Add as many particle rows as
needed to describe the full MM region.

.. caution:: In |PSIfour| previous to Spring 2022 and v1.6, setting an
   external potential like the above looked like ::

    Chrgfield = QMMM()
    Chrgfield.extern.addCharge(-0.834, 1.649232019048, 0.0, -2.356023604706)
    Chrgfield.extern.addCharge( 0.417, 0.544757019107, 0.0, -3.799961446760)
    Chrgfield.extern.addCharge( 0.417, 0.544757019107, 0.0, -0.912085762652)
    psi4.set_global_option_python('EXTERN', Chrgfield.extern)

    gradient('scf')

   The main differences are that (1) the specification of
   charge locations in the old way used the units of the active
   molecule, whereas the new way always uses Bohr and (2) the
   specification of the charge and locations in the old way used the
   :py:class:`psi4.driver.QMMM` class directly and added one charge
   per command, whereas the new way consolidates all into an array and
   passes it by keyword argument to the calculation.

   The successor to the :py:class:`psi4.driver.QMMM` class,
   :py:class:`psi4.driver.QMMMbohr`, is operable, but it is discouraged
   from being used directly.

To run a computation in a constant dipole field, the |scf__perturb_h|,
|scf__perturb_with| and |scf__perturb_dipole| keywords can be used.  As an
example, to add a dipole field of magnitude 0.05 a.u. in the y direction and
0.1 a.u. in the z direction, we can use the following keywords::

    set perturb_h true
    set perturb_with dipole
    set perturb_dipole [ 0, 0.05, 0.1 ]

Note that if any specified fields do not fall along a symmetry axis, the
symmetry of the calculation should be reduced accordingly; if in doubt run the
calculation in C1 symmetry.  For examples of SCF and MP2 calculations in an
external field, see :srcsample:`scf7` and :srcsample:`dfmp2-grad5`.

.. _`sec:scfdefault`:

Convergence and Algorithm Defaults
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _`table:conv_scf`:

.. table:: SCF algorithm and convergence criteria defaults by calculation type [#f1]_

    +--------------------+--------------------+----------------------+----------------------+---------------------+
    | *Ab Initio* Method | Calculation Type   | |scf__e_convergence| | |scf__d_convergence| | |globals__scf_type| |
    +====================+====================+======================+======================+=====================+
    | SCF of HF or DFT   | energy             | 6                    | 6                    | DF                  |
    +                    +--------------------+----------------------+----------------------+                     +
    |                    | optimization       | 8                    | 8                    |                     |
    +                    +--------------------+----------------------+----------------------+                     +
    |                    | frequency [#f7]_   | 8                    | 8                    |                     |
    +--------------------+--------------------+----------------------+----------------------+---------------------+
    | SCF of post-HF     | energy             | 8                    | 8                    | PK [#f3]_           |
    +                    +--------------------+----------------------+----------------------+                     +
    |                    | optimization       | 10                   | 10                   |                     |
    +                    +--------------------+----------------------+----------------------+                     +
    |                    | frequency [#f7]_   | 10                   | 10                   |                     |
    +                    +--------------------+----------------------+----------------------+                     +
    |                    | CC property [#f2]_ | 10                   | 10                   |                     |
    +--------------------+--------------------+----------------------+----------------------+---------------------+

.. _`table:conv_corl`:

.. table:: Post-SCF convergence criteria defaults by calculation type [#f4]_

    +--------------------+--------------------+----------------------+-------------------------+
    | *Ab Initio* Method | Calculation Type   | E_CONVERGENCE [#f5]_ | R_CONVERGENCE [#f6]_    |
    +====================+====================+======================+=========================+
    | post-HF of post-HF | energy             | 6                    |                         |
    +                    +--------------------+----------------------+-------------------------+
    |                    | optimization       | 8                    |                         |
    +                    +--------------------+----------------------+-------------------------+
    |                    | frequency [#f7]_   | 8                    |                         |
    +                    +--------------------+----------------------+-------------------------+
    |                    | CC property [#f2]_ | 8                    |                         |
    +--------------------+--------------------+----------------------+-------------------------+

.. rubric:: Footnotes

.. [#f1] Note that this table applies only the SCF module,
   not to the final convergence criteria for post-HF methods or to methods
   that use an alternate starting point, like MCSCF. SAPT computations, too,
   set tighter values.

.. [#f2] This applies to properties computed through the :py:func:`~psi4.driver.properties` function.

.. [#f3] Post-HF methods that do not rely upon the usual 4-index AO integrals use a
   density-fitted SCF reference. That is, for DF-MP2 and SAPT, the default |globals__scf_type| is DF.

.. [#f4] Note that this table applies to the final convergence criteria for
   all the post-SCF modules that define a |ccenergy__e_convergence| keyword.

.. [#f5] The E_CONVERGENCE keyword is implemented for most post-SCF modules.
   See a list beginning at |ccenergy__e_convergence|.

.. [#f6] The R_CONVERGENCE keyword places a convergence check on an internal
   residual error measure and is implemented for several post-SCF
   modules (see list beginning at |ccenergy__r_convergence|). It is defined
   according to the quantum chemical method and so its default value is set
   by each module individually.

.. [#f7] For frequency computations by finite difference of energies,
   convergence criteria are tightened further still to 10 for
   |scf__e_convergence| and |scf__d_convergence| for SCF of HF or DFT, 11
   for |scf__e_convergence| and |scf__d_convergence| for SCF of post-HF,
   and 10 for E_CONVERGENCE for post-HF of post-HF.

.. _`sec:scfrec`:

Recommendations
~~~~~~~~~~~~~~~

The SCF code is quite flexible and powerful.
We have tried as much as possible to keep the number of options to
a minimum and to allow all options to be used in the presence of all other
options. Below are some rough words of advice about using the SCF code for
practical calculations:

* For |scf__guess|, the ``SAD`` guess is usually your
  friend, even for open-shell systems (at the very least, it gets the right
  number of electrons, unlike some other programs). For instance, we have found
  that a simple SAD guess is often as good as doing a full SCF in a 3-21G basis
  and then performing a cast-up, at a fraction of the cost.  However, SAD and
  DOCC/SOCC arrays do not play very well together at the moment.
* For wall time, ``DF`` may be a factor of ten or more faster than the exact
  integral technologies available in |PSIfour|.
  Use ``DF`` unless you need absolute accuracy or do not
  have a -JKFIT auxiliary set for your orbital basis/atom type. Then use
  ``DIRECT``.
* Don't mess with the DIIS convergence options unless convergence is a problem.
  We have optimized the parameters for efficiency over a wide array of system
  types.
* Buy a developer a beer!

The "best-practice" input file for HF is::

    memory 1 GB  # As much as you've got, the DF algorithm can use

    molecule {
    O
    H 1 1.0
    H 1 1.0 2 104.5
    }

    set {
    basis cc-pvdz
    scf_type df
    guess sad
    ints_tolerance 1.0E-10  # Even this is epically tight, 1.0E-8 is OK
    }

    energy('scf')
