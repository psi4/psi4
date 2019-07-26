.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2019 The Psi4 Developers.
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

.. index:: CI

.. index::
   pair: CI; theory

.. _`sec:ci`:

CI: Configuration Interaction
=============================

.. codeauthor:: Daniel G. A. Smith, C. David Sherrill, and Matthew L. Leininger
.. sectionauthor:: Daniel G. A. Smith and C. David Sherrill

*Module:* :ref:`Keywords <apdx:detci>`, :ref:`PSI Variables <apdx:detci_psivar>`, :source:`DETCI <psi4/src/psi4/detci>`

Configuration interaction (CI) is one of the most general ways to
improve upon Hartree--Fock theory by adding a description of the
correlated motions of electrons.  Simply put, a CI wavefunction
is a linear combination of Slater determinants (or spin-adapted
configuration state functions), with the linear coefficients being
determined variationally via diagonalization of the Hamiltonian in the
given subspace of determinants.  For a "single-reference" CI based
on reference function :math:`| \Phi_0 \rangle`, we can write the CI expansion as
follows:

.. math:: | \Psi \rangle = c_0 | \Phi_0 \rangle
   + \sum_i^{\rm occ} \sum_a^{\rm vir} c_i^a | \Phi_i^a \rangle
   + \sum_{i<j}^{\rm occ} \sum_{a<b}^{\rm vir} c_{ij}^{ab} 
   | \Phi_{ij}^{ab} \rangle
   + \sum_{i<j<k}^{\rm occ} \sum_{a<b<c}^{\rm vir} c_{ijk}^{abc}
   | \Phi_{ijk}^{abc} \rangle + \cdots
   :label: CIexpansion

The simplest standard CI method that improves upon Hartree--Fock is a CI
that adds all singly :math:`| \Phi_i^a \rangle` and doubly 
:math:`| \Phi_{ij}^{ab} \rangle`
substituted determinants (CISD) to the reference determinant
:math:`| \Phi_0 \rangle`.  The CISD wavefunction has fallen out of favor
because truncated CI wavefunctions are not size-extensive, meaning
that their quality degrades for larger molecules.  MP2 is a less
expensive alternative giving results similar to those of CISD for small
molecules, but the quality of MP2 does not degrade for larger molecules.
Coupled-cluster singles and doubles (CCSD) is another size-extensive
alternative; it is only slightly more costly computationally than CISD,
but it typically provides significantly more accurate results.

The CI code in |PSIfour| is described in detail in 
[Sherrill:1999:CI]_.  For the reasons stated above, the CI code in
|PSIfour| is not optimized for CISD computations, and it uses data structures
that are particularly inefficient for CISD and may result in the program
running out of memory and crashing for CISD except on very small molecules.
Instead, DETCI was designed to be efficient
in handling more highly correlated CI wavefunctions that can be helpful in more 
challenging cases such as highly strained molecules or bond breaking reactions.  The CI
code is based on the fast, determinant-based string formalism
of Handy [Handy:1980]_.  It can solve for restricted active space
configuration interaction (RAS CI) wavefunctions as described by Olsen,
Roos, Jorgensen, and Aa. Jensen [Olsen:1988]_.  Excitation-class
selected multi-reference CI wavefunctions, such as second-order CI,
can be formulated as RAS CI's.  A RAS CI selects determinants for the
model space as those which have no more than :math:`n` holes in the lowest set
of orbitals (called RAS I) and no more than :math:`m` electrons in the highest
set of orbitals (called RAS III).  An intermediate set of orbitals, if
present (RAS II), has no restrictions placed upon it.  All determinants
satisfying these rules are included in the RAS CI.

The DETCI module is also very efficient at computing full configuration
interaction
wavefunctions, and it is used in this capacity in the complete-active-space
self-consistent-field (CASSCF) code.  It can also perform approximate
CASSCF computations in which one uses RAS restrictions on the CI excitations,
rather than doing a full CI in the active space.  This is called a 
RASSCF.  CASSCF and RASSCF computations are types of multi-configurational
self-consistent-field procedures, and are described in :ref:`sec:mcscf`.

As mentioned above, the DETCI module is designed for challenging
chemical systems for which simple CISD is not suitable.  Because
CI wavefunctions which go beyond CISD (such as RAS CI) are fairly complex,
typically the DETCI code will be used in cases where the
tradeoffs between computational expense and completeness of the
model space are nontrivial.  Hence, the user is advised to develop
a good working knowledge of multi-reference and RAS CI methods before
attempting to use the program for a production-level project.  This user's
manual will provide only an elementary introduction to the most
important keywords.  Additional information is available in the complete
list of keywords for DETCI provided in Appendix :ref:`apdx:detci`.

For single-reference CI computations, the easiest way to invoke a CI
computation with DETCI is simply to call :py:func:`~psi4.energy`, :py:func:`~psi4.optimize`, *etc.*,
with the common name for that CI wavefunction, like ``energy('cisd')`` 
for a CISD single-point energy.  The Python driver
recognizes ``cisd``, ``cisdt``, and ``cisdtq``.  Higher order
single-reference CI wavefunctions, like those including singles through
6-fold excitations, can be invoked using numbers, like ``ci6``.  A full
CI can be specified by ``fci``.  More complicated CI computations, like
RASCI, can be performed by setting the appropriate keywords and calling the
module generically like ``energy('detci')``.  The latter approach
will also work for any of the previously-mentioned CI wavefunctions for
which the driver has built-in shortcuts, so long as the relevant options
(especially |detci__ex_level|) are set appropriately.  Some
examples of single-refence CI, RASCI, and full CI computations are provided
in :source:`samples`.

.. _`table:ci_spaces`:

.. table:: Orbital spaces for CI computations

    +----------------------------+----------------------------+-------------------------------+
    | CI (e.g., CISD, FCI)       | RASCI                      | CASCI                         |
    +============================+============================+===============================+
    | |globals__frozen_uocc|     | |globals__frozen_uocc|     | |globals__frozen_uocc| [#f1]_ |
    +----------------------------+----------------------------+-------------------------------+
    | (all orbitals not in       | |globals__ras4|            |                               |
    + |globals__frozen_uocc|     +----------------------------+                               +
    | or |globals__frozen_docc|  | |globals__ras3|            |                               |
    + are included in CI)        +----------------------------+                               +
    |                            | |globals__ras2|            |                               |
    +                            +----------------------------+                               +
    |                            | |globals__ras1|            | |globals__active|             |
    +----------------------------+----------------------------+-------------------------------+
    | |globals__frozen_docc|     | |globals__frozen_docc|     | |globals__frozen_docc|        |
    +----------------------------+----------------------------+-------------------------------+

.. [#f1] |globals__frozen_uocc| is not required and will be inferred if 
   |globals__active| is provided.  However, if it is easier to specify
   |globals__frozen_uocc|, then this may be provided and |globals__active| can
   be inferred.

The table above shows the relevant orbitals spaces for CI computations (an
analogous :ref:`table <table:mcscf_spaces>` for MCSCF is also available).  
The third column of the
table refers to CASCI, in which a full CI is performed in some smaller
set of ``active`` orbitals; it is equivalent to CASSCF except without
any orbital optimization.  It can be invoked via ``energy('fci')``
with appropriate values selected for |globals__frozen_docc| and
|globals__active|.  For CI computations, there is no difference between
|globals__frozen_docc| and |globals__restricted_docc|, or between
|globals__frozen_uocc| and |globals__restricted_uocc|.  There are
differences between these keywords for :ref:`sec:mcscf`.

.. index:: CI; basic-keywords

Basic DETCI Keywords
~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/detci__reference.rst
.. include:: autodir_options_c/detci__r_convergence.rst
.. include:: autodir_options_c/detci__ex_level.rst
.. include:: autodir_options_c/detci__fci.rst
.. include:: autodir_options_c/globals__frozen_docc.rst
.. include:: autodir_options_c/globals__restricted_docc.rst
.. include:: autodir_options_c/globals__restricted_uocc.rst
.. include:: autodir_options_c/globals__frozen_uocc.rst
.. include:: autodir_options_c/detci__ci_maxiter.rst
.. include:: autodir_options_c/detci__num_roots.rst
.. include:: autodir_options_c/detci__icore.rst
.. include:: autodir_options_c/detci__diag_method.rst
.. include:: autodir_options_c/detci__opdm.rst
.. include:: autodir_options_c/detci__tdm.rst
.. include:: autodir_options_c/detci__dipmom.rst
.. include:: autodir_options_c/detci__mpn.rst

For larger computations, additional keywords may be required, as
described in the DETCI section of the Appendix :ref:`apdx:detci`.

.. index:: 
   pair: CI; arbitrary-order perturbation theory

.. _`sec:arbpt`:

Arbitrary Order Perturbation Theory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DETCI module is capable of computing energies for arbitrary
order |MollerPlesset| perturbation theory (MPn, for closed-shell
systems with an RHF reference) and for Z-averaged perturbation theory
(ZAPTn, open-shell systems with a ROHF reference). However, please
note that these computations are essentially doing high-order CI (up to
full CI) computations to obtain these results, and hence they will only
be possible for very small systems (generally a dozen electrons or less).

The simplest way to run high-order perturbation theory computations is to
call, *e.g.*, ``energy('mp10')`` to invoke a MP10 computation or
``energy('zapt25')`` to invoke a ZAPT25 computation.  This will
automatically set several additional user options to their appropriate
values.  The program uses the Wigner (2n+1) rule to obtain higher-order
energies from lower-order wavefunctions.

For the interested reader, the additional user options that are
automatically set up by the calls above are as follows.  A call like
``energy('mp10')`` sets |detci__mpn| to TRUE.
The program uses the Wigner (2n+1) rule by default
(|detci__mpn_wigner| = TRUE)
and figures out what order of wavefunction is
necessary to reach the desired order in the energy.  The program then
sets |detci__max_num_vecs| to the required order in the
wavefunction.
By default, the requested n-th order energy is saved as the current
energy to the process environment.   ZAPTN works essentially the same
way for an ROHF reference.

.. index:: 
   pair: CI; arbitrary-order coupled-cluster theory

Arbitrary Order Coupled-Cluster Theory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*This DETCI-based version of this feature is not yet released.  However,
the current version of the code does include an interface to*
:ref:`Kallay's MRCC <sec:mrcc>` *code.*

The DETCI module is also capable of computing arbitrary-order
coupled-cluster energies, using an approach similar to that of Hirata
and Bartlett [Hirata:2000:216]_, or of Olsen [Olsen:2000:7140]_.
Notably, the approach in DETCI also allows arbitrary-order 
*active space* coupled-cluster procedures.  The general algorithm
for doing this in DETCI is inefficient compared to optimized
lower-order coupled-cluster codes and should not be used for CCSD,
where the CCENERGY module is much more efficient.  For higher-order
CC (like CCSDT and beyond), the code is also not as efficient as the
MRCC code by K\ |a_acute|\ llay, to which |PSIfour| can interface (see Section
:ref:`sec:mrcc`); however, it may allow certain truncations of the model
space that might not be available presently in MRCC.  For very small
systems, the code can be useful for testing of, for example, CCSDTQ or
its active-space CCSDtq analog [Piecuch:1999:6103]_.

To perform arbitrary-order coupled-cluster, set the DETCI
option |detci__cc| to TRUE, and set
|detci__cc_ex_level| (note: not |detci__ex_level|)
to the desired coupled-cluster excitation level, and invoke 
``energy('detci')``.  Various other DETCI options have a similar
option for coupled-cluster, usually named beginning with CC.  The full
list of options is given in Appendix :ref:`apdx:detci`.

