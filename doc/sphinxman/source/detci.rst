.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2016 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This program is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU General Public License as published by
.. # the Free Software Foundation; either version 2 of the License, or
.. # (at your option) any later version.
.. #
.. # This program is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU General Public License for more details.
.. #
.. # You should have received a copy of the GNU General Public License along
.. # with this program; if not, write to the Free Software Foundation, Inc.,
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
|PSIfour| is not optimized for CISD computations.  Instead, emphasis
has been placed on developing a very efficient program to handle more
general CI wavefunctions which may be helpful in more challenging cases
such as highly strained molecules or bond breaking reactions.  The CI
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
satisfying these rules are included in the CI.

The DETCI module is also very efficient at computing full configuration
interaction
wavefunctions, and it is used in this capacity in the complete-active-space
self-consistent-field (CASSCF) code.  Use of DETCI for CASSCF
wavefunctions is described in another section of this manual.

As just mentioned, the DETCI module is designed for challenging
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

The division of the molecular orbitals into various subspaces such as
RAS spaces, or frozen *vs.* active orbitals, *etc.*, needs to be clear not
only to detci, but also at least to the transformation program
(and in the case of MCSCF, to other programs as well).  Thus, orbital
subspace keywords such as |globals__ras1|,
|globals__ras2|, |globals__ras3|, |globals__frozen_docc|, 
|globals__frozen_uocc|,
|globals__active|, etc., are set
in the global section of input so they may also be read by other modules.

For single-reference CI computations, the easiest way to invoke a CI
computation with DETCI is simply to call :py:func:`~psi4.energy`, :py:func:`~psi4.optimize`, *etc.*,
with the common name for that CI wavefunction, like ``energy('cisd')`` 
for a CISD single-point energy.  The Python driver
recognizes ``cisd``, ``cisdt``, and ``cisdtq``.  Higher order
single-refernce CI wavefunctions, like those including singles through
6-fold excitations, can be invoked using numbers, like ``ci6``.  A full
CI can be specifed by ``fci``.  More complicated CI computations, like
RASCI, can be performed by setting the appropriate keywords and calling the
module generically like ``energy('detci')``.  The latter approach
will also work for any of the previously-mentioned CI wavefunctions for
which the driver has built-in shortcuts, so long as the relevant options
(especially |detci__ex_level|) are set appropriately.  Some
examples of single-refence CI, RASCI, and full CI computations are provided
in :source:`samples`.

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
   pair: CI; multi-configurational self-consistent field

Multi-Configurational Self-Consistent Field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As the cost of Full CI scales exponentially with respect to the number of
active orbitals it is often advantageous to neglect orbitals that do not
exhibit strong correlation. These orbitals are variationally optimized
simultaneously with the CI coefficients and known as Multi-Configurational
Self-Consistent Field (MCSCF). The most commonly used MCSCF procedure is the
complete-active-space self-consistent-field (CASSCF) approach [Roos:1980]_,
which includes all possible determinants (with the proper symmetry) that can be
formed by distributing a set of active electrons among a set of active
orbitals. The MCSCF module performs CASSCF optimization of molecular orbitals
via a two-step procedure in which the CI coefficients and orbitals are
optimized in an alternating manner. The program uses a fairly simple
approximate orbital Hessian [Chaban:1997:88]_ and a Newton-Raphson update,
accelerated by Pulay's DIIS procedure [Pulay:1980]_. We have also implemented
the RASSCF method [Malmqvist:1990:RASSCF]_, which is another kind of MCSCF
which is typically less complete (and less expensive) than CASSCF.

Inactive orbitals in the MCSCF may be specified by the
|globals__restricted_docc| and |globals__restricted_uocc| keywords. These
orbitals will remain doubly-occupied or doubly-unoccupied, respectively, in the
MCSCF wavefunction.  However, the form of these orbitals will be optimized in
the MCSCF procedure.  It is also possible to literally freeze inactive orbitals
in their original (SCF) form using the |globals__frozen_docc| and
|globals__frozen_uocc| keywords.  This is not normally what one wishes to do in
an MCSCF computation (*e.g.*, it complicates the computation of gradients), but
it can make the computations faster and is helpful in some circumstances where
unphysical mixing of inactive and active occupied orbitals might occur.
Presently, it is not possible to mix the use of restricted and frozen orbitals
in |PSIfour|.

An illustrative CASSCF example is as follows::

    molecule {
    O
    H 1 1.00
    H 1 1.00 2 103.1
    }
    
    set {
        basis           6-31G**
        restricted_docc [1, 0, 0, 0]
        active          [3, 0, 1, 2]
    }
    energy('casscf')

This input will compute the CASSCF energy of water where the 1s Oxygen orbital
and several virtual orbitals are not included in the CI expansion, but are
still optimizd. The following is a full list of spaces within the various MCSCF
types.

.. _`table:mcscf_spaces`:

.. table:: Spaces available in the MCSCF program

    +----------------------------+----------------------------+----------------------------+
    | CI                         | RASSCF                     | CASSCF                     |
    +============================+============================+============================+
    |                            | |globals__frozen_uocc|     | |globals__frozen_uocc|     |
    +                            +----------------------------+----------------------------+
    | |globals__restricted_uocc| | |globals__restricted_uocc| | |globals__restricted_uocc| |
    +----------------------------+----------------------------+----------------------------+
    | |globals__active|          | |globals__ras4|            | |globals__active|          |
    +                            +----------------------------+                            +
    |                            | |globals__ras3|            |                            |
    +                            +----------------------------+                            +
    |                            | |globals__ras2|            |                            |
    +                            +----------------------------+                            +
    |                            | |globals__ras1|            |                            |
    +----------------------------+----------------------------+----------------------------+
    | |globals__restricted_docc| | |globals__restricted_docc| | |globals__restricted_docc| |
    +                            +----------------------------+----------------------------+
    |                            | |globals__frozen_docc|     | |globals__frozen_docc|     |
    +----------------------------+----------------------------+----------------------------+

Basic MCSCF Keywords
~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/detci__mcscf_e_convergence.rst
.. include:: autodir_options_c/detci__mcscf_r_convergence.rst
.. include:: autodir_options_c/detci__mcscf_type.rst
.. include:: autodir_options_c/detci__mcscf_algorithm.rst
.. include:: autodir_options_c/detci__mcscf_maxiter.rst
.. include:: autodir_options_c/detci__mcscf_rotate.rst
.. include:: autodir_options_c/detci__mcscf_diis_start.rst


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

