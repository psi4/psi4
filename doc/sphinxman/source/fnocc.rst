.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2017 The Psi4 Developers.
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

.. index:: Frozen natural orbital coupled cluster, FNO-CC

.. index::
   single: QCISD(T)
   single: MP4
   single: CEPA
   single: LCCSD
   single: Frozen Natural Orbitals
   single: FNO-QCISD(T)
   single: FNO-MP4
   single: FNO-CCSD(T)
   single: DF-CCSD(T)
   pair: CCSD(T); density-fitting

.. _`sec:fnocc`:

FNOCC: Frozen natural orbitals for CCSD(T), QCISD(T), CEPA, and MP4
===================================================================

.. codeauthor:: A. Eugene DePrince
.. sectionauthor:: A. Eugene DePrince

*Module:* :ref:`Keywords <apdx:fnocc>`, :ref:`PSI Variables <apdx:fnocc_psivar>`, :source:`FNOCC <psi4/src/psi4/fnocc>`

.. warning:: There is a known bug concerning the i7-5930 series combined
   with the Intel 15 compilers and MKL 11.2.3. When |PsiFour| is compiled
   under these conditions, parallel runs of the FNOCC code have experienced
   nonsensical CCSD correlation energies (often several Hartrees lower
   than the starting guess). At the moment, the only confirmed solutions
   are running serially, using a different BLAS implementation, or upgrading
   to Intel 16.0.2 and MKL 11.3.2.

Frozen natural orbitals (FNO)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The computational cost of the CCSD [Purvis:1982]_, CCSD(T)
[Raghavachari:1989]_, and related methods be reduced by constructing a
compact representation of the virtual space based on the natural orbitals
of second-order perturbation theory [Sosa:1989:148]_.  The most demanding
steps in the CCSD and (T) algorithms scale as :math:`{\cal{O}}(o^2v^4)`
and :math:`{\cal{O}}(o^3v^4)`, where :math:`o` and :math:`v` represent the
number of oribitals that are occupied and unoccupied (virtual) in the
reference function, respectively.  By reducing the the size of the virtual
space, the cost of evaluating these terms reduces by a factor of :math:`(v
/ v_{FNO})^4`, where :math:`v_{FNO}` represents the number of virtual
orbitals retained after the FNO truncation.

The general outline for the FNO procedure in |Psifour| is:

    (i)   construct the virtual-virtual block of the unrelaxed MP2 one-particle density matrix (OPDM) 
    (ii)  diagonalize this block of the OPDM to obtain a set of natural virtual orbitals
    (iii) based on some occupancy threshold, determine which orbitals are unimportant and may be discarded
    (iv)  project the virtual-virtual block of the Fock matrix onto the truncated space
    (v)   construct semicanonical orbitals by diagonalizing the virtual-virtual block of the Fock matrix
    (vi)  proceed with the QCISD(T) / CCSD(T) / MP4 computation in the reduced virtual space

A second-order correction based upon the MP2 energies in the full and
truncated spaces captures much of the missing correlation effects.  More
details on the implementation and numerical accuracy of FNO methods in
|Psifour| can be found in [DePrince:2013:293]_\.  FNO computations
are controlled through the keywords |fnocc__nat_orbs| and
|fnocc__occ_tolerance|, or by prepending a valid method name with "fno" in
the energy call as ::

    energy('fno-ccsd(t)')

If you wish to specify the number of active natural orbitals manually, use
the keyword |fnocc__active_nat_orbs|.  This keyword will override the 
keyword |fnocc__occ_tolerance|.

QCISD(T), CCSD(T), MP4, and CEPA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The FNOCC module in |Psifour| supports several related many-body quantum
chemistry methods, including the CCSD(T) and QCISD(T) methods, several
orders of many-body perturbation theory (MP2-MP4), and a family methods
related to the coupled electron pair approximation (CEPA).

Quadratic configuration interaction and coupled cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The quadratic configuration interaction singles doubles (QCISD) method of
Pople, Head-Gordon, and Raghavachari [Pople:1987:5968]_\  was originally
presented as a size-consistent extension of configuration interaction
singles doubles (CISD). The method can also be obtained as a
simplified version of the coupled cluster singles doubles (CCSD)
method [Purvis:1982]_\.  Consider the set of equations defining CCSD:

.. math::
    :label: CCSD

    \langle \Psi_0  | (H - E) (1 + T_1 + T_2 + \frac{1}{2}T_1^2)|\Psi_0\rangle = 0, \\
    \langle \Psi_i^a  | (H - E) (1 + T_1 + T_2 + \frac{1}{2}T_1^2+T_1T_2+\frac{1}{3!}T_1^3)|\Psi_0\rangle = 0, \\
    \langle \Psi_{ij}^{ab}  | (H - E) (1 + T_1 + T_2 + \frac{1}{2}T_1^2 + T_1T_2+\frac{1}{3!}T_1^3+\frac{1}{2}T_2^2+\frac{1}{2}T_1^2T_2+\frac{1}{4!}T_1^4)|\Psi_0\rangle = 0, \\


where we have chosen the intermediate normalization, 
:math:`\langle \Psi_0| \Psi \rangle = 1`, and the symbols :math:`T_1` 
and :math:`T_2` represent single and double excitation operators.  The 
QCISD equations can be obtained by omitting all but two terms that 
are nonlinear in :math:`T_1` and :math:`T_2`:

.. math::
    :label: QCISD

    \langle \Psi_0  | (H - E) (1 + T_1 + T_2)|\Psi_0\rangle = 0, \\
    \langle \Psi_i^a  | (H - E) (1 + T_1 + T_2 + T_1T_2)|\Psi_0\rangle = 0, \\
    \langle \Psi_{ij}^{ab}  | (H - E) (1 + T_1 + T_2 + \frac{1}{2}T_2^2)|\Psi_0\rangle = 0. \\

QCISD is slightly cheaper that CCSD computationally, but it retains the
:math:`{\cal{O}}(o^2v^4)` complexity of the original equations. Just as in
the familiar CCSD(T) method, the effects of connected triple excitations
may be included noniteratively to yield the QCISD(T) method.  Both the
QCISD(T) and CCSD(T) methods are implemented for closed-shell references
in the FNOCC module.

.. _`sec:fnompn`:

Many-body perturbation theory 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

QCI and CC methods are closely related to perturbation theory, and the
MP2, MP3, and MP4(SDQ) correlation energies can be obtained as a free
by-product of a CCSD or QCISD computation.  The following is an 
example of the results for a computation run with the call
``energy('fno-qcisd')`` to :py:func:`~psi4.energy`::

  QCISD iterations converged!

        OS MP2 FNO correction:                -0.000819116338
        SS MP2 FNO correction:                -0.000092278158
        MP2 FNO correction:                   -0.000911394496

        OS MP2 correlation energy:            -0.166478414245
        SS MP2 correlation energy:            -0.056669079827
        MP2 correlation energy:               -0.223147494072
      * MP2 total energy:                    -76.258836941658

        OS MP2.5 correlation energy:          -0.171225850256
        SS MP2.5 correlation energy:          -0.054028401038
        MP2.5 correlation energy:             -0.225254251294
      * MP2.5 total energy:                  -76.260943698880

        OS MP3 correlation energy:            -0.175973286267
        SS MP3 correlation energy:            -0.051387722248
        MP3 correlation energy:               -0.227361008515
      * MP3 total energy:                    -76.263050456101

        OS MP4(SDQ) correlation energy:       -0.180324322304
        SS MP4(SDQ) correlation energy:       -0.048798468084
        MP4(SDQ) correlation energy:          -0.230995119324
      * MP4(SDQ) total energy:               -76.266684566910

        OS QCISD correlation energy:          -0.181578117924
        SS QCISD correlation energy:          -0.049853548145
        QCISD correlation energy:             -0.231431666069
      * QCISD total energy:                  -76.267121113654

The first set of energies printed corresponds to the second-order FNO 
correction mentioned previously.  Results for many-body perturbation 
theory through partial fourth order are then provided.
The notation MP4(SDQ) indicates that we have included all contributions to
the correlation energy through fourth order, with the exception of that
from connected triple excitations.  

One need not run a full QCISD or CCSD computation to obtain these
perturbation theory results.  The keywords for invoking perturbation
theory computations are given below in
Table :ref:`FNOCC Methods <table:fnocc_methods>`.  Full MP4 correlation
energies are also available.

.. _`sec:fnocepa`:

Coupled electron pair approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Coupled-pair methods can be viewed as approximations to CCSD or as
size-extensive modifications of CISD.  The methods have the same
complexity as CISD, and solving the CISD or coupled-pair equations
requires fewer floating point operations than solving the CCSD.  CISD,
CCSD, and the coupled-pair methods discussed below all scale formally with
the sixth power of system size, and, as with the QCISD method, CEPA
methods retain :math:`{\cal{O}}(o^2v^4)` complexity of the CCSD equations.
For a detailed discussion of the properties of various coupled-pair
methods, see [Wennmohs:2008:217]_\.

What follows is a very basic description of the practical differences in
the equations that define each of the coupled-pair methods implemented in
|Psifour|.  We begin with the CISD wave function

.. math::
    :label: CIwfn

    | \Psi \rangle = | \Psi_0 \rangle + \sum_i^{occ} \sum_a^{vir} t_i^a | \Psi_i^a\rangle + \frac{1}{4}\sum_{ij}^{occ} \sum_{ab}^{vir} t_{ij}^{ab} | \Psi_{ij}^{ab}\rangle,

where we have chosen the intermediate normalization, :math:`\langle \Psi_0
| \Psi \rangle = 1`.  The CISD correlation energy is given by

.. math::
    :label: CIenergy
    
    E_c = \langle \Psi_0 | \hat{H} - E_0 | \Psi \rangle,

and the amplitudes can be determined by the solution to the coupled set of
equations:

.. math::
    :label: CIeqns
    
    0   &= \langle \Psi_{ij}^{ab} | \hat{H} - E_0 - E_c | \Psi \rangle, \\
    0   &= \langle \Psi_{i}^{a} | \hat{H} - E_0 - E_c | \Psi \rangle.

The CISD method is not size-extensive, but this problem can be overcome by
making very simple modifications to the amplitude equations.  We replace
the correlation energy, :math:`E_c`, with generalized shifts for the
doubles and singles equations, :math:`\Delta_{ij}` and :math:`\Delta_i`:

.. math::
    :label: CEPAeqns
    
    0   &= \langle \Psi_{ij}^{ab} | \hat{H} - E_0 - \Delta_{ij} | \Psi \rangle, \\
    0   &= \langle \Psi_{i}^{a} | \hat{H} - E_0 - \Delta_i | \Psi \rangle.

These shifts approximate the effects of triple and quadruple excitations.
The values for :math:`\Delta_{ij}` and :math:`\Delta_i`  used in several
coupled-pair methods are given in Table :ref:`CEPA Shifts
<table:cepa_shifts>`.  Note that these shifts are defined in a spin-free
formalism for closed-shell references only.

    .. _`table:cepa_shifts`:

    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | method                  | :math:`\Delta_{ij}`                                        |  :math:`\Delta_i`                            |
    +=========================+============================================================+==============================================+
    | cisd                    | :math:`E_c`                                                |  :math:`E_c`                                 |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | cepa(0)                 | 0                                                          |  0                                           |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | cepa(1)                 | :math:`\frac{1}{2}\sum_k(\epsilon_{ik}+\epsilon_{jk})`     | :math:`\sum_k \epsilon_{ik}`                 |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | cepa(3)                 | :math:`-\epsilon_{ij}+\sum_k(\epsilon_{ik}+\epsilon_{jk})` | :math:`-\epsilon_{ii}+2\sum_k \epsilon_{ik}` |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | acpf                    | :math:`\frac{2}{N} E_c`                                    | :math:`\frac{2}{N} E_c`                      |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | aqcc                    | :math:`[1-\frac{(N-3)(N-2)}{N(N-1)}]E_c`                   | :math:`[1-\frac{(N-3)(N-2)}{N(N-1)}]E_c`     |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+

.. comment    | dci                     | :math:`E_c`                                                |  NA                                          |
.. comment    +-------------------------+------------------------------------------------------------+----------------------------------------------+

The pair correlation energy, :math:`\epsilon_{ij}`, is simply a partial
sum of the correlation energy.  In a spin-free formalism, the pair energy
is given by

.. math::
   :label: pair_energy

   \epsilon_{ij} = \sum_{ab} v_{ij}^{ab} (2 t_{ij}^{ab} - t_{ij}^{ba})

Methods whose shifts (:math:`\Delta_{ij}` and :math:`\Delta_i`) do not
explicitly depend on orbitals :math:`i` or :math:`j` (CISD, CEPA(0), ACPF,
and AQCC) have solutions that render the energy stationary with respect
variations in the amplitudes.  This convenient property allows density
matrices and 1-electron properties to be evaluated without any additional
effort.  Note, however, that 1-electron properties are currently
unavailable when coupling these stationary CEPA-like methods with frozen
natural orbitals.

Density-fitted coupled cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Density fitting (DF) [or the resolution of the identity (RI)] and Cholesky
decomposition (CD) techniques are popular in quantum chemistry to avoid
the computation and storage of the 4-index electron repulsion integral
(ERI) tensor and even to reduce the computational scaling of some terms.
DF/CD-CCSD(T) computations are available in |Psifour|, with or without the
use of FNOs, through the FNOCC module.  The implementation and accuracy of
the DF/CD-CCSD(T) method are described in [DePrince:2013:2687]_\.

The DF-CCSD(T) procedure uses two auxiliary basis sets.  The first set is
that used in the SCF procedure, defined by the |scf__df_basis_scf|
keyword.  If this keyword is not specified, an appropriate -JKFIT set is
automatically selected.  This auxiliary set defines the ERIs used to
build the Fock matrix used in the DF-CCSD(T) procedure.  The second
auxiliary set is used to approximate all other ERIs in the DF-CCSD(T)
procedure. The choice of auxiliary basis is controlled by the keyword
|fnocc__df_basis_cc|.  By default, |fnocc__df_basis_cc| is the RI set
(optimized for DFMP2) most appropriate for use with the primary basis.
For example, if the primary basis is aug-cc-pVDZ, the default
|fnocc__df_basis_cc| will be aug-cc-pVDZ-RI.

Alternatively, the user can request that the DF-CCSD(T) procedure use a
set of vectors defined by the Cholesky decomposition of the ERI tensor as
the auxiliary basis. This feature is enabled by specifying |globals__cc_type| ``CD``.
CD methods can be enabled in the SCF
procedure as well, by specifying the |scf__scf_type| as ``CD``.  The
accuracy of the decomposition can be controlled through the keyword
|fnocc__cholesky_tolerance|.

.. comment This feature is enabled by specifying |fnocc__df_basis_cc| as "CHOLESKY".  

The following input file sets up a DF-CCSD(T)
computation using CD integrals ::

    molecule h2o {
        0 1
        O
        H 1 1.0
        H 1 1.0 2 104.5
    }
    
    set {
        scf_type cd
        cc_type cd
        basis aug-cc-pvdz
        freeze_core true
    }
    energy('ccsd(t)')

The resulting CCSD(T) correlation energy will be equivalent to that
obtained from a conventional computation if |fnocc__cholesky_tolerance| is
sufficiently small (*e.g.* ``1e-9``).

.. _`sec:fnogn`:

Gn theory
~~~~~~~~~

The FNOCC module contains all the components that comprise the Gn family
of composite methods.  Currently, only the G2 method is supported
[Curtiss:1991:7221]_\.  The G2 procedure may be called through the
:py:func:`~psi4.energy` wrapper: ::

    energy('gaussian-2')

Supported methods
~~~~~~~~~~~~~~~~~

The various methods supported by the FNOCC module in |Psifour| are detailed
in Table :ref:`FNOCC Methods <table:fnocc_methods>`.  Note that these methods
are implemented for closed-shell references only.  For open-shell references,
the calls ``energy('mp2.5')``, ``energy('mp3')``, and ``energy('mp4')`` will
default to implementations of these methods in :ref:`other modules <table:managedmethods>`.

    .. _`table:fnocc_methods`:

    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | name                    | calls method                                                | type select                               |
    +=========================+=============================================================+===========================================+
    | qcisd                   | quadratic configuration interaction singles doubles         | |globals__ci_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | qcisd(t)                | qcisd with perturbative triples                             | |globals__ci_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | mp2.5                   | average of second- and third-order perturbation theories    | |globals__mp_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | mp3                     | third-order perturbation theory                             | |globals__mp_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | mp4(sdq)                | fourth-order perturbation theory, minus triples contribution| |globals__mp_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | mp4                     | full fourth-order perturbation theory                       | |globals__mp_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | lccd                    | linear ccd                                                  | |globals__cc_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | cepa(0), lccsd          | coupled electron pair approximation, variant 0              | |globals__cc_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | cepa(1)                 | coupled electron pair approximation, variant 1              | |globals__cc_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | cepa(3)                 | coupled electron pair approximation, variant 3              | |globals__cc_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | acpf                    | averaged coupled-pair functional                            | |globals__cc_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | aqcc                    | averaged quadratic coupled-cluster                          | |globals__cc_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | cisd                    | configuration interaction with single and double excitations| |globals__ci_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-qcisd               | qcisd with frozen natural orbitals                          | |globals__ci_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-qcisd(t)            | qcisd(t) with frozen natural orbitals                       | |globals__ci_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-ccsd                | coupled cluster singles doubles with frozen natural orbitals| |globals__cc_type| CONV, DF, CD           |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-ccsd(t)             | ccsd with perturbative triples and frozen natural orbitals  | |globals__cc_type| CONV, DF, CD           |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-mp3                 | mp3 with frozen natural orbitals                            | |globals__mp_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-mp4(sdq)            | mp4(sdq) with frozen natural orbitals                       | |globals__mp_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-mp4                 | mp4 with frozen natural orbitals                            | |globals__mp_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-lccd                | linear ccd with frozen natural orbitals                     | |globals__cc_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-cepa(0), fno-lccsd  | cepa(0) with frozen natural orbitals                        | |globals__cc_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-cepa(1)             | cepa(1) with frozen natural orbitals                        | |globals__cc_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-cepa(3)             | cepa(3) with frozen natural orbitals                        | |globals__cc_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-acpf                | acpf with frozen natural orbitals                           | |globals__cc_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-aqcc                | aqcc with frozen natural orbitals                           | |globals__cc_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+
    | fno-cisd                | cisd with frozen natural orbitals                           | |globals__ci_type| CONV                   |
    +-------------------------+-------------------------------------------------------------+-------------------------------------------+

.. comment    | df-ccsd                 | ccsd with density fitting                                   |
.. comment    +-------------------------+-------------------------------------------------------------+
.. comment    | df-ccsd(t)              | ccsd(t) with density fitting                                |
.. comment    +-------------------------+-------------------------------------------------------------+
.. comment    | fno-df-ccsd             | ccsd with density fitting and frozen natural orbitals       |
.. comment    +-------------------------+-------------------------------------------------------------+
.. comment    | fno-df-ccsd(t)          | ccsd(t) with density fitting and frozen natural orbitals    |
.. comment    +-------------------------+-------------------------------------------------------------+
.. comment    | dci                     | configuration interaction with double excitations           |
.. comment    +-------------------------+-------------------------------------------------------------+
.. comment    | fno-dci                 | dci with frozen natural orbitals                            |
.. comment    +-------------------------+-------------------------------------------------------------+

.. index:: FNOCC; basic-keywords

Basic FNOCC Keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: /autodir_options_c/mints__basis.rst
.. include:: /autodir_options_c/globals__freeze_core.rst
.. include:: /autodir_options_c/fnocc__r_convergence.rst
.. include:: /autodir_options_c/fnocc__e_convergence.rst
.. include:: /autodir_options_c/fnocc__maxiter.rst
.. include:: /autodir_options_c/fnocc__diis_max_vecs.rst
.. include:: /autodir_options_c/fnocc__nat_orbs.rst
.. include:: /autodir_options_c/fnocc__occ_tolerance.rst
.. include:: /autodir_options_c/fnocc__triples_low_memory.rst
.. include:: /autodir_options_c/fnocc__cc_timings.rst
.. include:: /autodir_options_c/fnocc__df_basis_cc.rst
.. include:: /autodir_options_c/fnocc__cholesky_tolerance.rst
.. include:: /autodir_options_c/fnocc__cepa_no_singles.rst
.. include:: /autodir_options_c/fnocc__dipmom.rst

.. index:: FNOCC; advanced-keywords

Advanced FNOCC Keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. include:: /autodir_options_c/fnocc__scs_mp2.rst
.. include:: /autodir_options_c/fnocc__mp2_scale_os.rst
.. include:: /autodir_options_c/fnocc__mp2_scale_ss.rst
.. include:: /autodir_options_c/fnocc__scs_ccsd.rst
.. include:: /autodir_options_c/fnocc__cc_scale_os.rst
.. include:: /autodir_options_c/fnocc__cc_scale_ss.rst
.. include:: /autodir_options_c/fnocc__run_mp2.rst
.. include:: /autodir_options_c/fnocc__run_mp3.rst
.. include:: /autodir_options_c/fnocc__run_mp4.rst
.. include:: /autodir_options_c/fnocc__run_ccsd.rst
.. include:: /autodir_options_c/fnocc__run_cepa.rst
.. include:: /autodir_options_c/fnocc__compute_triples.rst
.. include:: /autodir_options_c/fnocc__compute_mp4_triples.rst
.. include:: /autodir_options_c/fnocc__dfcc.rst
.. include:: /autodir_options_c/fnocc__cepa_level.rst
