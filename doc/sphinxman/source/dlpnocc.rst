.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2026 The Psi4 Developers.
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

.. index::
   single: DLPNO-CC
   single: DLPNO-CCSDT
   single: DLPNO-CCSDT(Q)
   single: DLPNO-CCSDTQ

.. _`sec:dlpnocc`:

DLPNO Coupled-Cluster Methods Through CCSDTQ
============================================

.. codeauthor:: Andy Jiang 
.. sectionauthor:: Andy Jiang

*Module:* :ref:`Keywords <apdx:dlpno>`, :ref:`PSI Variables <apdx:dlpno_psivar>`, :source:`DLPNO coupled-cluster methods <psi4/src/psi4/dlpno>`

Introduction
------------

The CCSD(T) method is considered the "gold standard" method of quantum chemistry,
often yielding chemically accurate (< 1 kcal/mol) results relative to experiment or
full configuration interaction (FCI) computations. Unfortunately, the high 
:math:`\mathcal{O}(N^{7})` scaling of the canonical gold-standard CCSD(T) method
limits its applicability for larger molecules.

The domain-based local pair natural orbital (DLPNO)-CCSD(T) approach of Neese and coworkers 
[Riplinger:2013:034106]_ [Riplinger:2013:134101]_ [Riplinger:2016:024109]_ 
[Guo:2018:011101]_ uses localized orbitals to exploit molecular sparsity for a linear-scaling approximation 
to canonical CCSD(T). These errors are controllable through user-set parameters. Computations on large molecules, 
such as crambin (636 atoms) [Riplinger:2013:134101]_ and insulin [Jiang:2024:082502]_ have been performed
with this algorithm. This algorithm was originally implemented in the ORCA software package, but it
is now available in |PSIfour|. Psi4's in-tree implementation provides DLPNO-CCSD, DLPNO-CCSD(T0),
DLPNO-CCSD(T), DLPNO-CCSDT, DLPNO-CCSDT(Q0), DLPNO-CCSDT(Q), and DLPNO-CCSDTQ energies.

For a more comprehensive overview on local correlation and the DLPNO-CCSD(T) algorithm, the reader is referred
to the :ref:`DLPNO-MP2 <sec:dlpnomp2>` documentation, and the published work of Jiang et al. describing the
implementation of DLPNO-CCSD(T) in |Psifour| [Jiang:2024:082502]_.

An example input file for a DLPNO-CCSD computation is::
   
   memory 2 GB

   molecule h2o {
   0 1
   O
   H 1 1.0
   H 1 1.0 2 104.5
   symmetry c1
   }
   
   set basis cc-pvdz
   set scf_type df
   set freeze_core True
   set pno_convergence tight
   
   energy('dlpno-ccsd')

An example input file for a DLPNO-CCSD(T) computation is::

   memory 2 GB

   molecule h2o {
   0 1
   O
   H 1 1.0
   H 1 1.0 2 104.5
   symmetry c1
   }
   
   set basis cc-pvdz
   set scf_type df
   set freeze_core True
   set pno_convergence tight
   
   energy('dlpno-ccsd(t)') # dlpno-ccsd(t0) for the semicanonical (T0) computation

Higher-Order DLPNO Methods
--------------------------

Recent work [Jiang:2025:2386]_ [Jiang:2025:144102]_ [Jiang:2026:2825]_ has extended the
DLPNO framework to coupled-cluster methods beyond CCSD(T). The same locality principles
used for DLPNO-CCSD(T) are applied to full triples, perturbative quadruples, and full
quadruples, with local natural-orbital spaces constructed for each significant occupied
orbital tuple.

.. note::

   For DLPNO-CCSDT and higher methods, |PSIfour| sets
   |dlpno__t_cut_pairs_mp2| equal to |dlpno__t_cut_pairs|. This collapses the
   separate weak-pair interval used by the lower-order DLPNO methods: pairs below
   the common threshold are eliminated, while every retained pair is treated as
   strong. Triplets and quadruplets in these post-CCSD(T) methods are therefore
   constructed only from retained strong pairs, and quadruplet construction has no
   separate maximum-weak-pairs option.

Strong and Weak Pairs Versus Occupied Tuples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The words *strong* and *weak* describe two different screening concepts in the
DLPNO implementation:

* **Strong and weak pairs** refer to individual occupied-orbital pairs
  :math:`ij`, classified from their estimated pair-correlation energies. This
  distinction controls whether a pair is treated at the coupled-cluster level,
  retained only through a lower-level correction, or eliminated. As noted above,
  the separate weak-pair interval is collapsed for DLPNO-CCSDT and higher methods,
  so every pair retained by those methods is a strong pair.

* **Strong and weak triplets or quadruplets** refer instead to energy-ranked
  occupied-orbital tuples :math:`ijk` or :math:`ijkl`. After the semicanonical
  :math:`(T0)` or :math:`(Q0)` contributions are evaluated, tuples are sorted by
  the magnitudes of those contributions. Triplets are classified by the leading
  contribution to their signed total, while quadruplets are classified by the
  leading 90% of the summed absolute contribution so cancellation cannot make
  the running criterion non-monotonic. The remainder are labeled weak. These
  labels can select independently configurable TNO or QNO truncation thresholds;
  they do **not** imply that a tuple contains any weak occupied-orbital pairs.

For an iterative DLPNO-CCSD(T) calculation that is the final requested method,
|dlpno__t_cut_tno_strong| (default ``1.0e-8``) and
|dlpno__t_cut_tno_weak| (default ``1.0e-7``) are the absolute TNO occupation
cutoffs for the two energy-ranked triplet classes. If DLPNO-CCSDT or a still
higher method is requested, both values are superseded by
|dlpno__t_cut_tno_full| (default ``1.0e-7``) while constructing the iterative
:math:`(T)` amplitudes that precede the full-:math:`T_3` calculation.
For backward compatibility, ``T_CUT_TNO_STRONG_SCALE`` and
``T_CUT_TNO_WEAK_SCALE`` remain available: unless the corresponding absolute
cutoff is explicitly set, the effective cutoff is that scale multiplied by
|dlpno__t_cut_tno|. An explicitly set absolute cutoff takes precedence.

Analogously, when iterative DLPNO-CCSDT(Q) is the final requested method,
|dlpno__t_cut_qno_strong| (default ``3.33e-6``) and
|dlpno__t_cut_qno_weak| (default ``3.33e-6``) are the absolute QNO occupation
cutoffs for the energy-ranked quadruplet classes. The defaults are equal because
numerical tests found no meaningful accuracy benefit from tightening the strong
class alone; the separate options remain available for explicit tuning. If
DLPNO-CCSDTQ is requested, both are superseded by |dlpno__t_cut_qno_full|
(default ``3.33e-6``) while constructing the iterative :math:`(Q)` amplitudes
that precede the full-:math:`T_4` calculation. The semicanonical
``dlpno-ccsd(t0)`` and
``dlpno-ccsdt(q0)`` variants do not perform the corresponding iterative step.
These semicanonical-only controls apply only when the associated perturbative
method is the final requested method: DLPNO-CCSDT and higher always form the
iterative :math:`(T)` amplitudes, and DLPNO-CCSDTQ always forms the iterative
:math:`(Q)` amplitudes, before proceeding to full amplitudes of the next rank.

.. warning::

   User-specified ``T_CUT_TNO_STRONG`` or ``T_CUT_TNO_WEAK`` values are ignored
   when DLPNO-CCSDT or higher is requested; ``T_CUT_TNO_FULL`` controls both
   triplet classes instead. Likewise, user-specified ``T_CUT_QNO_STRONG`` or
   ``T_CUT_QNO_WEAK`` values are ignored for DLPNO-CCSDTQ, where
   ``T_CUT_QNO_FULL`` controls both quadruplet classes. |PSIfour| prints a runtime
   warning when explicitly set values are superseded in this way.

DLPNO-CCSDT
~~~~~~~~~~~

Singles and doubles amplitudes (:math:`t_i^a` and :math:`t_{ij}^{ab}`) in pair natural
orbital (PNO) spaces, together with triples amplitudes (:math:`t_{ijk}^{abc}`) in triples
natural orbital (TNO) spaces from a preceding DLPNO-CCSD(T) computation, provide the
initial guess for solving the CCSDT residual equations [Jiang:2025:2386]_. For tractability,
the triples amplitudes are commonly recomputed with the looser TNO occupation cutoff
|dlpno__t_cut_tno_full| (default ``1.0e-7``) than is used for the perturbative triples
correction. DLPNO-CCSDT is asymptotically linear-scaling, and near-linear behavior can
emerge at relatively small system sizes for sparse systems such as water clusters.

An example input file for a DLPNO-CCSDT computation is::

   memory 2 GB

   molecule h2o {
   0 1
   O
   H 1 0.96
   H 1 0.96 2 104.5
   symmetry c1
   }

   set basis cc-pvdz
   set scf_type df
   set freeze_core true
   set pno_convergence very_tight
   set t_cut_pairs 1.0e-8
   set t_cut_tno_full 1.0e-7

   energy('dlpno-ccsdt')

DLPNO-CCSDT(Q)
~~~~~~~~~~~~~~

The pair-natural-orbital framework can also be applied to interacting occupied-orbital
quadruplets (:math:`ijkl`). Quadruples natural orbitals (QNOs) are obtained by
diagonalizing an averaged pair-density matrix and are retained according to their
occupation numbers using |dlpno__t_cut_qno| (default ``3.33e-7``).

.. math::
   :label: Quadruples Density Matrix

   D_{ijkl} &= \frac{1}{6} [D_{ij} + D_{jk} + D_{ik} + D_{il} + D_{jl} + D_{kl}].

More information on the pair-density matrices can be found in the
:ref:`DLPNO-MP2 <sec:dlpnomp2>` documentation. The perturbative quadruples amplitudes
[Bomble:2005:054101]_ are evaluated in QNO spaces using singles, doubles, and triples
amplitudes from a preceding DLPNO-CCSDT computation [Jiang:2025:144102]_. The DLPNO
formulation reduces the asymptotic system-size scaling of the canonical
:math:`\mathcal{O}(N^9)` CCSDT(Q) correction to linear scaling.

An example input file for an iterative DLPNO-CCSDT(Q) computation is::

   memory 2 GB

   molecule h2o {
   0 1
   O
   H 1 0.96
   H 1 0.96 2 104.5
   symmetry c1
   }

   set basis cc-pvdz
   set scf_type df
   set freeze_core true
   set pno_convergence very_tight
   set t_cut_pairs 1.0e-8
   set t_cut_tno_full 1.0e-7
   set t_cut_qno 3.33e-7

   energy('dlpno-ccsdt(q)')  # Use dlpno-ccsdt(q0) for the semicanonical (Q0) correction

Representative Semibullvalene Calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Cope rearrangement of semibullvalene provides a representative small-molecule
application of post-CCSD(T) local correlation. Fishman *et al.* used this reaction to
benchmark canonical and LNO-based higher-order coupled-cluster methods
[Fishman:2026:3233]_. A corresponding frozen-core Psi4 calculation used the same
equilibrium and transition-state geometries and modified ``cc-pVDZ(p,p)`` primary basis,
``VERY_TIGHT`` PNO settings, 16 OpenMP threads, and an AMD EPYC 9534 processor. The
input, output, and timer files are available with the
`PR 3506 benchmark archive <https://github.com/psi4/psi4/pull/3506#issuecomment-5415257893>`__.

The table below compares post-CCSD(T) contributions to the barrier height. All values are
in :math:`\mathrm{kcal\,mol^{-1}}`; the error column is the absolute error in the combined
:math:`\mathrm{CCSDT(Q)}-\mathrm{CCSD(T)}` contribution relative to the canonical value.
Comparing increments, rather than total energies from different local-correlation models,
isolates the higher-order triples and quadruples effects of interest.

.. list-table:: Semibullvalene post-CCSD(T) barrier-height contributions.
   :header-rows: 1
   :widths: 27 18 18 22 15

   * - Method
     - :math:`\mathrm{CCSDT}-\mathrm{CCSD(T)}`
     - :math:`\mathrm{CCSDT(Q)}-\mathrm{CCSDT}`
     - :math:`\mathrm{CCSDT(Q)}-\mathrm{CCSD(T)}`
     - Absolute error
   * - Canonical [Fishman:2026:3233]_
     - 0.139
     - -0.268
     - -0.129
     - --
   * - LNO Loose [Fishman:2026:3233]_
     - 0.104
     - -0.220
     - -0.116
     - 0.013
   * - LNO Normal [Fishman:2026:3233]_
     - 0.152
     - -0.255
     - -0.103
     - 0.026
   * - Psi4 DLPNO ``VERY_TIGHT``
     - 0.110
     - -0.264
     - -0.154
     - 0.025

For this example, the DLPNO quadruples contribution differs from the canonical value by
only :math:`0.004\ \mathrm{kcal\,mol^{-1}}`. Its combined post-CCSD(T) error is
:math:`0.025\ \mathrm{kcal\,mol^{-1}}`, comparable to the
:math:`0.026\ \mathrm{kcal\,mol^{-1}}` LNO Normal error and somewhat larger than the
:math:`0.013\ \mathrm{kcal\,mol^{-1}}` LNO Loose error. The resulting Psi4
DLPNO-CCSDT(Q) barrier is :math:`6.545\ \mathrm{kcal\,mol^{-1}}`.

The calculation contains 16 atoms and 56 electrons. The 112-function primary basis gives
28 doubly occupied and 84 canonical virtual orbitals; freezing eight carbon core orbitals
leaves 20 active occupied orbitals. The correlation and SCF density-fitting bases contain
560 and 744 auxiliary functions, respectively. Local natural-orbital compression and the
reported resource estimates are summarized below. Parentheses give minimum--maximum
numbers of retained natural orbitals.

.. list-table:: Local-space sizes and resources for the semibullvalene calculations.
   :header-rows: 1
   :widths: 44 28 28

   * - Quantity
     - Equilibrium structure
     - Transition state
   * - Active occupied / canonical virtual orbitals
     - 20 / 84
     - 20 / 84
   * - Retained triplets
     - 1520 of 1520
     - 1520 of 1520
   * - Strong triplets
     - 619 (40.7%)
     - 672 (44.2%)
   * - Full-:math:`T_3` TNOs per triplet, average (range)
     - 32 (18--44)
     - 35 (18--53)
   * - Retained quadruplets after energy screening
     - 6670 of 8455 (78.9%)
     - 7048 of 8455 (83.4%)
   * - Strong quadruplets among retained quadruplets
     - 2055 (30.8%)
     - 2603 (36.9%)
   * - :math:`(Q0)` QNOs per quadruplet, average (range)
     - 32 (21--41)
     - 35 (21--48)
   * - Iterative-:math:`(Q)` QNOs per quadruplet, average (range)
     - 17 (12--24)
     - 19 (12--28)
   * - Estimated resident memory
     - 14.6 GB
     - 20.8 GB
   * - Estimated peak memory
     - 16.9 GB
     - 24.6 GB
   * - Wall time, 16 threads
     - 56 min 44 s
     - 80 min 48 s

The configured ``memory 200 GB`` value was an allocation ceiling, not the amount consumed.
The values above are Psi4's internal estimates; the archived timer files do not contain an
operating-system maximum-resident-set measurement. On both structures, nearly all wall time
is concentrated in three stages. Construction of the quadruples source and semicanonical
:math:`T_4` amplitudes (the ``gamma ijkl`` timer) required 1573 s (46%) and 2337 s (48%)
for the equilibrium and transition structures, respectively. The iterative CCSDT stage
required 1260 s (37%) and 1498 s (31%), with construction of the :math:`R_3` residual
accounting for most of that time. Finally, the iterative :math:`(Q)` stage required 481 s
(14%) and 918 s (19%); its repeated :math:`(Q)` energy contractions were the dominant
substep. SCF, PNO construction, lower-rank CCSD iterations, and TNO/QNO transformations
were comparatively inexpensive for this case.

.. note::

   This archived benchmark predates the replacement of QNO scale factors by separate
   absolute strong- and weak-quadruplet cutoffs. It effectively used
   ``T_CUT_QNO_STRONG = T_CUT_QNO_WEAK = 3.33e-6``. The current defaults now
   match those archived settings, so no explicit QNO-class overrides are needed
   to reproduce the tabulated energies and local-space sizes.

DLPNO-CCSDTQ
~~~~~~~~~~~~

DLPNO-CCSDTQ iteratively solves the full quadruples-amplitude equations in local QNO
spaces [Jiang:2026:2825]_. The full-quadruples QNO space is controlled by
|dlpno__t_cut_qno_full| (default ``3.33e-6``). Extended pair natural orbitals (XPNOs)
are additionally constructed to project selected quadruples amplitudes into extended
pair domains for the computationally dominant contractions. Their occupation cutoff is
controlled by |dlpno__t_cut_xpno|, which defaults to three times
|dlpno__t_cut_qno_full| (approximately ``1.0e-5``).

An example input file for a DLPNO-CCSDTQ computation is::

   memory 2 GB

   molecule h2o {
   0 1
   O
   H 1 0.96
   H 1 0.96 2 104.5
   symmetry c1
   }

   set basis cc-pvdz
   set scf_type df
   set freeze_core true
   set pno_convergence very_tight
   set t_cut_pairs 1.0e-8
   set t_cut_tno_full 1.0e-7
   set t_cut_qno 3.33e-7
   set t_cut_qno_full 3.33e-6
   set t_cut_xpno 1.0e-5

   energy('dlpno-ccsdtq')

Representative Benzene-Dimer DLPNO-CCSDTQ Calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A cofacial benzene sandwich dimer with a 3.9 Angstrom interplanar separation
provides a representative noncovalent test of the complete local hierarchy through
DLPNO-CCSDTQ. The archived calculations used the frozen-core ``cc-pVDZ`` basis,
``VERY_TIGHT`` PNO settings, Pipek--Mezey localized orbitals, direct SCF, energy
and residual convergence thresholds of ``1.0e-7``, and 48 OpenMP threads on an
AMD Genoa node. The dimer and counterpoise (CP) monomer calculations each used
228 primary basis functions; 114 of the CP-monomer functions were ghosts on the
partner fragment. The non-counterpoise (noCP) monomer used 114 primary basis
functions.

For the symmetric dimer, the interaction energies were evaluated as

.. math::

   \Delta E_\mathrm{int}^\mathrm{CP}
      = E_{AB}^{AB} - 2 E_A^{AB},
   \qquad
   \Delta E_\mathrm{int}^\mathrm{noCP}
      = E_{AB}^{AB} - 2 E_A^A.

Negative values are attractive. The absolute CP and noCP interaction energies
differ because the latter retain basis-set superposition error, but the focal-point
increments from the higher-order correlation models are much less sensitive to the
counterpoise convention.

.. list-table:: Benzene-dimer interaction energies and post-CCSD(T) focal-point increments.
   :header-rows: 1
   :widths: 46 18 18 18

   * - Quantity
     - CP
     - noCP
     - Absolute difference
   * - :math:`\Delta E_\mathrm{int}[\mathrm{CCSD(T)}]`
     - 0.207
     - -1.081
     - 1.288
   * - :math:`\Delta E_\mathrm{int}[\mathrm{CCSDT}]`
     - 0.256
     - -1.009
     - 1.265
   * - :math:`\Delta E_\mathrm{int}[\mathrm{CCSDT(Q)}]`
     - 0.217
     - -1.055
     - 1.271
   * - :math:`\Delta E_\mathrm{int}[\mathrm{CCSDTQ}]`
     - -0.183
     - -1.457
     - 1.273
   * - :math:`\mathrm{CCSDT}-\mathrm{CCSD(T)}`
     - 0.049
     - 0.071
     - 0.022
   * - :math:`\mathrm{CCSDT(Q)}-\mathrm{CCSDT}`
     - -0.039
     - -0.045
     - 0.006
   * - :math:`\mathrm{CCSDTQ}-\mathrm{CCSDT(Q)}`
     - -0.400
     - -0.402
     - 0.002
   * - :math:`\mathrm{CCSDTQ}-\mathrm{CCSD(T)}`
     - -0.391
     - -0.376
     - 0.015

All values in the preceding table are in :math:`\mathrm{kcal\,mol^{-1}}`.
In particular, the full-quadruples increment relative to perturbative quadruples
changes by only :math:`0.002\ \mathrm{kcal\,mol^{-1}}` between CP and noCP,
while the complete post-CCSD(T) contribution changes by
:math:`0.015\ \mathrm{kcal\,mol^{-1}}`. Slight differences from previously
circulated benchmark values are expected because those calculations used a
different cutoff set. The table above records the internally consistent CP/noCP
results from the supplied calculation archive.

The resource profile is summarized below. The memory values are Psi4's internal
peak estimates for each stage, rather than operating-system maximum-resident-set
measurements. The configured ``memory 2880 GB`` value was only an allocation
ceiling. All three calculations wrote the largest TNO- and QNO-basis density-fitting
integrals to disk.

.. list-table:: Local-space sizes and resources for the benzene-dimer calculations.
   :header-rows: 1
   :widths: 43 19 19 19

   * - Quantity
     - Dimer
     - CP monomer
     - noCP monomer
   * - Primary basis functions
     - 228
     - 228 (114 ghost)
     - 114
   * - Active occupied / canonical virtual orbitals
     - 30 / 186
     - 15 / 207
     - 15 / 93
   * - Retained triplets
     - 4132 of 4930
     - 665 of 665
     - 665 of 665
   * - Strong triplets
     - 568 (13.7%)
     - 270 (40.6%)
     - 270 (40.6%)
   * - Retained quadruplets
     - 8652 of 40020
     - 2552 of 2835
     - 2550 of 2835
   * - Strong quadruplets
     - 1625 (18.8%)
     - 752 (29.5%)
     - 752 (29.5%)
   * - Estimated DLPNO-CCSDT peak memory
     - 89.1 GB
     - 48.0 GB
     - 42.4 GB
   * - Estimated DLPNO-CCSDT(Q) peak memory
     - 113.5 GB
     - 41.4 GB
     - 40.7 GB
   * - Estimated DLPNO-CCSDTQ peak memory
     - 446.7 GB
     - 324.5 GB
     - 286.7 GB
   * - DLPNO-CCSDT wall time
     - 51 min 54 s
     - 6 min 39 s
     - 6 min 18 s
   * - ``gamma ijkl`` wall time
     - 38 min 14 s
     - 10 min 27 s
     - 10 min 20 s
   * - Iterative-(Q) wall time
     - 16 min 21 s
     - 5 min 25 s
     - 5 min 27 s
   * - DLPNO-CCSDTQ wall time
     - 4 d 2 h 25 min
     - 10 h 19 min
     - 9 h 19 min
   * - Total wall time
     - 4 d 4 h 13 min
     - 10 h 43 min
     - 9 h 42 min
   * - Full-quadruples iterations
     - 69
     - 24
     - 24

The dimer timing is conservative for this archived input: it explicitly set
|dlpno__extrapolate_t4| to ``FALSE``, excluding the :math:`T_4` amplitudes
from DIIS, whereas both monomer calculations used the default ``TRUE`` setting.
Consequently, the dimer required 69 full-quadruples iterations and spent about
98 hours in the DLPNO-CCSDTQ stage alone. Including :math:`T_4` in DIIS
normally improves convergence at the cost of larger DIIS vectors. The monomer
comparison also shows the cost of the ghost basis: the CP monomer required about
10% more wall time and a 13% larger estimated DLPNO-CCSDTQ peak than the noCP
monomer.

.. note::

   These calculations used ``T_CUT_TNO_STRONG = T_CUT_TNO_WEAK = 1.0e-7``,
   ``T_CUT_QNO_STRONG = T_CUT_QNO_WEAK = T_CUT_QNO_FULL = 3.33e-6``, and
   ``T_CUT_XPNO = 9.99e-6``. Because natural-orbital cutoffs control both the
   retained local spaces and the computational cost, timings and memory
   requirements should be interpreted as order-of-magnitude planning data for
   these settings and this hardware, not as universal performance guarantees.

PNO Convergence Settings
------------------------

Here we present a table of the PNO convergence settings, parameters, and recommended use cases. Most of these parameters and
settings are similar to what is found in ORCA, with two added parameters (|dlpno__t_cut_trace| and |dlpno__t_cut_energy|) to
increase the robustness of the PNO space. These added parameters truncate by percent recovery of the total occupation number,
as well as the percentage energy recovery of the PNOs compared to the non-truncated basis.

.. _`table:pno_convergence`:

.. table:: PNO convergence levels given in |Psifour|

   +--------------------------+------------+-------------+--------------+----------+-----------+-------------+---------------------------+
   | |dlpno__pno_convergence| | T_CUT_PNO  | T_CUT_TRACE | T_CUT_ENERGY | T_CUT_DO | T_CUT_MKN | T_CUT_PAIRS | Recommended Applications  |
   +==========================+============+=============+==============+==========+===========+=============+===========================+
   | Loose                    | 1.0e-6     | 0.9         | 0.9          | 2e-2     | 1e-3      | 1e-3        | High-throughput screening |
   +--------------------------+------------+-------------+--------------+----------+-----------+-------------+---------------------------+
   | Normal                   | 3.33e-7    | 0.99        | 0.99         | 1e-2     | 1e-3      | 1e-4        | Thermochemistry           |
   +--------------------------+------------+-------------+--------------+----------+-----------+-------------+---------------------------+
   | Tight                    | 1.0e-7     | 0.999       | 0.997        | 5e-3     | 1e-3      | 1e-5        | Non-covalent Interactions |
   +--------------------------+------------+-------------+--------------+----------+-----------+-------------+---------------------------+
   | Very_Tight               | 1.0e-8     | 0.999       | 0.997        | 5e-3     | 1e-4      | 1e-6        | Benchmarking, Focal Point |
   +--------------------------+------------+-------------+--------------+----------+-----------+-------------+---------------------------+

Practical Advice
----------------

* DLPNO-CCSD/(T) is almost always faster than the corresponding canonical CCSD/(T) computation,
  and computations involving DLPNO-CCSD/(T) are encouraged to be performed on large molecules 
  as a more accurate alternative to DFT.

* For most computations, |dlpno__pno_convergence| ``TIGHT`` is recommended, especially 
  those involving non-covalent interactions. For larger systems where ``TIGHT`` is too expensive, 
  ``NORMAL`` for |dlpno__pno_convergence| while setting |dlpno__t_cut_pairs| to ``1.0e-5`` is recommended. 
  This has been shown to yield errors on the order of kJ/mol for non-covalent interactions [Jiang:2024:082502]_.

* Based on user allocated memory, disk/core storage for various integrals and tensors
  for DLPNO-CCSD/(T) are automatically determined. There is no need to toggle with the disk/core
  options for the average user.

* DLPNO-CCSDT, DLPNO-CCSDT(Q), and DLPNO-CCSDTQ substantially extend the molecular sizes
  accessible to post-CCSD(T) benchmarks relative to their canonical counterparts. They remain
  considerably more demanding than DLPNO-CCSD(T), and the linear-scaling regime may not be reached
  at system sizes practical on a typical laboratory workstation.

* Memory capacity and available disk space often determine the largest feasible higher-order
  calculation. Resource requirements depend strongly on molecular sparsity, basis set, natural-orbital
  thresholds, and whether expensive intermediates are retained in memory or written to disk. For an
  iterative DLPNO-CCSDT(Q) calculation, |PSIfour| automatically moves the :math:`\Gamma_{ijkl}` source,
  :math:`T_{ijkl}` amplitudes, and reusable energy intermediates to disk when its memory estimate exceeds
  the available allocation. Expert users can force this behavior with
  |dlpno__write_quadruples_intermediates|.

* DLPNO-CCSDT, DLPNO-CCSDT(Q), and DLPNO-CCSDTQ require a |PSIfour| build with the optional
  Einsums tensor library enabled (``ENABLE_Einsums=ON``). DLPNO-MP2, DLPNO-CCSD,
  and DLPNO-CCSD(T) remain available when Einsums support is disabled.

* For higher-order calculations, the recommended starting settings are
  |dlpno__pno_convergence| ``VERY_TIGHT`` and |dlpno__t_cut_pairs| = ``1.0e-8``.
  These values are applied automatically by the higher-order drivers. Custom thresholds should be
  validated carefully for the intended application.

* DLPNO-CCSDT(Q) is particularly useful in composite or focal-point schemes, where the
  higher-order correction is formed as the difference between DLPNO-CCSDT(Q) and a canonical
  DF-CCSD(T) or rank-corrected DLPNO-CCSD(T) energy. Such combinations can provide favorable
  cancellation of local-correlation truncation errors [Jiang:2025:144102]_.

* In DLPNO methods, it is recommended to freeze core orbitals (by setting |globals__freeze_core|
  to ``True``), since core excitations are known to be more sensitive to PNO truncations than
  valence truncations. If a non-frozen core computation is requested, all PNOs corresponding to core-core
  or core-virtual pairs have cutoffs scaled by ``T_CUT_PNO_CORE_SCALE`` (default ``1.0e-2``).

* DLPNO does not yet use molecular point-group symmetry. The SCF reference may be computed in its detected
  point group, after which the converged wavefunction is transformed to C1 for the DLPNO calculation. A
  user-supplied higher-symmetry SCF reference wavefunction is transformed in the same way.

* At this time, all DLPNO coupled-cluster methods are available only for closed-shell RHF energy computations.

Computation Size Limits
-----------------------

* Since DLPNO-CCSD(T) is linear-scaling, with access to sufficient computing resources, a DLPNO-CCSD(T) computation 
  can be possible with any system. In fact, for larger systems, Hartree-Fock becomes the bottleneck (not the coupled-cluster)!
  Below, we tabulate the roughly projected limits (number of atoms) of DLPNO-CCSD(T) across different access to hardware. All these limits 
  are with a standard polarized double zeta basis set. For larger basis sets, divide by 3 for an increase in cardinality, 
  2 for full set of diffuse functions, and 1.5 for partial diffuse functions (e.g. if the size limit is 100
  for cc-pVDZ, expect 30 for cc-pVTZ, 50 for aug-cc-pVDZ, and 70 for jun-cc-pVDZ). Estimates for larger amounts of RAM should be
  taken with a grain of salt, as the computation may be theoretically possible with DLPNO-CCSD, but may be hindered by the cost of
  the preceeding Hartree-Fock computation (in time or memory). New approaches to make Hartree-Fock more efficient in |PSIfour| are
  currently under investigation.

.. _`table:size_limits`:

.. table:: Expected size limits (number of atoms) of DLPNO-CCSD(T) computation

   +------------------+------------+-------------+-------------+
   | Hardware         | RAM        | Normal      | Tight       |
   +==================+============+=============+=============+
   | Home Desktop     | 32 GB      | 90-100      | 80-90       |
   +------------------+------------+-------------+-------------+
   | Lab Workstation  | 64 GB      | 120-150     | 100-120     |
   +------------------+------------+-------------+-------------+
   | Lab Workstation  | 192 GB     | 300-400     | 200-300     |
   +------------------+------------+-------------+-------------+
   | Lab Cluster      | 512 GB     | 700+        | 400+        |
   +------------------+------------+-------------+-------------+
   | Lab/HPC Cluster  | 1 TB       | 1000+       | 700+        |
   +------------------+------------+-------------+-------------+
   | HPC Cluster      | 3 TB       | 1500+       | 1000+       |
   +------------------+------------+-------------+-------------+

Key Differences with DLPNO-CCSD(T) in ORCA
------------------------------------------

While the DLPNO-CCSD(T) formulation in |PSIfour| is heavily inspired by the original method
proposed by Neese and coworkers in ORCA [Riplinger:2013:034106]_ [Riplinger:2013:134101]_ [Riplinger:2016:024109]_ 
[Guo:2018:011101]_, |PSIfour| employs different algorithms for certain parts of the procedure. 
Both represent linear-scaling CCSD(T) algorithms solved in the the local pair natural orbital basis, with convergence
to canonical CCSD(T) results as the local tolerances are tightened. However, the manner in which the PNO spaces are truncated as
well as how the CCSD equations are solved are different. Notable differences in implementation between the two algorithms are
highlighted below:

* The most notable difference is that the DLPNO-CCSD equations in |PSIfour| utilize T1-dressed Hamiltonian and Fock matrix elements,
  which significantly simplifies the number of working equations as well as the number of mathematical operations involved in solving
  the CCSD residual equations. Because of the reduced number of intermediates that are required, we find that the DLPNO-CCSD code in 
  |PSIfour| potentially uses less RAM than the ORCA formulation for a given |dlpno__pno_convergence|.
  The runtimes are expected to be similar.

* |PSIfour| often recovers slightly more CCSD correlation energy than ORCA due to additional PNO cutoffs |dlpno__t_cut_trace| 
  and |dlpno__t_cut_energy| used in addition to the normal |dlpno__t_cut_pno|. These result in larger and more robust PNO spaces at a given 
  |dlpno__pno_convergence| in |PSIfour|. The difference is typically small in terms of absolute energies (on the order of ``1.0e-3`` to ``1.0e-4`` Hartrees),
  with agreement on the order of 99.95% or better. Due to the different T1 formulation, it is not possible to exactly match ORCA and |PSIfour| 
  DLPNO-CCSD energies by adjusting keywords.

* |PSIfour| also couples MP2 weak pair amplitudes to CCSD strong pair amplitudes in solving the CCSD residual equations. In our
  formulation, this does not add significantly more time and memory since the most expensive algorithmic steps result from self-coupling 
  terms, such as the particle-particle ladder (:math:`R_{ij}^{ab} \mathrel{+}= B^{Q}_{ac} t_{ij}^{cd} B^{Q}_{bd}`). This is likely an
  additional source of the increased recovery of correlation energy in |PSIfour| compared to ORCA.

* The more robust PNO space combined with strong pair/weak pair couplings often give |PSIfour| a slight edge in terms of relative energies
  compared to canonical CCSD(T) at a given |dlpno__pno_convergence|, with comparable runtimes. This is showcased by the comparing conformation
  energies on large water clusters (16-17 waters), as shown in [Jiang:2024:082502]_.

* Another difference with ORCA is the use of Full LMP2 prescreening across all PNO convergence levels, as opposed to only TightPNO in ORCA.
  For a more direct comparison, one can set ``UseFullLMP2Guess true`` in the corresponding DLPNO-CCSD(T) input file in ORCA. For the triples 
  computation, ORCA defaults to the semicanonical (T0) approach given in [Riplinger:2013:134101]_ when specifying ``DLPNO-CCSD(T)`` in the
  input file, while our code defaults to the iterative (T) computation given in [Guo:2018:011101]_ with an ``energy('dlpno-ccsd(t)')`` call.
  To perform the iterative (T) computation in ORCA, one needs to specify ``DLPNO-CCSD(T1)``. In our code, the semicanonical (T0) 
  implementation of triples can be requested through ``energy('dlpno-ccsd(t0)')``.

* In benchmarking studies, users are encouraged to use both ORCA and |PSIfour|'s implementation. The main advantage of the code in ORCA
  is the capability to run across different nodes through MPI, while our code is designed for an optimal performance on a single node
  through OpenMP. Both codes should converge to similar answers for relative energies with larger basis sets and tighter cutoffs, especially 
  when extrapolated to the complete basis set (CBS) limit.
