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
   single: DLPNO-CCSDT(Q)

.. _`sec:dlpnoccsdt_q`:

DLPNO-CCSDT(Q): Domain-Based Local Pair Natural Orbital CCSDT(Q)
================================================================

.. codeauthor:: Andy Jiang 
.. sectionauthor:: Andy Jiang

*Module:* :ref:`Keywords <apdx:dlpno>`, :ref:`PSI Variables <apdx:dlpno_psivar>`, :source:`DLPNOCCSDT(Q) <psi4/src/psi4/dlpno>`

Introduction
------------

Recent work [Jiang:2025:2386]_ [Jiang:2025:144102]_ has extended the domain-based local
pair natural orbital (DLPNO) framework to coupled-cluster methods beyond CCSD(T). The
same locality principles previously applied to DLPNO-CCSD(T)
[Riplinger:2013:034106]_ [Riplinger:2013:134101]_ [Riplinger:2016:024109]_
[Guo:2018:011101]_ [Jiang:2024:082502]_ are used for the higher-order methods described here.

Singles and doubles amplitudes (:math:`t_{i}^{a}` and :math:`t_{ij}^{ab}`) in the pair natural orbital (PNO) basis,
together with triples amplitudes (:math:`t_{ijk}^{abc}`) in the triples natural orbital (TNO) basis from a preceding
:ref:`DLPNO-CCSD(T) <sec:dlpnocc>` computation, provide the initial guess for solving the CCSDT residual equations
[Jiang:2025:2386]_. For tractability, the triples amplitudes are commonly recomputed with the looser TNO occupation
cutoff |dlpno__t_cut_tno_full| (default ``1.0e-7``) than is used for the perturbative triples correction.
DLPNO-CCSDT is asymptotically linear-scaling, and near-linear scaling behavior can emerge at relatively
small system sizes for sparse systems such as water clusters.

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
   set scf_type direct
   set freeze_core True
   set pno_convergence very_tight
   set t_cut_pairs 1.0e-8
   set t_cut_tno_full 1.0e-7
   
   energy('dlpno-ccsdt')

The pair natural orbital framework can also be applied to interacting occupied-orbital quadruplets (``ijkl``).
Quadruples natural orbitals (QNOs) are obtained by diagonalizing the averaged pair-density matrix, and QNOs are
retained according to their occupation numbers using |dlpno__t_cut_qno| (default ``3.33e-7``).

.. math::
   :label: Quadruples Density Matrix

   D_{ijkl} &= \frac{1}{6} [D_{ij} + D_{jk} + D_{ik} + D_{il} + D_{jl} + D_{kl}].

More information on the pair density matrices can be found in the :ref:`DLPNO-MP2 <sec:dlpnomp2>` documentation.
The (Q) quadruples amplitudes [Bomble:2005:054101]_ are evaluated in the QNO basis using singles, doubles, and
triples amplitudes from a preceding DLPNO-CCSDT computation [Jiang:2025:144102]_. The DLPNO formulation reduces
the asymptotic system-size scaling of the canonical :math:`\mathcal{O}(N^{9})` CCSDT(Q) correction to linear scaling.

An example input file for a DLPNO-CCSDT(Q) computation is::
   
   memory 2 GB

   molecule h2o {
   0 1
   O
   H 1 0.96
   H 1 0.96 2 104.5
   symmetry c1
   }
   
   set basis cc-pvdz
   set scf_type direct
   set freeze_core True
   set pno_convergence very_tight
   set t_cut_pairs 1.0e-8
   set t_cut_tno_full 1.0e-7
   set t_cut_qno 3.33e-7
   
   energy('dlpno-ccsdt(q)') # dlpno-ccsdt(q0) for the semicanonical (Q0) computation

Practical Advice
----------------

* DLPNO-CCSDT and DLPNO-CCSDT(Q) substantially extend the system sizes accessible to higher-order
  coupled-cluster benchmarks relative to their canonical counterparts.

* For most systems other than highly sparse clusters, the asymptotic linear-scaling regime may not be
  reached at system sizes practical on a typical laboratory workstation.

* These methods remain computationally demanding because they target high-accuracy approximations to
  canonical CCSDT and CCSDT(Q). Calculations on very large systems, such as proteins, are generally not
  practical with current computational resources.

* Memory capacity and available disk space often determine the largest feasible calculation. As a rough
  guide, a typical laboratory workstation can treat a benzene dimer with a polarized double-zeta basis,
  although resource requirements depend strongly on geometry, thresholds, and hardware.

* The recommended settings are ``VERY_TIGHT`` |dlpno__pno_convergence| and
  |dlpno__t_cut_pairs| = ``1.0e-8``. These values are applied automatically for the higher-order methods;
  custom thresholds should be validated carefully for the intended application.

* Users are encouraged to use DLPNO-CCSDT(Q) as part of composite or focal point schemes, with 
  the T(Q) correction being computed through the difference of the DLPNO-CCSDT(Q) energy and a
  canonical DF-CCSD(T) or rank-corrected DLPNO-CCSD(T) energy (either canonical MP2,
  tighter |dlpno__t_cut_tno| cutoff, or a CPS extrapolation). This approach has been shown to yield
  favorable error cancellation.

* In DLPNO methods, it is recommended to freeze core orbitals (by setting |globals__freeze_core|
  to ``True``), since core excitations are known to be more sensitive to PNO truncations than
  valence truncations.

* At this time, DLPNO-CCSDT/(Q) is only available for closed-shell RHF computations.