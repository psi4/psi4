.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2024 The Psi4 Developers.
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

In our recent works [Jiang:2025:2386]_ [Jiang:2025:144102]_, we have shown that the domain-based 
local pair natural orbital (DLPNO) approach successfully applied to gold-standard CCSD(T), 
[Riplinger:2013:034106]_ [Riplinger:2013:134101]_ [Riplinger:2016:024109]_ 
[Guo:2018:011101]_ [Jiang:2024:082502]_, can also be applied to higher levels of coupled-cluster theory.

Singles and doubles amplitudes (:math:`t_{i}^{a}` and :math:`t_{ij}^{ab}`) in the pair natural orbital (PNO) basis, 
and triples amplitudes (:math:`t_{ijk}^{abc}`) in the triples natural orbital (TNO) basis from a preceding 
:ref:`DLPNO-CCSD(T) <sec:dlpnoccsd_t>` computation are used as the starting guess as the solution of the
CCSDT residual equations [Jiang:2025:2386]_. The triples amplitudes are often recomputed at a looser TNO
tolerance than the one used for (T) |dlpno__t_cut_tno_full| (default ``1.0e-7``) for the sake of tractibility.
This DLPNO-CCSDT method is asymptotically linear-scaling, and near-linear scaling behavior can be observed
early for systems with significant sparsity, like water clusters. To run a DLPNO-CCSDT computation, a sample input
file is provided:

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

The pair natural orbital (PNO) framework can also be applied to interacting quadruplets (``ijkl``), computed
if the pairs ``ij``, ``jk``, ``ik``, ``il``, ``jl``, and ``kl`` are all interacting pairs. Quadruples natural
orbitals (QNOs), are formed as the eigenvalues of the quadruples density matrix, and then truncated based on
eigenvalues by the parameter |dlpno__t_cut_qno| (default ``3.33e-7``).

.. math::
   :label: Quadruples Density Matrix

   D_{ijkl} &= \frac{1}{6} [D_{ij} + D_{jk} + D_{ik} + D_{il} + D_{jl} + D_{kl}].

More information on the pair density matrices can be found in the :ref:`DLPNO-MP2 <sec:dlpnomp2>` documentation.
(Q) quadruples amplitudes [Bomble:2005:054101]_ are then solved in the QNO basis, using singles, doubles, and 
triples amplitudes from a preceding DLPNO-CCSDT computation, in our DLPNO-CCSDT(Q) algorithm [Jiang:2025:144102]_. 
Our algorithm brings the platinum-standard :math:`O(N^{9})` CCSDT(Q) algorithm, typically intractable for most molecules 
more than 10 atoms, down to asymptotic linear-scaling cost.

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

* DLPNO-CCSDT/(Q) makes the impossible possible! Platinum-standard theoretical benchmarks that would
  normally require millions or billions of years with canonical CCSDT(Q) are now possible with less
  than a couple weeks of computation time (or a day if you have access to lots of RAM/CPUs).

* For most systems (other than very sparse water clusters), asymptotic linear-scaling behavior 
  for DLPNO-CCSDT/(Q) will not be observed in the range of computations currently possible on the 
  typical lab workstation.

* DLPNO-CCSDT/(Q) computations can still take a while, after all, it is designed to accurately
  reproduce the result of one of the most expensive quantum chemistry algorithms in the world, to within
  a relative energy error of :math:`0.1 kcal mol^{-1}`. It probably is not (yet) possible to run a protein
  calculation with DLPNO-CCSDT/(Q).

* DLPNO-CCSDT/(Q) loves RAM and disk! It consumes computer resources like a growing adolescent. The computations
  you will be able to perform will be largely dictated by the computational resources you have access to!
  As a rule of thumb, a typical lab workstation can handle a benzene dimer computation with a polarized double
  zeta basis set!

* Always use ``VERY_TIGHT`` |dlpno__pno_convergence| for a DLPNO-CCSDT/(Q) computation, and |dlpno__t_cut_pairs|
  is set to ``1.0e-8`` to maintain the accuracy needed for typical applications. This is done by default, and users
  who attempt to use custom parameters for DLPNO-CCSDT/(Q) do so at their own peril!

* Users are encouraged to use DLPNO-CCSDT(Q) as part of composite or focal point schemes, with 
  the T(Q) correction being computed through the difference of the DLPNO-CCSDT(Q) energy and a
  canonical DF-CCSD(T) or rank-corrected DLPNO-CCSD(T) energy (either canonical MP2,
  tighter |dlpno__t_cut_tno| cutoff, or a CPS extrapolation). This approach has been shown to yield
  the best error cancellation!

* In DLPNO methods, it is recommended to freeze core orbitals (by setting |globals__freeze_core|
  to ``True``), since core excitations are known to be more sensitive to PNO truncations than
  valence truncations.

* At this time, DLPNO-CCSDT/(Q) is only available for closed-shell RHF computations.