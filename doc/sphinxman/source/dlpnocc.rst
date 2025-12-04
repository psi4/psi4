.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2025 The Psi4 Developers.
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
   
   .. _`sec:dlpnocc`:
   
=================================================
DLPNO-CC: Domain-Based Local Pair Natural Orbital
=================================================

.. index::
   single: DLPNO-CCSD(T)

.. _`sec:dlpnoccsd_t`:

DLPNO-CCSD(T): Domain-Based Local Pair Natural Orbital CCSD(T)
==============================================================

.. codeauthor:: Andy Jiang 
.. sectionauthor:: Andy Jiang

*Module:* :ref:`Keywords <apdx:dlpno>`, :ref:`PSI Variables <apdx:dlpno_psivar>`, :source:`DLPNOCCSD(T) <psi4/src/psi4/dlpno>`

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
to canonical CCSD(T). These errors are controllable through user set parameters. Computations on large molecules, 
such as crambin (636 atoms) [Riplinger:2013:134101]_ and insulin [Jiang:2024:082502]_ have been performed
with this algorithm. This algorithm was originally implemented in the ORCA software package, but it
is now available in |PSIfour|!

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

PNO Convergence Settings
------------------------

Here we present a table of the PNO convergence settings, paramaters, and recommended use cases. Most of these parameters and
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

* In DLPNO methods, it is recommended to freeze core orbitals (by setting |globals__freeze_core|
  to ``True``), since core excitations are known to be more sensitive to PNO truncations than
  valence truncations. If a non-frozen core computation is requested, all PNOs corresponding to core-core
  or core-virtual pairs have cutoffs scaled by ``T_CUT_PNO_CORE_SCALE`` (default ``1.0e-2``).

* Note that DLPNO does not yet have molecular point group symmetry implemented and will run in C1 symmetry.

* At this time, DLPNO-CCSD/(T) is only available for closed-shell RHF computations.

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