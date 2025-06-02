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

Practical Advice
----------------

* DLPNO-CCSD/(T) is almost always faster than the corresponding canonical CCSD/(T) computation,
  and computations involving DLPNO-CCSD/(T) are encouraged to be performed on large molecules 
  as a more accurate alternative to DFT.

* For most computations, |dlpno__pno_convergence| ``TIGHT`` is recommended, especially 
  those involving non-covalent interactions. For larger systems where ``TIGHT`` is too expensive, 
  ``NORMAL`` for |dlpno__pno_convergence| while setting |dlpno_t_cut_pairs| to ``1.0e-5`` is recommended. 
  This has been shown to yield errors on the order of kJ/mol for non-covalent interactions [Jiang:2024:082502]_.

* For the triples computation, the semicanonical (T0) approach given in [Riplinger:2013:134101]_ is NOT 
  recommended. As opposed to ORCA, the default computation for a ``dlpno-ccsd(t)`` call is |PSIfour| is 
  the iterative (T) computation given in [Guo:2018:011101]_ (i.e. DLPNO-CCSD(T1) in ORCA). The semicanonical (T0) 
  implementation of triples can be requested by the user through the ``energy('dlpno-ccsd(t0)')`` call.
  (NOTE: In ORCA, a DLPNO-CCSD(T) computation defaults to the semicanonical (T0) computation.)

* Due to differences in implementation and working equations, absolute DLPNO-CCSD(T) energies in |PSIfour| will not
  exactly match ORCA, as |PSIfours| |dlpno__pno_convergence| at a given level is designed to be more accurate than
  the same setting in ORCA (see results on water clusters in [Jiang:2024:082502]_).

* Based on user allocated memory, disk/core storage for various integrals and tensors
  for DLPNO-CCSD/(T) are automatically determined. There is no need to toggle with the disk/core
  options for the average user.

* In DLPNO methods, it is recommended to freeze core orbitals (by setting |globals__freeze_core|
  to ``True``), since core excitations are known to be more sensitive to PNO truncations than
  valence truncations.

* At this time, DLPNO-CCSD/(T) is only available for closed-shell RHF computations.