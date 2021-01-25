.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2021 The Psi4 Developers.
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
.. # what does "single" and "pair" mean below? I don't understand the indexing

.. include:: autodoc_abbr_options_c.rst

.. index::
   single: DLPNO-MP2
   pair: MP2; density-fitting

.. _`sec:dlpnomp2`:

DLPNO-MP2: Domain-Based Local Pair Natural Orbital MP2
======================================================

.. codeauthor:: Zach Glick 
.. sectionauthor:: Zach Glick

*Module:* :ref:`Keywords <apdx:dlpnomp2>`, :ref:`PSI Variables <apdx:dlpnomp2_psivar>`, :source:`DLPNOMP2 <psi4/src/psi4/dlpnomp2>`

Introduction
------------

The steep polynomial scaling (in both time and memory) of post-HF dynamic
correlation methods prohibits calculations on large systems, even for 
efficient codes like |Psi4four|'s DFMP2. This poor scaling is in part due to
the use of canonical HF orbitals, which are entirely delocalized across the 
molecule, to compute the correlation energy. Canonical orbitals are usually
used for their mathematical convencience. Another possible choice is localized
orbitals. any two orbitals localized to separate regions of a molecule
can be treated as non-interacting to a good approximation. Thus, when working
with localized orbitals, the number of interacting orbital pairs (and triples, 
quadruples, etc.) scales linearly with system size. If carefully implemented, 
programs that exploit this sparsity can be made to scale linearly (or else with 
lower order than their canonical counterparts) at the cost of of modest, 
controllable errors. This is the defining insight of DLPNO-MP2 and all "local
correlation" methods.

The DLPNO-MP2 code is a linear-scaling alternative to the :ref:`DF-MP2 <sec:dfmp2>`
code, and is intended for use with large systems for which DF-MP2 is intractable. 
When running DLPNO-MP2 with default settings, approximately 99.9% of the DF-MP2 
correlation energy is recovered. The general outline of the method is as follows:

    (1) Localize the canonical HF orbitals: 
        (a) Perform Foster-Boys localization on the active occupied MOs
        (b) Construct Projected AOs from the virtual MOs
    (3) Calculate three-index integrals in (sparse) LMO/PAO basis
    (3) Perform local density fitting
    (4) Transform local virtual space from PAOs to truncated PNOs 
    (5) Solve the iterative local MP2 equations in the LMO/PNO basis

An example input file is::

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
   
   energy('dlpno-mp2')

The only difference between this input and a DF-MP2 input is the ``energy('dlpno-mp2')``
call to :py:func:`~psi4.energy`. Note that the water molecule in this example is far too 
small for DLPNO-MP2 to be of any benefit relative to DF-MP2.


The theory, thresholds, and best practices associated with DLPNO-MP2 are presented below. 

.. index::
   pair: DLPNO-MP2; theory

Theory
------

See :ref:`DF-MP2 <sec:dfmp2>` for background on the theory of (non-local) density-fitted MP2.

TODO

Recommendations
---------------

Some practical notes on running the code:

* DLPNO-MP2 should be used for large calculations that cannot reasonable be performed
  with DF-MP2 (i.e. > X basis functions).

* The greater the spatial sparsity of the system, the smaller the pair domains and
  consequently the faster the calculation. A linear alkane will run much faster than 
  a globular protein, all else constant.

* Similar to molecular sparsity, the sparsity of the orbital basis affects the runtime.
  Using diffuse functions increases the size of the pair domains and will slow down
  the calculations.

* All aspects of DLPNO-MP2 run in core; no disk is required. An upside of this fact
  is that the code exhibits very good intra-node parallelism (since no I/O occurs).
  The amount of memory needed scales asymptotically linear with system size.

* DLPNO-MP2 is not symmetry aware. This should not be a concern for large systems in
  which symmetry is seldom present.

* The current code only works with RHF references.

* As with DF-MP2, freezing core orbitals is recommended for efficiency.
