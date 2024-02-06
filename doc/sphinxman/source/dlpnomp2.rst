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
   single: DLPNO-MP2

.. _`sec:dlpnomp2`:

DLPNO-MP2: Domain-Based Local Pair Natural Orbital MP2
======================================================

.. codeauthor:: Zach Glick 
.. sectionauthor:: Zach Glick

*Module:* :ref:`Keywords <apdx:dlpno>`, :ref:`PSI Variables <apdx:dlpno_psivar>`, :source:`DLPNOMP2 <psi4/src/psi4/dlpno>`

Introduction
------------

The steep polynomial scaling (in both time and memory) of post-HF dynamic
correlation methods prohibits calculations on large systems, even for efficient
codes like |PSIfours| :ref:`DF-MP2 <sec:dfmp2>`. This poor scaling is in part
due to the use of canonical HF orbitals, which are entirely delocalized across
the molecule. Canonical orbitals are commonly used because of mathematical
convenience. Another possible choice is localized orbitals. Any two orbitals
localized to separate regions of a molecule can be treated as non-interacting
to a good approximation. Thus, when working with localized orbitals, the number
of interacting orbital pairs (and triples, quadruples, etc.) scales linearly
with system size. If carefully implemented, programs that exploit this sparsity
can be made to scale linearly (or else with lower order than their canonical
counterparts) at the cost of of modest, controllable errors. This is the
defining insight of DLPNO-MP2 and all related "local correlation" methods.

The DLPNO-MP2 code is a linear-scaling alternative to the :ref:`DF-MP2 <sec:dfmp2>`
code, and is intended for use with large systems for which DF-MP2 is intractable.
When running DLPNO-MP2 with default settings, approximately 99.9% of the DF-MP2 
correlation energy is recovered. The general outline of the method is as follows:

    (1) Localize the active occupied MOs (with the Foster-Boys method)
    (2) Construct projected AOs (PAOs) from the virtual MOs
    (3) Calculate three-index integrals in the (sparse) LMO/PAO basis
    (4) Perform local density fitting to construct (sparse) exchange integrals
    (5) Transform local virtuals from PAOs to pair natural orbitals (PNOs), and truncate
    (6) Solve the iterative local MP2 equations in the LMO/PNO basis

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
   set pno_convergence normal
   
   energy('dlpno-mp2')

The main difference between this input and a DF-MP2 input is the ``energy('dlpno-mp2')``
call to :py:func:`~psi4.driver.energy`. The only other addition is the |dlpno__pno_convergence|
keyword, which determines the accuracy of the local approximations underlying 
the DLPNO-MP2 method. Note that the water molecule in this example is not large
enough for DLPNO-MP2 to be of any benefit relative to DF-MP2.

The theory of the DLPNO-MP2 method and practical recommendations for using the
code are presented below. 

.. index::
   pair: DLPNO-MP2; theory

Theory
------

See :ref:`DF-MP2 <sec:dfmp2>` for background on the theory of (non-local)
density-fitted MP2. |PSIfours| DLPNO-MP2 implementation is based on the 
manuscript by Pinski et al. [Pinski:2015:034108]_.

In DLPNO-MP2, as in all local MP2 methods, the second-order MBPT energy is determined 
variationally via the Hylleraas functional [Hylleraas:1930:209]_:

.. math::
   :label: Hylleraas

   E^{(2)} = 2 \langle \Psi_{0}^{(0)} | \hat{H} - E_{0}^{(0)} | \Psi_{0}^{(1)} \rangle - \langle \Psi_{0}^{(1)} | \hat{H}^{(0)} - E_{0}^{(0)} | \Psi_{0}^{(1)} \rangle = \min_{| \Psi_{0}^{(1)} \rangle}.

Determining the optimal :math:`| \Psi_{0}^{(1)} \rangle` entails iteratively
minimizing the following residual [Pulay:1986:357]_:

.. math::
   :label: Residual

   R_{ij}^{ab} = (ia|jb) + (\epsilon_a + \epsilon_b - f_{ii} - f_{jj})t_{ij}^{ab} - \sum_{k \ne j} f_{ik} \sum_{c,d} S_{ac}t_{kj}^{c,d}S_{db} - \sum_{k \ne i} f_{kj} \sum_{cd} S_{ac}t_{ik}^{cd}S_{db} 

where ``i``, ``j``, and ``k`` are (not necessarily canonical) occupied orbitals, ``a``,
``b``, ``c``, and ``d`` are virtual orbitals, :math:`f_{ij}` are fock matrix elements,
:math:`S_{ab}` are orbital overlaps, and finally :math:`t_{ij}^{ab}` are the MP2
amplitudes to be solved for. Virtual orbitals may be different for each pair
of occupied orbitals. For a given occupied orbital pair ``ij``, all virtuals are
orthogonal and canonical, but virtuals belonging to different pair domains
may not be orthogonal.

The following expression is used to evaluate the energy of a given set of amplitudes:

.. math::
   :label: Energy

   E^{(2)} &= \sum_{i,j} e_{ij}, \\
   e_{ij} &= \sum_{a, b}((ia|jb) + R_{ij}^{ab})(2t_{ij}^{ab} - t_{ij}^{ba}).

The error in :math:`E^{(2)}` scales quadratically with the error in the amplitudes.

No local approximations have been made so far, and this iterative approach can
be used to exactly determine :math:`E^{(2)}` with :math:`{\cal O}(N^5)` cost.
In DLPNO-MP2, the first local approximation is to screen distant, non-interacting
orbital pairs ``ij``. Orbital pairs are screened if below both an overlap criteria:

.. math::
   :label: Differential Overlap Integral

   DOI_{ij} \equiv \sqrt{\int d\mathbf{r} | \chi_{i}(\mathbf{r}) | ^{2} | \chi_{j}(\mathbf{r}) | ^{2}}, 

and a pair energy estimate: 

.. math::
   :label: Dipole Approximation

   e_{ij}^{approx} = -\frac{4}{R^{6}} \sum_{a_{i} \in [i],b_{j} \in [j]} \frac{ (2 \langle i | \mathbf{r} | a_{i} \rangle \langle j | \mathbf{r} | b_{j} \rangle)^{2}}{\epsilon_{a_{i}} + \epsilon_{b_{j}} - f_{ii} - f_{jj}},

in which small domains of virtual orbitals are used for each local MO. As a
result, an asymptotically linear number of ``ij`` pairs enter the local MP2
equations, and the approximate pair energy of neglected pairs is added to
the final energy.

The second major local approximation in DLPNO-MP2 is the truncation of the virtual
space. Initially, exchange integrals are calculated in the LMO/PAO basis using the
standard density fitting approach:

.. math::
   :label: Exchange

   (ia|jb) = \sum_{K,L \in [ij]} (ia|K)[\mathbf{J}^{-1}]_{KL}(L|jb)

This is done with linear scaling effort by exploiting the locality of the LMOs, PAOs,
and auxiliary basis functions.
Solving the iterative local MP2 equations in the LMO/PAO basis requires large PAO
domains to achieve reasonable accuracy. Instead, the virtual space is transformed
into the much more compact pair natural orbital representation. The (approximate)
PNOs diagonalize the virtual-virtual block of the (approximate) MP2 density matrix:

.. math::
   :label: PNOs

   D_{ij}^{ab} = \frac{1}{1 + \delta_{ij}}[\tilde{t}_{ij}^{\dagger}t_{ij} + \tilde{t}_{ij}t_{ij}^{\dagger}]^{ab}

which is constructed from semicanonical amplitudes:

.. math::
   :label: Amplitudes

   t_{ij}^{ab}         &= - \frac{(iajb)}{\epsilon_{a} + \epsilon_{b} - \epsilon_{i} - \epsilon_{j}}, \\
   \tilde{t}_{ij}^{ab} &= 2t_{ij}^{ab} - t_{ij}^{ba}.

PNOs with small occupation numbers are discarded, and the local MP2 equations are
solved in the LMO/PNO basis.

Recommendations
---------------

Some practical notes on running the code:

* DLPNO-MP2 is not a drop-in replacement for DF-MP2. Instead, it should be used for
  large calculations that cannot reasonably be performed with DF-MP2. The crossover
  point between DF-MP2 and DLPNO-MP2 depends on details of both the calculation and
  the hardware, but can be as low as 2,000 basis functions.

* The accuracy of DLPNO-MP2 (relative to DF-MP2) can be controlled with the
  |dlpno__pno_convergence| keyword according to recommendation by Liakos et al.
  [Liakos:2015:1525]_. For non-covalent interactions ``TIGHT`` is highly recommended. 
  
* The greater the spatial sparsity of a molecular system, the smaller the pair
  domains and consequently the faster the calculation. DLPNO-MP2 is much faster
  for linear alkanes than for globular proteins, all else constant.

* Similar to molecular sparsity, the sparsity of the orbital basis affects runtime.
  Diffuse functions increase the size of the pair domains and therefore lead to 
  longer calculations.

* All aspects of DLPNO-MP2 run in core; no disk is required. As a result, the
  code exhibits very good intra-node parallelism, and benefits from many threads.
  The amount of memory needed scales linearly with system size.

* DLPNO-MP2 is not symmetry aware. This should not be a concern for large systems in
  which symmetry is seldom present.

* As with DF-MP2, freezing core orbitals (by setting |globals__freeze_core|
  to ``True``) is recommended for efficiency. In DLPNO methods, this is also
  recommended for accuracy, since core excitations are known to exhibit
  greater errors relative to valence excitations.

* At the moment, the DLPNO-MP2 code is only compatible with with RHF references.
