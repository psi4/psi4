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

.. _`sec:introduction`:

============
Introduction
============

Overview
========

|PSIfour| provides a wide variety of quantum chemical methods using
state-of-the-art numerical methods and algorithms.  Several parts of
the code feature shared-memory parallelization to run efficiently on
multi-core machines (see Sec. :ref:`sec:threading`).
An advanced parser written in Python allows the user
input to have a very simple style for routine computations, but it can also
automate very complex tasks with ease. 

|PSIfour| is, in many ways, a whole new package compared to Psi3.
While some libraries and modules remain the same, the majority of the code has
been rewritten from scratch based on a powerful set of new libraries written
in C++.  A totally new Python front-end makes |PSIfour| incredibly user-friendly
and automates many common tasks such as basis set extrapolation, composite
methods, running the same computation on every molecule in a test set, etc.
Density-functional theory, absent in Psi3, is quite efficient
in |PSIfour|, with many functionals available.  Density fitting is ubiquitous in
|PSIfour|, leading to some of the most efficient MP2 and CCSD(T) code available.
|PSIfour| also introduces extensive,
powerful features for energy component analysis of non-covalent interactions
via symmetry-adapted perturbation theory.  Orbital-optimized versions of
perturbation theory and coupled-cluster methods, and their analytic gradients,
have also been added.  Through external libraries, |PSIfour| gains access to implicit
solvent (PCM) capabilities, density-matrix renormalization group CI, effective
fragment potentials, Grimme dispersion corrections, and high-order
coupled-cluster theory.

In this section, we provide an overview of some of the features of
|PSIfour| along with the prerequisite steps for running calculations.
Sec. :ref:`Tutorial <sec:tutorial>` provides a brief tutorial to help new users
get started.  Section :ref:`Psithon <sec:psithonInput>` offers further details into the
structure of |PSIfour| input files and how Python can be mixed with
quantum chemistry directives in |PSIfour|. Section :ref:`Psithon Functions <sec:psithonFunc>`
provides more detail on the Python functions provided by |PSIfour|
and discusses some of the higher-level functions such as counterpoise
correction, complete-basis-set extrapolation, and running computations
on an entire database of molecules at a time.  Later sections deal with
the different types of computations which can be done using |PSIfour|
(*e.g.*, Hartree |--| Fock, MP2, coupled-cluster) and general procedures
such as geometry optimization and vibrational frequency analysis.
The :ref:`Appendices <sec:appendices>` include a complete description of all possible input
keywords for each module, as well as tables of available basis sets and
a listing of the sample input files available under :source:`samples`.
The user is urged to examine this directory of sample inputs, as
most common types of computations are represented there.
For the latest |PSIfour| documentation, check 
`www.psicode.org <http://www.psicode.org>`_.

Citing |PSIfour|
================

Overall |PSIfour| Package
^^^^^^^^^^^^^^^^^^^^^^^^^

The following citation should be used in any publication utilizing the
|PSIfour| program package:

* "Psi4: An open-source *ab initio* electronic structure program,"
  J. M. Turney, A. C. Simmonett, R. M. Parrish, E. G. Hohenstein, F.
  Evangelista, J. T. Fermann, B. J. Mintz, L. A. Burns, J. J. Wilke, M. L.
  Abrams, N. J.  Russ, M. L. Leininger, C. L. Janssen, E. T. Seidl, W. D.
  Allen, H. F.  Schaefer, R. A. King, E. F. Valeev, C. D. Sherrill, and T.
  D. Crawford, *WIREs Comput. Mol. Sci.* **2**, 556 (2012).
  (doi: `10.1002/wcms.93 <http://dx.doi.org/10.1002/wcms.93>`_).


Depending on the particular modules used, the user may also wish to
cite some of the following references for theoretical, algorithmic,
or implementation contributions specific to |PSIfour| (in addition to
appropriate references for the underlying theory, which are not necessarily
included in the list below).

Density Cumulant Functional Theory (DCFT)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _`intro:dcftcitations`:

|PSIfour| features several formulations of newly-developed density cumulant
functional theory (DCFT). The theory and benchmark of this theory are
discussed in the following papers:

DC-06 (also known as DCFT-06):

* "Density Cumulant Functional Theory: First Implementation and
  Benchmark Results for the DCFT-06 Model," A. C. Simmonett,
  J. J. Wilke, H. F. Schaefer, and W. Kutzelnigg, *J. Chem. Phys.*
  **133**, 174122 (2010).
  (doi: `10.1063/1.3503657 <http://dx.doi.org/10.1063/1.3503657>`_).

* "Analytic gradients for density cumulant functional theory: The 
  DCFT-06 model," A. Yu. Sokolov, J. J. Wilke, A. C. Simmonett, 
  and H. F. Schaefer, *J. Chem. Phys.* **137**, 054105 (2012).
  (doi: `10.1063/1.4739423 <http://dx.doi.org/10.1063/1.4739423>`_).

DC-12:

* "Density cumulant functional theory: The DC-12 method, an improved 
  description of the one-particle density matrix," A. Yu. Sokolov, 
  A. C. Simmonett, and H. F. Schaefer, *J. Chem. Phys.*  **138**, 024107 
  (2013).
  (doi: `10.1063/1.4773580 <http://dx.doi.org/10.1063/1.4773580>`_).

ODC-06 and ODC-12:

* "Orbital-optimized density cumulant functional theory," A. Yu. Sokolov, and
  H. F. Schaefer, *J. Chem. Phys.*  **139**, 204110 (2013). 
  (doi: `10.1063/1.4833138 <http://dx.doi.org/10.1063/1.4833138>`_).

ODC-13:

* "Density cumulant functional theory from a unitary transformation:
  N-representability, three-particle correlation effects, and application
  to O4+," A. Yu. Sokolov, H. F. Schaefer, and W. Kutzelnigg,
  *J. Chem. Phys.*  **141**, 074111 (2014). 
  (doi: `10.1063/1.4892946 <http://dx.doi.org/10.1063/1.4892946>`_).

Configuration Interaction (CI)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PSI has a highly optimized code for full configuration interaction
and highly correlated configuration interaction, as described in

* "The Configuration Interaction Method: Advances in Highly 
  Correlated Approaches," C. D. Sherrill and H. F. Schaefer, in
  *Adv. Quantum Chem.*, vol. 34, P.-O. L\ |o_dots|\ wdin, Ed.
  (Academic Press, New York, 1999), pp. 143-269.
  (doi: `10.1016/S0065-3276(08)60532-8
  <http://dx.doi.org/10.1016/S0065-3276(08)60532-8>`_).

Coupled Cluster (CC)
^^^^^^^^^^^^^^^^^^^^

A general discussion of coupled cluster theory is given in

* "An Introduction to Coupled Cluster Theory for Computational
  Chemists," T. D. Crawford and H. F. Schaefer, *Rev. Comp. Chem.* 
  **14**, 33-136 (2000).
  (doi: `10.1002/9780470125915.ch2
  <http://dx.doi.org/10.1002/9780470125915.ch2>`_).

Implementation of frozen natural orbital (FNO) coupled cluster theory
in PSI and its performance for non-covalent interactions is discussed
in

* "Accurate Noncovalent Interaction Energies Using Truncated Basis Sets 
  Based on Frozen Natural Orbitals," A. E. DePrince and C. D. Sherrill,
  *J. Chem. Theory Comput.* **9**, 293-299 (2013).
  (doi: `10.1021/ct300780u <http://dx.doi.org/10.1021/ct300780u>`_).

Implementation of density-fitted (DF) and Cholesky decomposition (CD)
coupled cluster in PSI, and its performance for non-covalent interactions
and reaction energies, is discussed in

* "Accuracy and Efficiency of Coupled-Cluster Theory Using
  Density Fitting / Cholesky Decomposition, Frozen Natural Orbitals,
  and a T1-Transformed Hamiltonian," A. E. DePrince and C. D. Sherrill,
  *J. Chem. Theory Comput.* **9**, 2687-2696 (2013).
  (doi: `10.1021/ct400250u <http://dx.doi.org/10.1021/ct400250u>`_).

Implementation of the asymmetric triples correction for the density-fitted 
and cholesky-decomposed coupled-cluster singles and doubles method

* "A noniterative asymmetric triple excitation correction for the density-fitted
  coupled-cluster singles and doubles method: Preliminary applications," 
  U. Bozkaya,   *J. Chem. Phys.* **144**, 144108 (2016).
  (doi: `10.1063/1.4945706 <http://dx.doi.org/10.1063/1.4945706>`_).

 
Mukherjee State-Specific Multi-Reference Coupled Cluster (Mk-MRCC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
|PSIfour| features production-level Mukherjee-style state-specific 
coupled-cluster theory, including perturbative triples and also associated
multi-reference perturbation theories.  The theory and |PSIfour| 
implementation of these methods is discussed in the following papers.

General Mk-MRCC

* "Coupling Term Derivation and General Implementation of
  State-Specific Multireference Coupled-Cluster Theories,"
  F. A. Evangelista, W. D. Allen, and H. F. Schaefer, 
  *J. Chem. Phys.* **127**, 024102 (2007).
  (doi: `10.1063/1.2743014 <http://dx.doi.org/10.1063/1.2743014>`_).

Mk-MRCCSD(T)

* "Perturbative Triples Corrections in State-Specific Multireference
  Coupled Cluster Theory,"
  F. A. Evangelista, E. Prochnow, J. Gauss, and H. F. Schaefer,
  *J. Chem. Phys.* **132**, 074107 (2010).
  (doi: `10.1063/1.3305335 <http://dx.doi.org/10.1063/1.3305335>`_).

Mk-MRCCSDT(-n)

* "Triple Excitations in State-Specific Multireference Coupled
  Cluster Theory: Application of Mk-MRCCSDT and Mk-MRCCSDT-n Methods to
  Model Systems," F. A. Evangelista, A. C. Simmonett, W. D. Allen,
  H. F. Schaefer, and J. Gauss, *J. Chem. Phys.* **128**, 124104
  (2008).
  (doi: `10.1063/1.2834927 <http://dx.doi.org/10.1063/1.2834927>`_).

Mk-MRPT2

* "A Companion Perturbation Theory for State-specific
  Multireference Coupled Cluster Methods,"
  F. A. Evangelista, A. C. Simmonett, H. F. Schaefer, D. Mukherjee, and
  W. D. Allen,
  *Phys. Chem. Chem. Phys.* **11**, 4728-4741 (2009).
  (doi: `10.1039/b822910d <http://dx.doi.org/10.1039/b822910d>`_).

Symmetry-Adapted Perturbation Theory (SAPT)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|PSIfour| features an extremely efficient code to perform wavefunction-based
Symmetry Adapted Perturbation Theory (SAPT). A good review article for this 
method is as follows:

* "Perturbation Theory Approach to Intermolecular Potential Energy
  Surfaces of van der Waals Complexes," B. Jeziorski, R. Moszynski,
  and K. Szalewicz, *Chem. Rev.* **94**, 1887-1930 (1994).   
  (doi: `10.1021/cr00031a008 <http://dx.doi.org/10.1021/cr00031a008>`_).

|PSIfour| benefits enormously from the introduction of density fitting (DF)
into SAPT. There are several SAPT truncations available in |PSIfour|. For 
guidance on which one to choose, see the SAPT section of the manual
and refer to the following systematic study:

* "Levels of  Symmetry Adapted Perturbation Theory (SAPT). I. Efficiency and
  Performance for Interaction Energies,'' T. M. Parker, L. A. Burns, R. M.
  Parrish, A. G. Ryno, and C. D. Sherrill, *J. Chem. Phys.* **140**, 
  094106 (2014).
  (doi: `10.1063/1.4867135 <http://dx.doi.org/10.1063/1.4867135>`_).

The theory and implementation of DF-SAPT is discussed 
in the following papers for various levels of SAPT.

DF-SAPT0

* "Large-scale Symmetry-adapted Perturbation Theory Computations via
  Density Fitting and Laplace Transformation Techniques: Investigating the
  Fundamental Forces of DNA-Intercalator Interactions," E. G. Hohenstein,
  R. M. Parrish, C. D. Sherrill, J. M. Turney, and H. F. Schaefer, *J.
  Chem. Phys.* **135**, 174017 (2011).
  (doi: `10.1063/1.3656681 <http://dx.doi.org/10.1063/1.3656681>`_).

* "Density Fitting and Cholesky Decomposition Approximations
  in Symmetry-Adapted Perturbation Theory: Implementation and Application
  to Probe the Nature of :math:`\pi - \pi` Interactions in Linear Acenes,"
  E. G. Hohenstein and C. D. Sherrill, *J. Chem. Phys.* **132**,
  184111 (2010).
  (doi: `10.1063/1.3426316 <http://dx.doi.org/10.1063/1.3426316>`_).

SAPT2

* "Density Fitting of Intramonomer Correlation Effects in
  Symmetry-Adapted Perturbation Theory,"
  E. G. Hohenstein and C. D. Sherrill, *J. Chem. Phys.* **133**,
  014101 (2010).
  (doi: `10.1063/1.3451077 <http://dx.doi.org/10.1063/1.3451077>`_).

SAPT2+, SAPT2+(3), SAPT2+3

* "Wavefunction Methods for Noncovalent Interactions," E. G.
  Hohenstein and C. D. Sherrill, *WIREs: Comput. Mol. Sci.* **2**,
  304-326 (2012).
  (doi: `10.1002/wcms.84 <http://dx.doi.org/10.1002/wcms.84>`_).

* "Density Fitting of Intramonomer Correlation Effects in
  Symmetry-Adapted Perturbation Theory,"
  E. G. Hohenstein and C. D. Sherrill, *J. Chem. Phys.* **133**,
  014101 (2010).
  (doi: `10.1063/1.3451077 <http://dx.doi.org/10.1063/1.3451077>`_).

* "Efficient Evaluation of Triple Excitations in Symmetry-Adapted
  Perturbation Theory via MP2 Natural Orbitals," E. G. Hohenstein
  and C. D. Sherrill, *J. Chem. Phys.* **133**, 104107 (2010).
  (doi: `10.1063/1.3479400 <http://dx.doi.org/10.1063/1.3479400>`_).


SAPT2+(CCD), SAPT2+(3)(CCD), and SAPT2+3(CCD)

* "Tractability Gains in Symmetry-Adapted Perturbation Theory Including
  Coupled Double Excitations: CCD+ST(CCD) Dispersion with Natural Orbital
  Truncations,'' R. M. Parrish, E. G. Hohenstein, and C. D. Sherrill, 
  *J. Chem. Phys.* **139**, 174102 (2013).
  (doi: `10.1063/1.4826520 <http://dx.doi.org/10.1063/1.4826520>`_).

* "Wavefunction Methods for Noncovalent Interactions," E. G.
  Hohenstein and C. D. Sherrill, *WIREs: Comput. Mol. Sci.* **2**,
  304-326 (2012).
  (doi: `10.1002/wcms.84 <http://dx.doi.org/10.1002/wcms.84>`_).

* "Density Fitting of Intramonomer Correlation Effects in
  Symmetry-Adapted Perturbation Theory,"
  E. G. Hohenstein and C. D. Sherrill, *J. Chem. Phys.* **133**,
  014101 (2010).
  (doi: `10.1063/1.3451077 <http://dx.doi.org/10.1063/1.3451077>`_).

F/I-SAPT

* "Chemical Assignment of Symmetry-Adapted Perturbation Theory Interaction
  Energy Components: The Functional-Group SAPT Partition,"
  R. M. Parrish, T. M. Parker, and C. D. Sherrill,
  *J. Chem. Theory Comput.* **10**, 4417 (2014).
  (doi: `10.1021/ct500724p <http://dx.doi.org/10.1021/ct500724p>`_).

* "Communication: Practical Intramolecular Symmetry Adapted Perturbation Theory
  via Hartree-Fock Embedding,"
  R. M. Parrish, J. F. Gonthier, C. Corminboeuf, and C. D. Sherrill,
  *J. Chem. Phys.* **143**, 051103 (2015).
  (doi: `10.1063/1.4927575 <http://dx.doi.org/10.1063/1.4927575>`_)

Orbital-Optimized Post-Hartree |--| Fock Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Orbital-optimized second-order perturbation theory (OMP2)

* "Quadratically convergent algorithm for orbital optimization in the 
  orbital-optimized coupled-cluster doubles method and in orbital-optimized 
  second-order |MollerPlesset| perturbation theory,"
  U. Bozkaya, J. M. Turney, Y. Yamaguchi, H. F. Schaefer, and C. D. Sherrill,
  *J. Chem. Phys.* **135**, 104103 (2011).
  (doi: `10.1063/1.3631129 <http://dx.doi.org/10.1063/1.3631129>`_).

* "Analytic energy gradients for the orbital-optimized second-order 
  |MollerPlesset| perturbation theory," U. Bozkaya and
  C. D. Sherrill, *J. Chem. Phys.* **138**, 184103 (2013).
  (doi: `10.1063/1.4803662 <http://dx.doi.org/10.1063/1.4803662>`_).

* "Orbital-Optimized Second-Order Perturbation Theory with Density-Fitting
  and Cholesky Decomposition Approximations: An Efficient Implementation,"
  U. Bozkaya,   *J. Chem. Theory Comput.* **10**, 2371 (2014).
  (doi: `10.1021/ct500231c <http://dx.doi.org/10.1021/ct500231c>`_).

Orbital-optimized third-order perturbation theory (OMP3)

* "Orbital-Optimized Third-Order |MollerPlesset| Perturbation
  Theory and Its Spin-Component and Spin-Opposite Scaled Variants: Application 
  to Symmetry Breaking Problems," U. Bozkaya,
  *J. Chem. Phys.* **135**, 224103 (2011).
  (doi: `10.1063/1.3665134 <http://dx.doi.org/10.1063/1.3665134>`_).

* "Assessment of Orbital-Optimized Third-Order |MollerPlesset|
  Perturbation Theory and Its Spin-Component and Spin-Opposite Scaled Variants 
  for Thermochemistry and Kinetics," E. Soydas and U. Bozkaya,  
  *J. Chem. Theory Comput.* **9**, 1452 (2013).
  (doi: `10.1021/ct301078q <http://dx.doi.org/10.1021/ct301078q>`_).

* "Analytic energy gradients for the orbital-optimized third-order |MollerPlesset|
  Perturbation Theory," U. Bozkaya,  
  *J. Chem. Phys.* **139**, 104116 (2013).
  (doi: `10.1063/1.4820877 <http://dx.doi.org/10.1063/1.4820877>`_).

Orbital-optimized linearized coupled-cluster doubles method (OLCCD)

* "Orbital-optimized coupled-electron pair theory and its analytic gradients: 
  Accurate equilibrium geometries, harmonic vibrational frequencies, and hydrogen transfer 
  reactions," U. Bozkaya and C. D. Sherrill,
  *J. Chem. Phys.* **139**, 054104 (2013).
  (doi: `10.1063/1.4816628 <http://dx.doi.org/10.1063/1.4816628>`_).

Orbital-optimized MP2.5 (OMP2.5)

* "Orbital-optimized MP2.5 and its analytic gradients: Approaching CCSD(T)
  quality for noncovalent interactions," U. Bozkaya and C. D. Sherrill, 
  *J. Chem. Phys.* **141**, 204105 (2014).
  (doi: `10.1063/1.4902226 <http://dx.doi.org/10.1063/1.4902226>`_).

Extended Koopmans' Theorem

* "The extended Koopmans' theorem for orbital-optimized methods: Accurate 
  computation of ionization potentials," U. Bozkaya,  *J. Chem. Phys.* 
  **139**, 154105 (2013).
  (doi: `10.1063/1.4825041 <http://dx.doi.org/10.1063/1.4825041>`_).

* "Accurate Electron Affinities from the Extended Koopmans' Theorem Based on Orbital-Optimized Methods,"
  U. Bozkaya,   *J. Chem. Theory Comput.* **10**, 2041 (2014).
  (doi: `10.1021/ct500186j <http://dx.doi.org/10.1021/ct500186j>`_).

Density-Fitted and Cholesky-Decomposed Orbital-optimized second-order perturbation theory (DF-OMP2)

* "Orbital-Optimized Second-Order Perturbation Theory with Density-Fitting 
  and Cholesky Decomposition Approximations: An Efficient Implementation," 
  U. Bozkaya,   *J. Chem. Theory Comput.* **10**, 2371 (2014).
  (doi: `10.1021/ct500231c <http://dx.doi.org/10.1021/ct500231c>`_).

* "Analytic Energy Gradients and Spin Multiplicities for Orbital-Optimized 
  Second-Order Perturbation Theory with Density-Fitting Approximation: An 
  Efficient Implementation," U. Bozkaya, *J. Chem. Theory Comput.* **10**, 4389 (2014).
  (doi: `10.1021/ct500634s <http://dx.doi.org/10.1021/ct500634s>`_).

Density-Fitted and Cholesky-Decomposed Orbital-optimized MP3 and MP2.5 (DF-OMP3 and DF-OMP2.5)

* "Orbital-Optimized MP3 and MP2.5 with Density-Fitting 
  and Cholesky Decomposition Approximations," 
  U. Bozkaya,   *J. Chem. Theory Comput.* **12**, 1179 (2016).
  (doi: `10.1021/acs.jctc.5b01128 <http://dx.doi.org/10.1021/acs.jctc.5b01128>`_).

Density-Fitted and Cholesky-Decomposed Orbital-Optimized Linearized Coupled-Cluster Doubles Method (DF-OLCCD)

* "Orbital-optimized linearized coupled-cluster doubles with density-fitting 
  and Cholesky decomposition approximations: an efficient implementation," 
  U. Bozkaya,   *Phys. Chem. Chem. Phys.* **18**, 11362 (2016).
  (doi: `10.1039/c6cp00164e <http://dx.doi.org/10.1039/c6cp00164e>`_).


Second-Order Algebraic-Diagrammatic Construction [ADC(2)]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

General ADC(2) theory

* "Intermediate state representation approach to physical properties of 
  electronically excited molecules,"
  J. Schirmer, and A. B. Trofimov, *J. Chem. Phys.* **120**,
  11449-11464 (2004).
  (doi: `10.1063/1.1752875 <http://dx.doi.org/10.1063/1.1752875>`_).

Theory of "Partially-renormalized" CIS(D) and ADC(2) [PR-CIS(D) and PR-ADC(2)]
and their implementation in |PSIfour|

* "Excited State Calculation for Free-Base and Metalloporphyrins with
  the Partially Renormalized Polarization Propagator Approach,"
  M. Saitow and Y. Mochizuki, *Chem. Phys. Lett.* **525**, 144-149
  (2012).  
  (doi: `10.1016/j.cplett.2011.12.063 
  <http://dx.doi.org/10.1016/j.cplett.2011.12.063>`_).

Density Matrix Renormalization Group (DMRG)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* "CheMPS2: a free open-source spin-adapted implementation of the density 
  matrix renormalization group for ab initio quantum chemistry,"
  S. Wouters, W. Poelmans, P. W. Ayers and D. Van Neck,
  *Comput. Phys. Commun.* **185** (6), 1501-1514 (2014).
  (doi: `10.1016/j.cpc.2014.01.019 <http://dx.doi.org/10.1016/j.cpc.2014.01.019>`_).

* "The density matrix renormalization group for ab initio quantum chemistry,"
  S. Wouters and D. Van Neck, *Eur. Phys. J. D* **68** (9), 272 (2014).
  (doi: `10.1140/epjd/e2014-50500-1 <http://dx.doi.org/10.1140/epjd/e2014-50500-1>`_).

Scalar Relativistic Corrections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

General theory for the exact two-component approach (X2C)

* "Analytic energy gradients for the spin-free exact two-component theory
  using an exact block diagonalization for the one-electron Dirac
  Hamiltonian,"
  L. Cheng and J. Gauss, *J. Chem. Phys.* **135**, 084114 (2011).
  (doi: `10.1063/1.3624397 <http://dx.doi.org/10.1063/1.3624397>`_).

Implementation within Psi4

* "Predicting Near Edge X-ray Absorption Spectra with the Spin-Free
  Exact-Two-Component Hamiltonian and Orthogonality Constrained Density
  Functional Theory,"
  P. Verma, W. D. Derricotte and F. A. Evangelista,
  *J. Chem. Theory Comput.* (2015).
  (doi: `10.1021/acs.jctc.5b00817 <http://dx.doi.org/10.1021/acs.jctc.5b00817>`_).

.. index:: architectures
.. index:: compilers

Supported Architectures
=======================

The majority of |PSIfour| was developed on Mac and Linux machines; in
principle, it should work on any Unix system. The latest version of the
|PSIfour| program package may be obtained at `psicode.org
<http://psicode.org>`_. The package is available as a binary for Linux,
macOS, or Windows Subsystem for Linux (:ref:`Installing from Binary
<sec:conda>`) or as source code (zipped archive or git repository from
`https://github.com/psi4/psi4 <http://www.github.com/psi4/psi4>`_).

|PSIfour| has been successfully compiled using Intel, GCC, and Clang
compilers. :ref:`Compiler requirements <faq:approvedcxx>` are primarily
C++11 compliance (*i.e.*, GCC version 4.9 or above).
For some architectures, a :ref:`precompiled binary
<sec:conda>` is available. See :ref:`Compiling and Installing
<sec:installFile>` for details.


Capabilities
============

|PSIfour| can perform *ab initio* computations employing basis
sets of contracted Gaussian-type functions of virtually arbitrary
orbital quantum number. Many parts of |PSIfour| can recognize and
exploit the largest Abelian subgroup of the molecular point group.
Table :ref:`Methods <table:methods>` displays the range of theoretical methods
available in |PSIfour|.
For more details, see Tables :ref:`Energy <table:energy_gen>`, 
:ref:`Energy (DFT) <table:energy_dft>`, :ref:`Energy (MRCC) <table:energy_mrcc>`,
:ref:`Energy (CFOUR) <table:energy_cfour>`, :ref:`Gradient <table:grad_gen>`, 
:ref:`Gradient (CFOUR) <table:grad_cfour>`, and :ref:`Frequency <table:freq_gen>`.

.. _`table:methods`:

.. table:: Summary of theoretical methods available in |PSIfour|

    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | Method                  | Energy                             | Gradient                      | Parallelism                 | Variants               |
    +                         +----------------+-------------------+--------------+----------------+                             +-----------+------------+
    |                         | Reference      | Type              | Reference    | Type           |                             | OO [#f9]_ | FNO [#f1]_ |
    +=========================+================+===================+==============+================+=============================+===========+============+
    | HF, HF-3c               | RHF/UHF/ROHF   | CONV/DF/CD        | RHF/UHF/ROHF | CONV/DF        | threaded                    |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | DFT, PBEh-3c            | RKS/UKS        | CONV/DF/CD [#f7]_ | RKS/UKS      | DF [#f4]_      | threaded                    |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | EFP [#f5]_              | RHF            |                   | ---          | ---            |                             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | DCFT                    | RHF/UHF        | CONV/DF           | RHF/UHF      | CONV           | partially threaded          |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | MP2 [#f10]_             | RHF/UHF/ROHF   | CONV/DF/CD        | RHF/UHF      | CONV/DF        | threaded [#f3]_             | E/G       |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | MP3 [#f10]_             | RHF/UHF        | CONV/DF/CD        | RHF/UHF      | CONV/DF        | threaded [#f3]_             | E/G       | E          |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | MP2.5 [#f10]_           | RHF/UHF        | CONV/DF/CD        | RHF/UHF      | CONV/DF        | threaded [#f3]_             | E/G       |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | MP4 [#f10]_             | RHF            | CONV              | ---          | ---            | threaded [#f3]_             |           | E          |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | MP(n)                   | RHF/ROHF       | CONV              | ---          | ---            | partially threaded          |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | ZAPT(n)                 | RHF/ROHF       | CONV              | ---          | ---            | partially threaded          |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | LCCD [#f10]_            | RHF/UHF        | CONV/DF/CD        | RHF/UHF      | CONV/DF        | threaded [#f3]_             | E/G       | E          |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | LCCSD, CEPA(0) [#f10]_  | RHF            | CONV              | ---          | ---            | threaded [#f3]_             |           | E          |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | CEPA(n), n=0,1,3        | RHF            | CONV              | ---          | ---            | threaded [#f3]_             |           | E          |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | CCD                     | RHF            | DF/CD             | RHF          | DF             | threaded [#f3]_             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | CC2                     | RHF/UHF/ROHF   | CONV              | RHF          | CONV            | threaded [#f3]_             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | CCSD [#f10]_            | RHF/UHF/ROHF   | CONV/DF/CD [#f8]_ | RHF/UHF/ROHF | CONV/DF [#f8]_ | threaded [#f3]_             |           | E [#f2]_   |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | CCSD(T) [#f10]_         | RHF/UHF/ROHF   | CONV/DF/CD [#f8]_ | RHF/UHF      | CONV           | threaded (pthreads) [#f3]_  |           | E [#f2]_   |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | CCSD(AT) [#f10]_        | RHF            | CONV/DF/CD        | ---          | ---            | threaded [#f3]_             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | CC3                     | RHF/UHF/ROHF   | CONV              | ---          | ---            | threaded (pthreads)         |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | ACPF/AQCC               | RHF            | CONV              | ---          | ---            | threaded [#f3]_             |           | E          |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | CISD [#f10]_            | RHF/ROHF       | CONV              | ---          | ---            | partially threaded          |           | E          |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | QCISD                   | RHF            | CONV              | ---          | ---            | threaded [#f3]_             |           | E          |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | QCISD(T)                | RHF            | CONV              | ---          | ---            | threaded [#f3]_             |           | E          |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | CI(n)                   | RHF/ROHF       | CONV              | ---          | ---            | partially threaded          |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | FCI                     | RHF/ROHF       | CONV              | ---          | ---            | partially threaded          |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | Mk-MRPT2                | RHF/ROHF/TCSCF | CONV              | ---          | ---            | threaded [#f3]_             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | Mk-MRCCSD               | RHF/ROHF/TCSCF | CONV              | ---          | ---            | threaded [#f3]_             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | Mk-MRCCSD(T)            | RHF/ROHF/TCSCF | CONV              | ---          | ---            | threaded [#f3]_             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | RAS-CI                  | RHF/ROHF       | CONV              | ---          | ---            | partially threaded          |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | CASSCF, RASSCF          | RHF/ROHF       | CONV/DF           | ---          | ---            | partially threaded          |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | SAPT                    | RHF            | CONV              | ---          | ---            | threaded                    |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | ADC(2)                  | RHF            | CONV              | ---          | ---            | threaded [#f3]_             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | EOM-CC2                 | RHF            | CONV              | ---          | ---            | threaded [#f3]_             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | EOM-CCSD                | RHF/UHF/ROHF   | CONV              | RHF/UHF/ROHF | CONV           | threaded [#f3]_             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | EOM-CC3                 | RHF/UHF/ROHF   | CONV              | ---          | ---            | threaded [#f3]_             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | DMRG-CI                 | Y              | CONV              | ---          | ---            |                             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+
    | DMRG-SCF                | Y              | CONV              | ---          | ---            |                             |           |            |
    +-------------------------+----------------+-------------------+--------------+----------------+-----------------------------+-----------+------------+

Geometry optimization can be performed using either analytic gradients
or energy points. Likewise, vibrational frequencies can be 
computed by analytic second derivatives, by finite
differences of analytic gradients, or by finite differences of energies.
|PSIfour| can also compute an extensive list of one-electron properties.

.. index::
   single: contact
   single: bugs

Technical Support
=================

The |PSIfour| package is distributed for free and without any guarantee of
reliability, accuracy, or suitability for any particular purpose. No
obligation to provide technical support is expressed or implied. As time
allows, the developers will attempt to answer inquiries on the `forum
<http://forum.psicode.org>`_ or `GitHub
<https://github.com/psi4/psi4/issues/new>`_. For bug reports,
specific and detailed information, with example inputs, would be
appreciated.

Where-to-post summary:[#f6]_

* How do I? -- `ask the forum <http://forum.psicode.org>`_

* I got this error, why? -- `ask the forum <http://forum.psicode.org>`_

* I got this error and I'm sure it's a bug -- `file a GitHub issue <https://github.com/psi4/psi4/issues/new>`_

* Can I open a discussion on this bit of code? -- `file a GitHub issue <https://github.com/psi4/psi4/issues/new>`_

* I have an idea/request and a plan -- `file a GitHub issue <https://github.com/psi4/psi4/issues/new>`_

* I have an idea/request -- `ask the forum <http://forum.psicode.org>`_

* Why do you? -- `ask the forum <http://forum.psicode.org>`_

* When will you? -- `ask the forum <http://forum.psicode.org>`_

* I have an experience that can improve the build documentation -- `inform the forum <http://forum.psicode.org>`_ or :source:`add to the documentation itself <doc/sphinxman/source>`

* Anything you want to share privately -- `crawdad@vt.edu <mailto:crawdad@vt.edu>`_ or `sherrill@gatech.edu <mailto:sherrill@gatech.edu>`_


.. rubric:: Footnotes

.. [#f1] Frozen natural orbital variant available. In particular, RHF available as CONV.
.. [#f2] Frozen natural orbital variant available. In particular, RHF available as CONV/DF.
.. [#f3] threading through BLAS routines only
.. [#f4] DFT gradients only implemented for SCF type DF. LRC-DFT gradients not implemented yet. DH-DFT gradients not implemented.
.. [#f5] Both EFP/EFP and QM/EFP energies are available.
.. [#f6] Adapted from `here <https://groups.google.com/forum/#!topic/google-collections-users/m8FnCcmtC88>`_.
.. [#f7] DH-DFT only available with DF-MP2.
.. [#f8] Not all combinations of reference and algorithm available. In particular, non-RHF references only available as CONV. For CCSD(T), gradients are not available with frozen core.
.. [#f9] Orbital-optimized variant available. In particular, all references available as CONV/DF.
.. [#f10] Capabilities breakdown in great detail can be found :ref:`here <table:managedmethods>`.

.. toctree::
   :hidden:

   mrcc_table_energy
   cfour_table_energy
   cfour_table_grad
