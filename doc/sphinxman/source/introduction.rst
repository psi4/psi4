
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
(e.g., Hartree |--| Fock, MP2, coupled-cluster) and general procedures
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

Overall PSI4 Package
^^^^^^^^^^^^^^^^^^^^

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

* "Density Cumulant Functional Theory: First Implementation and
  Benchmark Results for the DCFT-06 Model," A. C. Simmonett,
  J. J. Wilke, H. F. Schaefer, and W. Kutzelnigg, *J. Chem. Phys.*
  **133**, 174122 (2010).
  (doi: `10.1063/1.3503657 <http://dx.doi.org/10.1063/1.3503657>`_).

* "Analytic gradients for density cumulant functional theory: The 
  DCFT-06 model," A. Yu. Sokolov, J. J. Wilke, A. C. Simmonett, 
  and H. F. Schaefer, *J. Chem. Phys.* **137**, 054105 (2012).
  (doi: `10.1063/1.4739423 <http://dx.doi.org/10.1063/1.4739423>`_).

* "Density cumulant functional theory: The DC-12 method, an improved 
  description of the one-particle density matrix," A. Yu. Sokolov, 
  A. C. Simmonett, and H. F. Schaefer, *J. Chem. Phys.*  **138**, 024107 
  (2013).
  (doi: `10.1063/1.4773580 <http://dx.doi.org/10.1063/1.4773580>`_).

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
Symmetry Adapted Perturbation Theory (SAPT).  A good review article for this 
method is as follows:

* "Perturbation Theory Approach to Intermolecular Potential Energy
  Surfaces of van der Waals Complexes," B. Jeziorski, R. Moszynski,
  and K. Szalewicz, *Chem. Rev.* **94**, 1887-1930 (1994).   
  (doi: `10.1021/cr00031a008 <http://dx.doi.org/10.1021/cr00031a008>`_).

|PSIfour| benefits enormously from the introduction of density fitting (DF)
into SAPT.  There are several SAPT truncations available in PSI.  For 
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


Orbital-Optimized Post-Hartree-Fock Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Orbital-optimized second-order perturbation theory (OMP2)

* "Quadratically convergent algorithm for orbital optimization in the 
  orbital-optimized coupled-cluster doubles method and in orbital-optimized 
  second-order M\ |o_slash|\ ller--Plesset perturbation theory," 
  U. Bozkaya, J. M. Turney, Y. Yamaguchi, H. F. Schaefer, and C. D. Sherrill,
  *J. Chem. Phys.* **135**, 104103 (2011).
  (doi: `10.1063/1.3631129 <http://dx.doi.org/10.1063/1.3631129>`_).

* "Analytic energy gradients for the orbital-optimized second-order 
  M\ |o_slash|\ ller--Plesset perturbation theory," U. Bozkaya and 
  C. D. Sherrill, *J. Chem. Phys.* **138**, 184103 (2013).
  (doi: `10.1063/1.4803662 <http://dx.doi.org/10.1063/1.4803662>`_).

* "Orbital-Optimized Second-Order Perturbation Theory with Density-Fitting
  and Cholesky Decomposition Approximations: An Efficient Implementation,"
  U. Bozkaya,   *J. Chem. Theory Comput.* **10**, 2371 (2014).
  (doi: `10.1021/ct500231c <http://dx.doi.org/10.1021/ct500231c>`_).

Orbital-optimized third-order perturbation theory (OMP3)

* "Orbital-Optimized Third-Order M\ |o_slash|\ ller--Plesset Perturbation 
  Theory and Its Spin-Component and Spin-Opposite Scaled Variants: Application 
  to Symmetry Breaking Problems," U. Bozkaya,
  *J. Chem. Phys.* **135**, 224103 (2011).
  (doi: `10.1063/1.3665134 <http://dx.doi.org/10.1063/1.3665134>`_).

* "Assessment of Orbital-Optimized Third-Order M\ |o_slash|\ ller--Plesset 
  Perturbation Theory and Its Spin-Component and Spin-Opposite Scaled Variants 
  for Thermochemistry and Kinetics," E. Soydas and U. Bozkaya,  
  *J. Chem. Theory Comput.* **9**, 1452 (2013).
  (doi: `10.1021/ct301078q <http://dx.doi.org/10.1021/ct301078q>`_).

* "Analytic energy gradients for the orbital-optimized third-order M\ |o_slash|\ ller--Plesset 
  Perturbation Theory," U. Bozkaya,  
  *J. Chem. Phys.* **139**, 104116 (2013).
  (doi: `10.1063/1.4820877 <http://dx.doi.org/10.1063/1.4820877>`_).

Orbital-optimized coupled electron pair approximation (OCEPA)

* "Orbital-optimized coupled-electron pair theory and its analytic gradients: 
  Accurate equilibrium geometries, harmonic vibrational frequencies, and hydrogen transfer 
  reactions," U. Bozkaya and C. D. Sherrill,
  *J. Chem. Phys.* **139**, 054104 (2013).
  (doi: `10.1063/1.4816628 <http://dx.doi.org/10.1063/1.4816628>`_).

Orbital-optimized MP2.5 (OMP2.5)

* "Orbital-Optimized Third-Order M\ |o_slash|\ ller--Plesset Perturbation 
  Theory and Its Spin-Component and Spin-Opposite Scaled Variants: Application 
  to Symmetry Breaking Problems," U. Bozkaya,
  *J. Chem. Phys.* **135**, 224103 (2011).
  (doi: `10.1063/1.3665134 <http://dx.doi.org/10.1063/1.3665134>`_).

* U. Bozkaya and C. D. Sherrill, (unpublished).

Extended Koopmans' Theorem

* "The extended Koopmans' theorem for orbital-optimized methods: Accurate computation of ionization potentials," 
   U. Bozkaya,  *J. Chem. Phys.* **139**, 154105 (2013).
  (doi: `10.1063/1.4825041 <http://dx.doi.org/10.1063/1.4825041>`_).

* "Accurate Electron Affinities from the Extended Koopmans' Theorem Based on Orbital-Optimized Methods,"
  U. Bozkaya,   *J. Chem. Theory Comput.* **10**, 2041 (2014).
  (doi: `10.1021/ct500186j <http://dx.doi.org/10.1021/ct500186j>`_).


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


.. index:: architectures
.. index:: compilers

Supported Architectures
=======================

The majority of |PSIfour| was developed on Mac and Linux machines.  In
principle, it should work on any Unix system; however, we have not tested
extensively on systems other than Mac and Linux.  There is not a Windows
version of |PSIfour|.

|PSIfour| has been successfully compiled using Intel, GCC, and Clang
compilers.  For the Intel compilers, use versions 11 or 12.1 (we have had
trouble with version 12.0). See Sec. :ref:`Compiling and Installing
<sec:installFile>` for details.


Capabilities
============

|PSIfour| can perform *ab initio* computations employing basis
sets of contrated Gaussian-type functions of virtually arbitrary
orbital quantum number.  Many parts of |PSIfour| can recognize and
exploit the largest Abelian subgroup of the molecular point group.
Table :ref:`Methods <table:methods>` displays the range of theoretical methods
available in |PSIfour|.
For more details, see Tables :ref:`Energy <table:energy_gen>`, 
:ref:`Energy (DFT) <table:energy_dft>`, :ref:`Energy (MRCC) <table:energy_mrcc>`,
:ref:`Energy (CFOUR) <table:energy_cfour>`, :ref:`Gradient <table:grad_gen>`, 
:ref:`Gradient (CFOUR) <table:grad_cfour>`, and :ref:`Frequency <table:freq_gen>`.

.. _`table:methods`:

.. table:: Summary of theoretical methods available in |PSIfour|

    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | Method                  | Energy    | Gradient  | Reference            | Parallelism                 |
    +=========================+===========+===========+======================+=============================+
    | SCF (HF and DFT)        | Y         | Y [#f4]_  | RHF/ROHF/UHF/RKS/UKS | threaded                    |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | DF-SCF (HF and DFT)     | Y         | Y [#f4]_  | RHF/ROHF/UHF/RKS/UKS | threaded                    |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | DCFT                    | Y         | Y         | UHF                  | partially threaded          |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | MP2                     | Y         | Y         | RHF/ROHF/UHF         | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | DF-MP2                  | Y         | Y [#f2]_  | RHF/ROHF/UHF         | threaded                    |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | MP3                     | Y         | Y         | RHF/UHF              | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | MP2.5                   | Y         | Y         | RHF/UHF              | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | MP4                     | Y         | ---       | RHF                  | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | MP(n)                   | Y         | ---       | RHF/ROHF             | partially threaded          |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | ZAPT(n)                 | Y         | ---       | RHF/ROHF             | partially threaded          |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | OMP2                    | Y         | Y         | RHF/ROHF/UHF/RKS/UKS | partially threaded          |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | OMP3                    | Y         | Y         | RHF/ROHF/UHF/RKS/UKS | partially threaded          |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | OMP2.5                  | Y         | Y         | RHF/ROHF/UHF/RKS/UKS | partially threaded          |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | OCEPA                   | Y         | Y         | RHF/ROHF/UHF/RKS/UKS | partially threaded          |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | CEPA(0)                 | Y         | Y         | RHF/UHF              | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | DF-OMP2                 | Y         | Y         | RHF/ROHF/UHF/RKS/UKS | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | CEPA(n), n=0,1,3        | Y         | ---       | RHF                  | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | ACPF/AQCC               | Y         | ---       | RHF                  | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | QCISD                   | Y         | ---       | RHF                  | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | QCISD(T)                | Y         | ---       | RHF                  | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | CC2                     | Y         | ---       | RHF/ROHF/UHF         | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | CCSD                    | Y         | Y         | RHF/ROHF/UHF         | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | DF-CCSD                 | Y         | ---       | RHF                  | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | CCSD(T)                 | Y         | Y [#f1]_  | RHF/ROHF/UHF         | threaded (pthreads)         |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | DF-CCSD(T)              | Y         | ---       | RHF                  | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | CC3                     | Y         | ---       | RHF/ROHF/UHF         | threaded (pthreads)         |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | Mk-MRPT2                | Y         | ---       | RHF/ROHF/TCSCF       | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | Mk-MRCCSD               | Y         | ---       | RHF/ROHF/TCSCF       | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | Mk-MRCCSD(T)            | Y         | ---       | RHF/ROHF/TCSCF       | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | CI(n)                   | Y         | ---       | RHF/ROHF             | partially threaded          |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | RAS-CI                  | Y         | ---       | RHF/ROHF             | partially threaded          |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | SAPT                    | Y         | ---       | RHF                  | threaded                    |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | CIS/RPA/TDHF            | Y         | ---       |                      |                             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | ADC(2)                  | Y         | ---       | RHF                  | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+
    | EOM-CCSD                | Y         | Y         | RHF/ROHF/UHF         | threaded [#f3]_             |
    +-------------------------+-----------+-----------+----------------------+-----------------------------+

..    | %HF DBOC                | Y         | N         |
    +-------------------------+-----------+-----------+
    | %TCSCF                  | Y         | Y         |
    +-------------------------+-----------+-----------+
    | %CASSCF                 | Y         | Y         |
    +-------------------------+-----------+-----------+
    | %RASSCF                 | Y         | Y         |
    +-------------------------+-----------+-----------+
    | %RHF MP2-R12            | Y         | N         |
    +-------------------------+-----------+-----------+
    | %RAS-CI DBOC            | Y         | N         |
    +-------------------------+-----------+-----------+

Geometry optimization can be performed using either analytic gradients
or energy points.  Likewise, vibrational frequencies can be 
computed by analytic second derivatives, by finite
differences of analytic gradients, or by finite differences of energies.
|PSIfour| can also compute an extensive list of one-electron properties.

.. index::
   single: contact
   single: bugs

Technical Support
=================

The |PSIfour| package is
distributed for free and without any guarantee of reliability,
accuracy, or suitability for any particular purpose.  No obligation
to provide technical support is expressed or implied.  As time
allows, the developers will attempt to answer inquiries directed to
`crawdad@vt.edu <mailto:crawdad@vt.edu>`_
or `sherrill@gatech.edu <mailto:sherrill@gatech.edu>`_.
For bug reports, specific and detailed information, with example
inputs, would be appreciated.  Questions or comments regarding
this user's manual may be sent to 
`sherrill@gatech.edu <mailto:sherrill@gatech.edu>`_.

Alternatively, bug reports and comments can be submitted to the `Issue
tracker on GitHub <https://github.com/psi4/psi4.0b4/issues/new>`_ . This site
is viewable by all, but reporting bugs requires signing up for a `free
GitHub account <https://github.com/signup/free>`_.


.. rubric:: Footnotes

.. [#f1] UHF-CCSD(T) gradients only, as of |version|
.. [#f2] RHF and UHF reference are available, however the latter one should be requsted from DFOCC module.  DF-MP2 is recommended as a faster alternative.
.. [#f3] threading through BLAS routines only
.. [#f4] DFT gradients only implemented for SCF type DF. LRC-DFT gradients not implemented yet. 

.. toctree::
   :hidden:

   mrcc_table_energy
   cfour_table_energy
   cfour_table_grad
