
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
structure of |PSIfour| input files, and how Python can be mixed with
quantum chemistry directives in |PSIfour|. Section :ref:`Psithon Functions <sec:psithonFunc>`
provides more detail on the Python functions provided by |PSIfour|
and discusses some of the higher-level functions such as counterpoise
correction, complete-basis-set extrapolation, and running computations
on an entire database of molecules at a time.  Later sections deal with
the different types of computations which can be done using |PSIfour|
(e.g., Hartree |--| Fock, MP2, coupled-cluster) and general procedures
such as geometry optimization and vibrational frequency analysis.
The Appendix includes a complete description of all possible input
keywords for each module, as well as tables of available basis sets and
a listing of the sample input files available under :source:`samples`.

The user is urged to examine this directory of sample inputs, as
most common types of computations are represented there.
For the latest |PSIfour| documentation, check 
`www.psicode.org <http://www.psicode.org>`_.

Citing |PSIfour|
================

The following citation should be used in any publication utilizing the
|PSIfour| program package:

.. rubric:: `PSI4 <http://dx.doi.org/10.1002/wcms.93>`_

* "Psi4: an open-source *ab initio* electronic structure program,"
  J. M. Turney, A. C. Simmonett, R. M. Parrish, E. G. Hohenstein, F.
  Evangelista, J. T. Fermann, B. J. Mintz, L. A. Burns, J. J. Wilke, M. L.
  Abrams, N. J.  Russ, M. L. Leininger, C. L. Janssen, E. T. Seidl, W. D.
  Allen, H. F.  Schaefer, R. A. King, E. F. Valeev, C. D. Sherrill, and T.
  D. Crawford, *WIREs Comput. Mol. Sci.*, (2011) (doi: 10.1002/wcms.93)


Depending on the particular modules used, the user may also wish to
cite some of the following references for theoretical, algorithmic,
or implementation contributions specific to |PSIfour| (in addition to
appropriate references for the underlying theory):

.. rubric:: DCFT 

* "Density Cumulant Functional Theory: First Implementation and
  Benchmark Results for the DCFT-06 Model," A. C. Simmonett,
  J. J. Wilke, H. F. Schaefer, and W. Kutzelnigg, *J. Chem. Phys.*
  **133**, 174122 (2010).

.. rubric:: ADC(2)

* "Intermediate state representation approach to physical properties of electronically excited molecules,"
  J. Schirmer, and A. B. Trofimov, *J. Chem. Phys.* **120**,
  11449-11464 (2004).

.. rubric:: PR-CIS(D) and PR-ADC(2)

* "Excited State Calculation for Free-Base and Metalloporphyrins with
  the Partially Renormalized Polarization Propagator Approach,"
  M. Saitow and Y. Mochizuki, *Chem. Phys. Lett.*, in press.

.. rubric:: CI

* "The Configuration Interaction Method: Advances in Highly 
  Correlated Approaches," C. D. Sherrill and H. F. Schaefer, in
  *Adv. Quantum Chem.*, vol. 34, P.-O. L\ |o_dots|\ wdin, Ed.
  (Academic Press, New York, 1999), pp. 143-269.

.. rubric:: CC

* "An Introduction to Coupled Cluster Theory for Computational
  Chemists," T. D. Crawford and H. F. Schaefer, *Rev. Comp. Chem.* 
  **14**, 33-136 (2000).

.. rubric:: Mk-MRCCSD

* "Coupling Term Derivation and General Implementation of
  State-Specific Multireference Coupled-Cluster Theories,"
  F. A. Evangelista, W. D. Allen, and H. F. Schaefer, 
  *J. Chem. Phys.* **127**, 024102 (2007).

.. rubric:: Mk-MRCCSDT(-n)

* "Triple Excitations in State-Specific Multireference Coupled
  Cluster Theory: Application of Mk-MRCCSDT and Mk-MRCCSDT-n Methods to
  Model Systems," F. A. Evangelista, A. C. Simmonett, W. D. Allen,
  H. F. Schaefer, and J. Gauss, *J. Chem. Phys.* **128**, 124104
  (2008).

.. rubric:: SAPT (General)

All capabilities of the SAPT module are based on Symmetry Adapted
Perturbation Theory.  A good review article for this method is as
follows:

* "Perturbation Theory Approach to Intermolecular Potential Energy
  Surfaces of van der Waals Complexes," B. Jeziorski, R. Moszynski,
  and K. Szalewicz, *Chem. Rev.* **94**, 1887-1930 (1994).   

The particular implementation and algorithms for various orders of SAPT
available in |PSIfour| are provided below.

.. rubric:: SAPT0 

* "Large-scale Symmetry-adapted Perturbation Theory Computations via
  Density Fitting and Laplace Transformation Techniques: Investigating the
  Fundamental Forces of DNA-Intercalator Interactions," E. G. Hohenstein,
  R. M. Parrish, C. D. Sherrill, J. M. Turney, and H. F. Schaefer, *J.
  Chem. Phys.* **135**, 174017 (2011).

* "Density Fitting and Cholesky Decomposition Approximations
  in Symmetry-Adapted Perturbation Theory: Implementation and Application
  to Probe the Nature of :math:`\pi - \pi` Interactions in Linear Acenes,"
  E. G. Hohenstein and C. D. Sherrill, *J. Chem. Phys.* **132**,
  184111 (2010).

.. rubric:: SAPT2, SAPT2+, SAPT2+(3), SAPT2+3

* "Density Fitting of Intramonomer Correlation Effects in
  Symmetry-Adapted Perturbation Theory,"
  E. G. Hohenstein and C. D. Sherrill, *J. Chem. Phys.* **133**,
  014101 (2010).

* "Wavefunction Methods for Noncovalent Interactions," E. G.
  Hohenstein and C. D. Sherrill, *WIREs: Comput. Mol. Sci.*, in press.

.. rubric:: Using Natural Orbitals in SAPT

* "Efficient Evaluation of Triple Excitations in Symmetry-Adapted
  Perturbation Theory via MP2 Natural Orbitals," E. G. Hohenstein
  and C. D. Sherrill, *J. Chem. Phys.* **133**, 104107 (2010).

.. _`sec:installation`:

Obtaining and Installing |PSIfour|
==================================

The latest version of the |PSIfour| program package may be obtained at
`www.psicode.org <http://www.psicode.org>`_.  The
source code is available as a gzipped tar archive (named, for example,
``psi4.X.tar.gz``, and binaries may be available for certain architectures.
For detailed installation and testing instructions, please refer to the
installation instructions at the |PSIfour| website above or to the file
:source:`INSTALL` distributed with the package. Additional compilation
hints may be found at `Psi Compiling <http://sirius.chem.vt.edu/trac/wiki/CompilingPsi>`_.

To avoid dependency problems, a script is provided at
:source:`lib/scripts/psi4depend-v2.sh` that will download, configure, and
install all dependencies that the Psi4 developers currently recommend.
This script simply needs one to edit the top lines to tell it where to
install all the files to and for the user to put this into their path. The
script will print out the lines that you need to add to your
:envvar:`path` and :envvar:`LD_LIBRARY_PATH`. This script will build the
following: autoconf-2.68, automake-1.11, gcc-4.1.2, mpich2-1.2.1,
Python-2.6.6, Boost-1.48.0 (as well as some dependencies). After running
the script, proceed with building |PSIfour| as described in
:source:`INSTALL`.

.. index:: architectures
.. index:: compilers

Supported Architectures
=======================

The majority of |PSIfour| was developed on Mac and Linux machines.  In
principle, it should work on any Unix system; however, we have not tested
extensively on systems other than Mac and Linux.  There is not a Windows
version of |PSIfour|.

|PSIfour| has been successfully compiled using Intel, GCC, and Clang
compilers.  For the Intel compilers, use versions 11 or
12.1 (we have had trouble with version 12.0).  


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
:ref:`Gradient <table:grad_gen>`, and :ref:`Frequency <table:freq_gen>`.

.. _`table:methods`:

.. comment table:: Summary of theoretical methods available in |PSIfour|

.. comment    +-------------------------+-----------+-----------+
.. comment    | Method                  | Energy    | Gradient  |
.. comment    +=========================+===========+===========+
.. comment    | RHF/ROHF/UHF SCF        | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF/ROHF/UHF DF-SCF     | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | CIS/RPA/TDHF            | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | UHF DCFT                | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF SAPT                | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF MP2                 | Y         | Y         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | UHF/ROHF MP2            | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF DF-MP2              | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF ADC(2)              | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF/ROHF CI(n)          | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF/ROHF RAS-CI         | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF/ROHF MP(n)          | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF/ROHF ZAPT(n)        | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF/UHF/ROHF CC2        | Y         | N         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF/UHF/ROHF CCSD       | Y         | Y         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF/UHF/ROHF CCSD(T)    | Y         | Y [#f1]_  |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF/UHF/ROHF EOM-CCSD   | Y         | Y         |
.. comment    +-------------------------+-----------+-----------+
.. comment    | RHF/UHF/ROHF CC3        | Y         | N         |
.. comment    +-------------------------+-----------+-----------+

.. table:: Summary of theoretical methods available in |PSIfour|

    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | Method                  | Energy    | Gradient  | Reference       | Parallelism                 | Verified By   |
    +=========================+===========+===========+=================+=============================+===============+
    | SCF                     | Y         | ---       | RHF/ROHF/UHF    | threaded                    | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | DF-SCF                  | Y         | ---       | RHF/ROHF/UHF    | threaded                    | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | CIS/RPA/TDHF            | Y         | ---       |                 |                             | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | DCFT                    | Y         | ---       | UHF             |                             | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | SAPT                    | Y         | ---       | RHF             | threaded                    | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | MP2                     | Y         | Y [#f2]_  | RHF/ROHF/UHF    |                             | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | DF-MP2                  | Y         | ---       | RHF             | threaded                    | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | ADC(2)                  | Y         | ---       | RHF/ROHF        |                             | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | CI(n)                   | Y         | ---       | RHF/ROHF        | threaded (pthreads)         |               |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | RAS-CI                  | Y         | ---       | RHF/ROHF        |                             | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | MP(n)                   | Y         | ---       | RHF/ROHF        |                             | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | ZAPT(n)                 | Y         | ---       | RHF/ROHF        |                             | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | CC2                     | Y         | ---       | RHF/ROHF/UHF    |                             | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | CCSD                    | Y         | Y         | RHF/ROHF/UHF    | threaded [#f3]_             |               |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | CCSD(T)                 | Y         | Y [#f1]_  | RHF/ROHF/UHF    | threaded (pthreads) [#f2]_  |               |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | CC3                     | Y         | ---       | RHF/ROHF/UHF    |                             | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+
    | EOM-CCSD                | Y         | Y         | RHF/ROHF/UHF    |                             | UNVERIFIED    |
    +-------------------------+-----------+-----------+-----------------+-----------------------------+---------------+

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

Geometry optimization (currently restricted to true minima on the potential
energy surface) can be performed using either analytic gradients
or energy points.  Likewise, vibrational frequencies can be 
computed by analytic second derivatives, by finite
differences of analytic gradients, or by finite differences of energies.
|PSIfour| can also compute an extensive list of one-electron properties.

.. index:: contact

Technical Support
=================

The |PSIfour| package is
distributed for free and without any guarantee of reliability,
accuracy, or suitability for any particular purpose.  No obligation
to provide technical support is expressed or implied.  As time
allows, the developers will attempt to answer inquiries directed to
`crawdad@vt.edu <mailto:crawdad@vt.edu>`_.
For bug reports, specific and detailed information, with example
inputs, would be appreciated.  Questions or comments regarding
this user's manual may be sent to 
`sherrill@gatech.edu <mailto:sherrill@gatech.edu>`_.


.. rubric:: Footnotes

.. [#f1] UHF-CCSD(T) gradients only, as of |version|
.. [#f2] RHF reference only
.. [#f3] threading through BLAS routines only


