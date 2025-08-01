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
   single: Orbital-Optimized Methods, OMP2
   single: Orbital-Optimized Methods, OMP3
   single: Orbital-Optimized Methods, OMP2.5
   single: Orbital-Optimized Methods, OLCCD

.. index::
   pair: Orbital-Optimized Methods; theory
   pair: OMP2; theory
   pair: OMP3; theory
   pair: OLCCD; theory

.. _`sec:occ_oo`:

OCC: Orbital-Optimized Coupled-Cluster and |MollerPlesset| Perturbation Theories
================================================================================

.. codeauthor:: Ugur Bozkaya
.. sectionauthor:: Ugur Bozkaya

*Module:* :ref:`Keywords <apdx:occ>`, :ref:`PSI Variables <apdx:occ_psivar>`, :source:`OCC <psi4/src/psi4/occ>`

*Module:* :ref:`Keywords <apdx:dfocc>`, :ref:`PSI Variables <apdx:dfocc_psivar>`, :source:`DFOCC <psi4/src/psi4//dfocc>`

Introduction
~~~~~~~~~~~~

Orbital-optimized methods have several advantages over their non-optimized counterparts. 
Once the orbitals are optimized, the wave function will obey the Hellmann--Feynman theorem
for orbital rotation parameters. Therefore, there is no need for orbital response terms 
in the evaluation of analytic gradients. In other words, it is unnecessary to solve the 
first order coupled-perturbed CC and many-body perturbation theory (MBPT) equations. 
Further, computation of one-electron properties is easier because there are no response contributions to the particle 
density matrices (PDMs). Moreover, active space approximations can be readily incorporated into the CC methods 
[Krylov:2000:vod]_. Additionally, orbital-optimized coupled-cluster avoids spurious second-order 
poles in its response function, and its transition dipole moments are gauge invariant [Pedersen:1999:od]_.

Another advantage is that the orbital-optimized methods do not suffer from artifactual symmetry-breaking 
instabilities [Crawford:1997:instability]_, [Sherrill:1998:od]_, [Bozkaya:2011:omp2]_, and [Bozkaya:2011:omp3]_.
Furthermore, Kurlancheek and Head-Gordon [Kurlancek:2009]_ demonstrated that first order properties such as 
forces or dipole moments are discontinuous along nuclear coordinates when such a symmetry breaking occurs. 
They also observed that although the energy appears well behaved, the MP2 method can have natural occupation 
numbers greater than 2 or less than 0, hence may violate the N-representability condition. They further 
discussed that the orbital response equations generally have a singularity problem at the unrestriction point 
where spin-restricted orbitals become unstable to unrestriction. This singularity yields to extremely large or 
small eigenvalues of the one-particle density matrix (OPDM). These abnormal eigenvalues may lead to unphysical 
molecular properties such as vibrational frequencies. However, orbital-optimized MP2 (also MP3)
will solve this N-representability problem by disregarding orbital response contribution of one-particle
density matrix. 

Although the performance of coupled-cluster singles and doubles (CCSD) and orbital-optimized 
CCD (OD) is similar, the situation is different in the case of triples corrections, especially at stretched 
geometries [Bozkaya:2012:odtl]_. Bozkaya and Schaefer demonstrated that orbital-optimized coupled cluster based 
triple corrections, especially those of asymmetrics, provide significantly better potential energy curves than 
CCSD based triples corrections.  

A lot of the functionality in OCC has been enabled with Density Fitting (DF) and Cholesky 
Decomposition (CD) techniques, which can greatly speed up calculations and reduce memory
requirements for typically negligible losses in accuracy.

**NOTE**: As will be discussed later, all methods with orbital-optimization functionality have non-orbital 
optimized counterparts. Consequently, there arise two possible ways to call density-fitted MP2. In most
cases, users should prefer the DF-MP2 code described in the :ref:`DF-MP2 <sec:dfmp2>` section because it is
faster. If gradients are needed (like in a geometry optimization), then the procedures outlined hereafter
should be followed.
In general, choose the desired method, reference, and ERI type (*e.g.*,
``set reference uhf``, ``set mp2_type df``, ``opt('mp2')``) and the most
efficient module will be selected automatically, according to
:ref:`Cross-module Redundancies <table:managedmethods>`.

Thus, there arise a few categories of method, each with corresponding input keywords:

* Orbital-optimized MP and CC methods with conventional integrals (:ref:`OMP Methods <sec:occ_oo_mtds>` OCC keywords)
* Orbital-optimized MP and CC methods with DF and CD integrals (:ref:`OMP Methods <sec:occ_oo_mtds>` DFOCC keywords)
* Non-orbital-optimized MP and CC methods with conventional integrals (:ref:`MP/CC Methods <sec:occ_nonoo>` OCC keywords)
* Non-orbital-optimized MP and CC methods with DF and CD integrals (:ref:`MP/CC Methods <sec:occ_nonoo>` DFOCC keywords)

Theory
~~~~~~

What follows is a very basic description of orbital-optimized |MollerPlesset| perturbation
theory as implemented in |Psifour|.  We will follow our previous presentations ([Bozkaya:2011:omp2]_,
[Bozkaya:2011:omp3]_, and [Bozkaya:2012:odtl]_)

The orbital variations may be expressed by means of an exponential unitary operator

.. math::
   \widetilde{\hat{p}}^{\dagger} &= e^{\hat{K}} \hat{p}^{\dagger} e^{-\hat{K}}\\
   \widetilde{\hat{p}} &= e^{\hat{K}} \ \hat{p} \ e^{-\hat{K}} \\
   | \widetilde{p} \rangle &= e^{\hat{K}} \ | p \rangle

where :math:`\hat{K}` is the orbital rotation operator


.. math::
   \hat{K} &= \sum_{p,q}^{} K_{pq} \ \hat{E}_{pq} = \sum_{p>q}^{} \kappa_{pq} \ \hat{E}_{pq}^{-} \\
   \hat{E}_{pq}  &= \hat{p}^{\dagger} \hat{q} \\
   \hat{E}_{pq}^{-} &= \hat{E}_{pq} \ - \ \hat{E}_{qp} \\
   {\bf K} &= Skew({\bf \kappa}) 

The effect of the orbital rotations on the MO coefficients can be written as

.. math::
   {\bf C({\bf \kappa})} = {\bf C^{(0)}} \ e^{{\bf K}}

where :math:`{\bf C^{(0)}}` is the initial MO coefficient matrix and :math:`{\bf C({\bf \kappa})}` is the new
MO coefficient matrix as a function of :math:`{\bf \kappa}`.
Now, let us define a variational energy functional (Lagrangian) as a function of :math:`{\bf \kappa}`

* OMP2

.. math::
   \widetilde{E}({\bf \kappa}) &= \langle 0| \hat{H}^{\kappa} | 0 \rangle \\
   &+  \langle 0| \big(\hat{W}_{N}^{\kappa}\hat{T}_{2}^{(1)}\big)_{c} | 0 \rangle \\
   &+  \langle 0| \{\hat{\Lambda}_{2}^{(1)} \ \big(\hat{f}_{N}^{\kappa} \hat{T}_{2}^{(1)} 
   \ + \ \hat{W}_{N}^{\kappa} \big)_{c}\}_{c} | 0 \rangle
    
* OMP3

.. math::
   \widetilde{E}({\bf \kappa}) &= \langle 0| \hat{H}^{\kappa} | 0 \rangle \\
   &+ \langle 0| \big(\hat{W}_{N}^{\kappa}\hat{T}_{2}^{(1)}\big)_{c} | 0 \rangle 
   \ + \ \langle 0| \big(\hat{W}_{N}^{\kappa}\hat{T}_{2}^{(2)}\big)_{c} | 0 \rangle \\
   &+  \langle 0| \{\hat{\Lambda}_{2}^{(1)} \ \big(\hat{f}_{N}^{\kappa} \hat{T}_{2}^{(1)} 
   \ + \ \hat{W}_{N}^{\kappa} \big)_{c}\}_{c} | 0 \rangle \\
   &+ \langle 0| \{\hat{\Lambda}_{2}^{(1)} \ \big(\hat{f}_{N}^{\kappa} \hat{T}_{2}^{(2)} 
   \ + \ \hat{W}_{N}^{\kappa}\hat{T}_{2}^{(1)} \big)_{c}\}_{c} | 0 \rangle \\
   &+ \langle 0| \{\hat{\Lambda}_{2}^{(2)} \ \big(\hat{f}_{N}^{\kappa} \hat{T}_{2}^{(1)} 
   \ + \ \hat{W}_{N}^{\kappa} \big)_{c}\}_{c} | 0 \rangle  
    
* OLCCD

.. math::
   \widetilde{E}({\bf \kappa}) &= \langle 0| \hat{H}^{\kappa} | 0 \rangle 
   \ + \ \langle 0| \big(\hat{W}_{N}^{\kappa}\hat{T}_{2}\big)_{c} | 0 \rangle \\
   &+ \langle 0| \{\hat{\Lambda}_{2} \ \big(\hat{W}_{N}^{\kappa} \ + \ \hat{H}_{N}^{\kappa}\hat{T}_{2} \big)_{c}\}_{c}  | 0 \rangle

where subscript c means only connected diagrams are allowed, and 
:math:`\hat{H}^{\kappa}`, :math:`\hat{f}_{N}^{\kappa}`, and :math:`\hat{W}_{N}^{\kappa}` defined as

.. math::
   \hat{H}^{\kappa} &=  e^{-\hat{K}} \hat{H} e^{\hat{K}} \\
   \hat{f}_{N}^{\kappa} &=  e^{-\hat{K}} \hat{f}_{N}^{d} e^{\hat{K}} \\
   \hat{W}_{N}^{\kappa} &=  e^{-\hat{K}} \hat{W}_{N} e^{\hat{K}} 

where :math:`\hat{f}_{N}`, and :math:`\hat{W}_{N}` are the one- and two-electron components of normal-ordered Hamiltonian. Then, 
first and second derivatives of the energy with respect to the :math:`{\bf \kappa}` parameter at :math:`{\bf \kappa} = 0`

.. math::
   w_{pq} = \frac{\partial \widetilde{E}}{\partial \kappa_{pq}} 

.. math::
   A_{pq,rs} = \frac{\partial^2 \widetilde{E}}{\partial \kappa_{pq} \partial \kappa_{rs}} 

Then the energy can be expanded up to second-order as follows

.. math::
   \widetilde{E}^{(2)}({\bf \kappa}) = \widetilde{E}^{(0)} + {\bf \kappa^{\dagger} w}  + \frac{1}{2}~{\bf \kappa^{\dagger} A \kappa}

where :math:`{\bf w}` is the MO gradient vector, :math:`{\bf \kappa}` is the MO rotation vector,
and :math:`{\bf A}` is the MO Hessian matrix. Therefore, minimizing the energy with respect to :math:`{\bf \kappa}`
yields

.. math::
   {\bf \kappa} = -{\bf A^{-1}w}

This final equation corresponds to the usual Newton-Raphson step.

* OREMP

The REMP hybrid perturbation theory is a constrained mixture of the |MollerPlesset| perturbation theory and the
Retaining the Excitation degree perturbation theory([Fink:2006:RE]_, [Behnle:2019:REMP]_).
The mixing ratio is determined by the parameter :math':`A`:

.. math::
   \widehat{H}^{(0)}_\text{REMP} = (1-A)\widehat{H}^{(0)}_\text{RE} + A\widehat{H}^{(0)}_\text{MP}

Technically, the second order of RE corresponds to LCCD for RHF and UHF references. REMP2 and its orbital-optimized variant OREMP2
are thus straightforward to implement in a (O)LCCD program by appropriate scaling of residual vector contributions and density matrices.


* **OMP2** [Bozkaya:2011:omp2]_ and [Bozkaya:2013:omp2grad]_

* **OMP3** [Bozkaya:2011:omp3]_ , [Bozkaya:2013:omp3]_, and [Bozkaya:2013:omp3grad]_

* **OMP2.5** [Bozkaya:2011:omp3]_

* **OLCCD** [Bozkaya:2013:ocepa]_

* **LCCD** [Bozkaya:2013:ocepa]_

* **OREMP2** [Behnle:2021:OREMP]_, and [Behnle:2022:OREMP]_

Convergence Problems
~~~~~~~~~~~~~~~~~~~~

For problematic open-shell systems, we recommend to use the ROHF or DFT orbitals as an initial guess for orbital-optimized methods. Both ROHF and 
DFT orbitals may provide better initial guesses than UHF orbitals, hence convergence may be significantly speeded up with ROHF or DFT orbitals. 
In order to use ROHF orbitals, simply ``set reference rohf``. For DFT orbitals, ``set reference uks`` and ``set dft_functional b3lyp``. Of
course users can use any DFT functional available in |PSIfour|.

.. _`sec:occ_oo_mtds`:

Methods
~~~~~~~

The various orbital-optimized methods supported by the OCC/DFOCC
modules in |Psifour| are summarized in Table :ref:`OCC OO
Methods <table:occ_oo_calls>` and detailed in Table :ref:`OCC
OO Capabilities <table:occ_stdsuite_oo>`. Note that while two
separate libraries OCC (conventional integrals ``CONV``) and DFOCC
(density-fitted ``DF`` and Cholesky-decomposed ``CD``) together provide
the methods described on this page, they are controlled through one
|globals__qc_module| value ``OCC``. Without ``set qc_module occ``,
these methods may default to implementations in :ref:`other modules
<table:managedmethods>` based on efficiency considerations.

.. _`table:occ_oo_calls`:

.. table:: Orbital-optimized theoretical methods accessible through OCC/DFOCC

   +-------------------------+---------------------------------------------------------------+-------+
   | name                    | calls method                                                  |  OO   |
   +=========================+===============================================================+=======+
   | omp2                    | orbital-optimized second-order MP perturbation theory         |  E/G  |
   +-------------------------+---------------------------------------------------------------+-------+
   | omp2.5                  | orbital-optimized average of MP2 and MP3                      |  E/G  |
   +-------------------------+---------------------------------------------------------------+-------+
   | omp3                    | orbital-optimized third-order MP perturbation theory          |  E/G  |
   +-------------------------+---------------------------------------------------------------+-------+
   | oremp2                  | orbital-optimized second-order REMP hybrid PT                 |  E/G  |
   +-------------------------+---------------------------------------------------------------+-------+
   | olccd                   | orbital-optimized linear coupled cluster doubles              |  E/G  |
   +-------------------------+---------------------------------------------------------------+-------+

.. include:: autodoc_capabilities_occ_oo.rst


.. _`table:occ_scsoo_calls`:

.. table:: Spin-Component-Scaled Orbital-Optimized MP capabilities of OCC/DFOCC modules

    +-------------------------+--------------------------------------------------------------+----------------------+----------------------+
    | name                    | calls method                                                 |  Energy              | Gradient             |
    +=========================+==============================================================+======================+======================+
    | scs-omp3                | Spin-Component Scaled Orbital-Optimized MP3                  | RHF/UHF/ROHF/RKS/UKS | ---                  |
    +-------------------------+--------------------------------------------------------------+----------------------+----------------------+
    | sos-omp3                | Spin-Opposite Scaled Orbital-Optimized MP3                   | RHF/UHF/ROHF/RKS/UKS | ---                  |
    +-------------------------+--------------------------------------------------------------+----------------------+----------------------+
    | scs(n)-omp3             | A special version of SCS-OMP3 for nucleobase interactions    | RHF/UHF/ROHF/RKS/UKS | ---                  |
    +-------------------------+--------------------------------------------------------------+----------------------+----------------------+
    | scs-omp3-vdw            | A special version of SCS-OMP3 (from ethene dimers)           | RHF/UHF/ROHF/RKS/UKS | ---                  |
    +-------------------------+--------------------------------------------------------------+----------------------+----------------------+
    | sos-pi-omp3             | A special version of SOS-OMP3 for :math:`\pi`-systems        | RHF/UHF/ROHF/RKS/UKS | ---                  |
    +-------------------------+--------------------------------------------------------------+----------------------+----------------------+
    | scs-omp2                | Spin-Component Scaled Orbital-Optimized MP2                  | RHF/UHF/ROHF/RKS/UKS | ---                  |
    +-------------------------+--------------------------------------------------------------+----------------------+----------------------+
    | sos-omp2                | Spin-Opposite Scaled Orbital-Optimized MP2                   | RHF/UHF/ROHF/RKS/UKS | ---                  |
    +-------------------------+--------------------------------------------------------------+----------------------+----------------------+
    | scs(n)-omp2             | A special version of SCS-OMP2 for nucleobase interactions    | RHF/UHF/ROHF/RKS/UKS | ---                  |
    +-------------------------+--------------------------------------------------------------+----------------------+----------------------+
    | scs-omp2-vdw            | A special version of SCS-OMP2 (from ethene dimers)           | RHF/UHF/ROHF/RKS/UKS | ---                  |
    +-------------------------+--------------------------------------------------------------+----------------------+----------------------+
    | sos-pi-omp2             | A special version of SOS-OMP2 for :math:`\pi`-systems        | RHF/UHF/ROHF/RKS/UKS | ---                  |
    +-------------------------+--------------------------------------------------------------+----------------------+----------------------+

.. comment    | scs-ocepa               | Spin-Component Scaled Orbital-Optimized CEPA                 | RHF/UHF/ROHF/RKS/UKS | ---                  |
.. comment    | sos-ocepa               | Spin-Opposite Scaled Orbital-Optimized CEPA                  | RHF/UHF/ROHF/RKS/UKS | ---                  |
.. comment    | scs-mi-omp2             | A special version of SCS-OMP2 (from S22 database)            |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
.. comment    | scs-mi-omp3             | A special version of SCS-OMP3 (from S22 database)            |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |

.. index:: OMP2; setting keywords
.. index:: OMP3; setting keywords
.. index:: OMP2.5; setting keywords
.. index:: OLCCD; setting keywords
.. index:: OREMP; setting keywords

Basic OCC Keywords
~~~~~~~~~~~~~~~~~~

.. include:: /autodir_options_c/occ__e_convergence.rst
.. include:: /autodir_options_c/occ__r_convergence.rst
.. include:: /autodir_options_c/occ__rms_mograd_convergence.rst
.. include:: /autodir_options_c/occ__max_mograd_convergence.rst
.. include:: /autodir_options_c/occ__mo_maxiter.rst
.. include:: /autodir_options_c/occ__wfn_type.rst
.. include:: /autodir_options_c/occ__orb_opt.rst

Advanced OCC Keywords
~~~~~~~~~~~~~~~~~~~~~

.. include:: /autodir_options_c/occ__opt_method.rst
.. include:: /autodir_options_c/occ__orth_type.rst
.. include:: /autodir_options_c/occ__mp2_os_scale.rst
.. include:: /autodir_options_c/occ__mp2_ss_scale.rst
.. include:: /autodir_options_c/occ__mp2_sos_scale.rst
.. include:: /autodir_options_c/occ__mp2_sos_scale2.rst
.. include:: /autodir_options_c/occ__nat_orbs.rst
.. include:: /autodir_options_c/occ__occ_orbs_print.rst
.. include:: /autodir_options_c/occ__tpdm_abcd_type.rst
.. include:: /autodir_options_c/occ__do_diis.rst
.. include:: /autodir_options_c/occ__do_level_shift.rst

Basic DFOCC Keywords
~~~~~~~~~~~~~~~~~~~~

.. include:: /autodir_options_c/dfocc__e_convergence.rst
.. include:: /autodir_options_c/dfocc__r_convergence.rst
.. include:: /autodir_options_c/dfocc__rms_mograd_convergence.rst
.. include:: /autodir_options_c/dfocc__max_mograd_convergence.rst
.. include:: /autodir_options_c/dfocc__mo_maxiter.rst
.. include:: /autodir_options_c/dfocc__orb_opt.rst

Advanced DFOCC Keywords
~~~~~~~~~~~~~~~~~~~~~~~

.. include:: /autodir_options_c/dfocc__opt_method.rst
.. include:: /autodir_options_c/dfocc__hess_type.rst
.. include:: /autodir_options_c/dfocc__mo_diis_num_vecs.rst
.. include:: /autodir_options_c/dfocc__orth_type.rst
.. include:: /autodir_options_c/dfocc__do_diis.rst
.. include:: /autodir_options_c/dfocc__do_level_shift.rst



.. _`sec:occ_nonoo`:

Conventional (Non-OO) Coupled-Cluster and |MollerPlesset| Perturbation Theories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The various non-orbital-optimized methods supported by the OCC/DFOCC
modules in |Psifour| are summarized in Table :ref:`OCC non-OO
Methods <table:occ_nonoo_calls>` and detailed in Table :ref:`OCC
non-OO Capabilities <table:occ_stdsuite_nonoo>`. Note that while two
separate libraries OCC (conventional integrals ``CONV``) and DFOCC
(density-fitted ``DF`` and Cholesky-decomposed ``CD``) together provide
the methods described on this page, they are controlled through one
|globals__qc_module| value ``OCC``. Without ``set qc_module occ``,
these methods may default to implementations in :ref:`other modules
<table:managedmethods>` based on efficiency considerations.

Starting in v1.4, MP2.5 and MP3 default to the density-fit algorithm. Set |globals__mp_type| to ``CONV`` to get previous behavior.

Publications resulting from the use of the non-OO CC codes should cite the following publications:

* **MP2** [Bozkaya:2011:omp2]_, [Bozkaya:2013:omp2grad]_, and [Bozkaya:2014:dfomp2grad]_

* **MP3** [Bozkaya:2011:omp3]_ , [Bozkaya:2013:omp3]_, [Bozkaya:2013:omp3grad]_, [Bozkaya:2016:dfomp3]_, and [Bozkaya:2018:dfomp3grad]_

* **MP2.5** [Bozkaya:2011:omp3]_, [Bozkaya:2016:dfomp3]_, and [Bozkaya:2018:dfomp3grad]_

* **REMP** [Behnle:2019:REMP]_, [Behnle:2022:OREMP]_

* **LCCD** [Bozkaya:2013:ocepa]_ and [Bozkaya:2016:dfolccd]_

* **CCD** [Bozkaya:2020:dfoccd]_

* **CCSD** [Bozkaya:2016:dfccsdgrad]_

* **CCSD(T)** [Bozkaya:2017:dfccsdtgrad]_

* **CCSD(AT)** [Bozkaya:2016:dfccsdat]_


.. _`table:occ_nonoo_calls`:

.. table:: Non-OO theoretical methods accessible through OCC/DFOCC

   +-------------------------+---------------------------------------------------------------+-------+-------+
   | name                    | calls method                                                  | plain |  FNO  |
   +=========================+===============================================================+=======+=======+
   | mp2                     | second-order MP perturbation theory                           |  E/G  |  n/a  |
   +-------------------------+---------------------------------------------------------------+-------+-------+
   | mp2.5                   | average of MP2 and MP3                                        |  E/G  |       |
   +-------------------------+---------------------------------------------------------------+-------+-------+
   | mp3                     | third-order MP perturbation theory                            |  E/G  |       |
   +-------------------------+---------------------------------------------------------------+-------+-------+
   | remp2                   | second-order retaining-the-excitation-degree MP hybrid PT     |  E    |       |
   +-------------------------+---------------------------------------------------------------+-------+-------+
   | lccd                    | linear coupled cluster doubles                                |  E/G  |       |
   +-------------------------+---------------------------------------------------------------+-------+-------+
   | ccd                     | coupled cluster doubles                                       |  E/G  |       |
   +-------------------------+---------------------------------------------------------------+-------+-------+
   | ccsd                    | coupled cluster singles and doubles                           |  E/G  |       |
   +-------------------------+---------------------------------------------------------------+-------+-------+
   | ccsd(t)                 | coupled cluster singles and doubles with perturbative triples |  E/G  |       |
   +-------------------------+---------------------------------------------------------------+-------+-------+
   | a-ccsd(t)               | CCSD with asymmetric perturbative triples                     |  E    |       |
   +-------------------------+---------------------------------------------------------------+-------+-------+

.. include:: autodoc_capabilities_occ_nonoo.rst

