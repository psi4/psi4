
.. include:: autodoc_abbr_options_c.rst

.. index:: 
   single: Density-Fitted Orbital-Optimized Methods, DF-OMP2

.. index::
   pair: Orbital-Optimized Methods; theory
   pair: DF-OMP2; theory

.. _`sec:dfocc`:

DF-OCC: Density-Fitted Orbital-Optimized Coupled-Cluster and M\ |o_slash|\ ller--Plesset Perturbation Theories 
==============================================================================================================

.. codeauthor:: Ugur Bozkaya
.. sectionauthor:: Ugur Bozkaya

*Module:* :ref:`Keywords <apdx:dfocc>`, :ref:`PSI Variables <apdx:dfocc_psivar>`, :source:`DFOCC <src/bin/dfocc>`

Introduction
~~~~~~~~~~~~

Orbital-optimized methods have several advantages over non-optimized counterparts. 
Once the orbitals are optimized, the wave function will obey the Hellmann-Feynman theorem 
for orbital rotation parameters. Therefore, there is no need for orbital response terms 
in the evaluation of analytic gradients. In other words, it is unnecessary to solve the 
first order coupled-perturbed CC and many-body perturbation theory (MBPT) equations. 
Further, computation of one-electron properties is easier because there are no response contributions to the particle 
density matrices (PDMs). Moreover, active space approximations can be readily incorporated into the CC methods 
[Krylov:2000:vod]_. Additionally, orbital-optimized coupled-cluster avoids spurious second-order 
poles in its response function, and its transition dipole moments are gauge invarianti [Pedersen:1999:od]_.

Another advantage is that the orbital-optimized methods does not suffer from artifactual symmetry-breaking 
instabilities [Crawford:1997:instability]_, [Sherrill:1998:od]_, [Bozkaya:2011:omp2]_, and [Bozkaya:2011:omp3]_.
Further, Kurlancheek and Head-Gordon [Kurlancek:2009]_ demonstrated that first order properties such as 
forces or dipole moments are discontinuous along nuclear coordinates when such a symmetry breaking occurs. 
They also observed that although the energy appears well behaved, the MP2 method can have natural occupation 
numbers greater than 2 or less than 0, hence may violate the N-representability condition. They further 
discussed that the orbital response equations generally have a singularity problem at the unrestriction point 
where spin-restricted orbitals become unstable to unrestriction. This singularity yields to extremely large or 
small eigenvalues of the one-particle density matrix (OPDM). These abnormal eigenvalues may lead to unphysical 
molecular properties such as vibrational frequencies. However, orbital optimized MP2 (hence Orbital optimized MP3) 
will solve this N-representability problem by disregarding orbital response contribution of one-partical 
density matrix. 

Although the performance of coupled-cluster singles and doubles (CCSD) and orbital-optimized 
CCD (OD) is similar, the situation is different in the case of triples corrections, especially at stretched 
geometries [Bozkaya:2012:odtl]_. Bozkaya and Schaefer demonstrated that orbital-optimized coupled cluster based 
triple corrections, especially those of asymmetrics, provide significantly better potential energy curves than 
CCSD based triples corrections.  

Theory 
~~~~~~

What follows is a very basic description of orbital-optimized M\ |o_slash|\ ller--Plesset perturbation 
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

Publications resulting from the use of the OMP2 code should cite the following publications: 

[Bozkaya:2011:omp2]_ and [Bozkaya:2013:omp2grad]_.

Convergence Problems
~~~~~~~~~~~~~~~~~~~~

For problematic open-shell systems, we recommend to use the ROHF or DFT orbitals as an initial guess for orbital-optimized methods. Both ROHF and 
DFT orbitals may provide better initial guesses than UHF orbitals, hence convergence may be significantly speeded up with ROHF or DFT orbitals. 
In order to use ROHF orbitals we can simply use "reference rohf" option. For DFT orbitals one should use "reference uks" and "dft_functional b3lyp" options. Of 
course users can use any DFT functional available in Psi4. 


Methods
~~~~~~~

Density-fitted conventional and orbital-optimized MP2 methods currently supported in |Psifour| are outlined in Table :ref:`DF-OMP2 Methods <table:dfomp2_calls>`.

    .. _`table:dfomp2_calls`:

    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | Name                    | Calls Method                                                 |  Energy | Gradient | Reference              |
    +=========================+==============================================================+=========+==========+========================+
    | ri-mp2                  | Density-Fitted MP2                                           |    Y    |     Y    | RHF/ROHF/UHF           |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | cd-mp2                  | Cholesky-Decomposed MP2                                      |    Y    |     N    | RHF/ROHF/UHF           |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | df-omp2                 | Density-Fitted Orbital-Optimized MP2                         |    Y    |     Y    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | cd-omp2                 | Cholesky-Decomposed Orbital-Optimized MP2                    |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+

.. index:: DF-OMP2; setting keywords

Basic Keywords
~~~~~~~~~~~~~~

.. include:: /autodir_options_c/dfocc__e_convergence.rst
.. include:: /autodir_options_c/dfocc__r_convergence.rst
.. include:: /autodir_options_c/dfocc__rms_mograd_convergence.rst
.. include:: /autodir_options_c/dfocc__max_mograd_convergence.rst
.. include:: /autodir_options_c/dfocc__mo_maxiter.rst
.. include:: /autodir_options_c/dfocc__orb_opt.rst

Advanced Keywords
~~~~~~~~~~~~~~~~~
.. include:: /autodir_options_c/dfocc__opt_method.rst
.. include:: /autodir_options_c/dfocc__hess_type.rst
.. include:: /autodir_options_c/dfocc__mo_diis_num_vecs.rst
.. include:: /autodir_options_c/dfocc__orth_type.rst
.. include:: /autodir_options_c/dfocc__do_diis.rst
.. include:: /autodir_options_c/dfocc__do_level_shift.rst


.. _`sec:dfconvocc`:

DF-OCC: Conventional M\ |o_slash|\ ller--Plesset Perturbation Theories 
======================================================================

*Module:* :ref:`Keywords <apdx:dfocc>`, :ref:`PSI Variables <apdx:dfocc_psivar>`, :source:`DFOCC <src/bin/dfocc>`

|PSIfour| also has a density-fitted MP2 algorithm for RHF, UHF, and ROHF
energies in the DFOCC module. 

    .. _`table:dfnonoo`:

    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | Name                    | Calls Method                                                 |  Energy | Gradient | Reference              |
    +=========================+==============================================================+=========+==========+========================+
    | ri-mp2                  | DF-MP2                                                       |    Y    |     Y    | RHF/ROHF/UHF           |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | cd-mp2                  | MP2                                                          |    Y    |     N    | RHF/ROHF/UHF           |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+

