
.. include:: autodoc_abbr_options_c.rst

.. index:: 
   single: Orbital-Optimized Methods, OMP2
   single: Orbital-Optimized Methods, OMP3
   single: Orbital-Optimized Methods, OMP2.5
   single: Orbital-Optimized Methods, OCEPA

.. index::
   pair: Orbital-Optimized Methods; theory
   pair: OMP2; theory
   pair: OMP3; theory
   pair: OCEPA; theory

.. _`sec:occ`:

OCC: Orbital-Optimized Coupled-Cluster and M\ |o_slash|\ ller--Plesset Perturbation Theories 
============================================================================================

.. codeauthor:: Ugur Bozkaya
.. sectionauthor:: Ugur Bozkaya

*Module:* :ref:`Keywords <apdx:occ>`, :ref:`PSI Variables <apdx:occ_psivar>`, :source:`OCC <src/bin/occ>`

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

Another advantage is that the orbital-optimized methods does not suffer from the artifactual symmetry-breaking 
instabilities [Sherrill:1998:od]_, [Bozkaya:2011:omp2]_, and [Bozkaya:2011:omp3]_.
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
    
* OCEPA

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

Publications resulting from the use of the OMP2 code should cite the following publication(s): 

[Bozkaya:2011:omp2]_.

Publications resulting from the use of the OMP3 code should cite the following publications: 

[Bozkaya:2011:omp3]_ and [Bozkaya:2013:omp3]_.

Publications resulting from the use of the OMP2.5 code should cite the following publications: 

[Bozkaya:2011:omp3]_ and [Bozkaya:2013:omp3]_.

Publications resulting from the use of the OCEPA code should cite the following publication(s): 

[Bozkaya:2011:omp2]_.

Publications resulting from the use of the MP2 code should cite the following publication(s): 

[Bozkaya:2011:omp2]_.

Publications resulting from the use of the MP3 code should cite the following publications: 

[Bozkaya:2011:omp3]_ and [Bozkaya:2013:omp3]_.

Publications resulting from the use of the MP2.5 code should cite the following publications: 

[Bozkaya:2011:omp3]_ and [Bozkaya:2013:omp3]_.

Publications resulting from the use of the CEPA0 code should cite the following publication(s): 

[Bozkaya:2011:omp2]_.



Methods
~~~~~~~

The orbital-optimized MP2 methods currently supported in |Psifour| are outlined in Table :ref:`OMP2 Methods <table:omp2_calls>`.

    .. _`table:omp2_calls`:

    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | Name                    | Calls Method                                                 |  Energy | Gradient | Reference              |
    +=========================+==============================================================+=========+==========+========================+
    | conv-mp2                | MP2                                                          |    Y    |     Y    | RHF/UHF                |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | omp2                    | Orbital-Optimized MP2                                        |    Y    |     Y    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | scs-omp2                | Spin-Component Scaled Orbital-Optimized MP2                  |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | sos-omp2                | Spin-Opposite Scaled Orbital-Optimized MP2                   |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | scsn-omp2               | A special version of SCS-OMP2 for nucleobase interactions    |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | scs-mi-omp2             | A special version of SCS-OMP2 (from S22 database)            |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | scs-omp2-vdw            | A special version of SCS-OMP2 (from ethene dimers)           |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | sos-pi-omp2             | A special version of SOS-OMP2 for :math:`\pi`-systems        |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+

The orbital-optimized MP3 methods currently supported in |Psifour| are outlined in Table :ref:`OMP3 Methods <table:omp3_calls>`.

    .. _`table:omp3_calls`:

    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | Name                    | Calls Method                                                 |  Energy | Gradient | Reference              |
    +=========================+==============================================================+=========+==========+========================+
    | mp3                     | MP3                                                          |    Y    |     Y    | RHF/UHF                |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | omp3                    | Orbital-Optimized MP3                                        |    Y    |     Y    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | scs-omp3                | Spin-Component Scaled Orbital-Optimized MP3                  |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | sos-omp3                | Spin-Opposite Scaled Orbital-Optimized MP3                   |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | scsn-omp3               | A special version of SCS-OMP3 for nucleobase interactions    |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | scs-mi-omp3             | A special version of SCS-OMP3 (from S22 database)            |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | scs-omp3-vdw            | A special version of SCS-OMP3 (from ethene dimers)           |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | sos-pi-omp3             | A special version of SOS-OMP3 for :math:`\pi`-systems        |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+

The orbital-optimized MP2.5 methods currently supported in |Psifour| are outlined in Table :ref:`OMP2.5 Methods <table:omp2_5_calls>`.

    .. _`table:omp2_5_calls`:

    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | Name                    | Calls Method                                                 |  Energy | Gradient | Reference              |
    +=========================+==============================================================+=========+==========+========================+
    | mp2.5                   | MP2.5                                                        |    Y    |     Y    | RHF/UHF                |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | omp2.5                  | Orbital-Optimized MP2.5                                      |    Y    |     Y    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+


The orbital-optimized CEPA methods currently supported in |Psifour| are outlined in Table :ref:`OCEPA Methods <table:ocepa_calls>`.

    .. _`table:ocepa_calls`:

    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | Name                    | Calls Method                                                 |  Energy | Gradient | Reference              |
    +=========================+==============================================================+=========+==========+========================+
    | ocepa                   | Orbital-Optimized CEPA                                       |    Y    |     Y    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | scs-ocepa               | Spin-Component Scaled Orbital-Optimized CEPA                 |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | sos-ocepa               | Spin-Opposite Scaled Orbital-Optimized CEPA                  |    Y    |     N    | RHF/ROHF/UHF/RKS/UKS   |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+
    | cepa0                   | CEPA(0) (identical to Linearized CCD)                        |    Y    |     Y    | RHF/UHF                |
    +-------------------------+--------------------------------------------------------------+---------+----------+------------------------+


.. index:: OMP2; setting keywords
.. index:: OMP3; setting keywords
.. index:: OMP2.5; setting keywords
.. index:: OCEPA; setting keywords

Basic Keywords
~~~~~~~~~~~~~~

.. include:: /autodir_options_c/occ__e_convergence.rst
.. include:: /autodir_options_c/occ__r_convergence.rst
.. include:: /autodir_options_c/occ__rms_mograd_convergence.rst
.. include:: /autodir_options_c/occ__max_mograd_convergence.rst
.. include:: /autodir_options_c/occ__mo_maxiter.rst
.. include:: /autodir_options_c/occ__wfn_type.rst
.. include:: /autodir_options_c/occ__orb_opt.rst

Advanced Keywords
~~~~~~~~~~~~~~~~~
.. include:: /autodir_options_c/occ__opt_method.rst
.. include:: /autodir_options_c/occ__mo_diis_num_vecs.rst
.. include:: /autodir_options_c/occ__lineq_solver.rst
.. include:: /autodir_options_c/occ__orth_type.rst
.. include:: /autodir_options_c/occ__mp2_os_scale.rst
.. include:: /autodir_options_c/occ__mp2_ss_scale.rst
.. include:: /autodir_options_c/occ__mp2_sos_scale.rst
.. include:: /autodir_options_c/occ__mp2_sos_scale2.rst
.. include:: /autodir_options_c/occ__nat_orbs.rst
.. include:: /autodir_options_c/occ__occ_orbs_print.rst
.. include:: /autodir_options_c/occ__tpdm_abcd_type.rst
.. include:: /autodir_options_c/occ__do_diis.rst

