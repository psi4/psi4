
Theory, Usage, and Notes
------------------------

Coupled-pair methods
^^^^^^^^^^^^^^^^^^^^^

.. codeauthor:: A. Eugene DePrince
.. sectionauthor:: A. Eugene DePrince

Coupled-pair methods can be viewed as approximations to coupled-cluster (CC) theory or as size-extensive modifications
of truncated configuration interaction (CI) theory.  The methods have the same complexity as CI
with single and double excitations (CISD), and solving the CISD or coupled-pair equations requires
about half the floating point operations required to solve the CC with singles and doubles (CCSD) equations.  CISD,
CCSD, and the coupled-pair methods discussed below all scale formally with the sixth power of system size.  For a
detailed discussion of the properties of various coupled-pair methods, see Ref. [Wennmohs:2008]_.

What follows is a very basic description of the practical differences in the equations that define each of the
coupled-pair methods implemented in |Psifour|.  We begin with the CISD wave function

.. math::
    :label: CIwfn

    | \Psi \rangle = | \Psi_0 \rangle + \sum_i^{occ} \sum_a^{vir} t_i^a | \Psi_i^a\rangle + \frac{1}{4}\sum_{ij}^{occ} \sum_{ab}^{vir} t_{ij}^{ab} | \Psi_{ij}^{ab}\rangle,

where we have chosen the intermediate normalization, :math:`\langle \Psi_0 | \Psi \rangle = 1`.
For a closed-shell reference, the CISD correlation energy is given by

.. math::
    :label: CIenergy
    
    E_c = \frac{1}{4}\sum_{ijab}\langle \Psi_{ij}^{ab} | \hat{H} - E_0 | \Psi \rangle,

and the amplitudes can be determined by the solution to the coupled set of eqations:

.. math::
    :label: CIeqns
    
    0   &= \langle \Psi_{ij}^{ab} | \hat{H} - E_0 - E_c | \Psi \rangle, \\
    0   &= \langle \Psi_{i}^{a} | \hat{H} - E_0 - E_c | \Psi \rangle.

The CISD method is not size-extensive, but this problem can be overcome by making very simple modifications to the amplitude
equations.  With malice and forethought, we replace the correlation energy, :math:`E_c`, with generalized shifts for
the doubles and singles equations, :math:`\Delta_{ij}` and :math:`\Delta_i`:

.. math::
    :label: CEPAeqns
    
    0   &= \langle \Psi_{ij}^{ab} | \hat{H} - E_0 - \Delta_{ij} | \Psi \rangle, \\
    0   &= \langle \Psi_{i}^{a} | \hat{H} - E_0 - \Delta_i | \Psi \rangle.

The shifts, :math:`\Delta_{ij}` and :math:`\Delta_i`, are
chosen to approximate (with varying degrees of accuracy) the effects of triple and quadruple excitations.  
The values for :math:`\Delta_{ij}` and :math:`\Delta_i`  used in several coupled-pair methods are given in Table 
:ref:`CEPA Shifts <table:cepa_shifts>`.  Note that these shifts are defined in a spin-free formalism 
for closed-shell references only.  

    .. _`table:cepa_shifts`:

    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | method                  | :math:`\Delta_{ij}`                                        |  :math:`\Delta_i`                            |
    +=========================+============================================================+==============================================+
    | sdci                    | :math:`E_c`                                                |  :math:`E_c`                                 |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | dci                     | :math:`E_c`                                                |  NA                                          |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | cepa(0)                 | 0                                                          |  0                                           |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | cepa(1)                 | :math:`\frac{1}{2}\sum_k(\epsilon_{ik}+\epsilon_{jk})`     | :math:`\sum_k \epsilon_{ik}`                 |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | cepa(2)                 | :math:`\epsilon_{ij}`                                      | :math:`\epsilon_{ii}`                        |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | cepa(3)                 | :math:`-\epsilon_{ij}+\sum_k(\epsilon_{ik}+\epsilon_{jk})` | :math:`-\epsilon_{ii}+2\sum_k \epsilon_{ik}` |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | acpf                    | :math:`\frac{2}{N} E_c`                                    | :math:`\frac{2}{N} E_c`                      |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+
    | aqcc                    | :math:`[1-\frac{(N-3)(N-2)}{N(N-1)}]E_c`                   | :math:`[1-\frac{(N-3)(N-2)}{N(N-1)}]E_c`     |
    +-------------------------+------------------------------------------------------------+----------------------------------------------+

The pair correlation energy, :math:`\epsilon_{ij}`, is simply a partial sum of the correlation energy.  In a spin-free formalism,
the pair energy is given by

.. math::
   :label: pair_energy

   \epsilon_{ij} = \sum_{ab} t_{ij}^{ab} (2 v_{ij}^{ab} - v_{ij}^{ba})

Methods whose shifts (:math:`\Delta_{ij}` and :math:`\Delta_i`) do not explicitly depend on orbitals :math:`i` or :math:`j` 
(CISD, CEPA(0), ACPF, and AQCC)
have solutions that render the energy stationary with respect variations in the amplitudes.  This convenient property allows
density matrices and 1-electron properties to be evaluated without any additional effort.  

Example input
^^^^^^^^^^^^^

The following is a minimal input that will perform a CEPA(1) computation on a 
water molecule described by the cc-pVDZ basis. ::

	sys.path.insert(0, '/Users/deprince/psi4/tests')
	import plugin_cepa
	molecule water_cluster {
		0 1
                O
                H 1 1.0
                H 1 1.0 2 104.5
	}
	set basis cc-pVDZ
	energy('cepa(1)')

Note that we have included the path to the plugin directory (here, /Users/deprince/psi4/tests/) and imported the plugin. These commands are necesarry to call the CEPA(1) procedure via :py:func:`~driver.energy`.  The coupled-pair methods currently supported in |Psifour| are outlined in Table :ref:`CEPA Methods <table:cepa_calls>`.

    .. _`table:cepa_calls`:

    +-------------------------+--------------------------------------------------------------+---------+-------------+------------------------+
    | name                    | calls method                                                 |  energy | derivatives | 1-electron properties  |
    +=========================+==============================================================+=========+=============+========================+
    | cepa(0)                 | coupled electron pair approximation, variant 0               |    Y    |     N       |          Y             |
    +-------------------------+--------------------------------------------------------------+---------+-------------+------------------------+
    | cepa(1)                 | coupled electron pair approximation, variant 1               |    Y    |     N       |          N             |
    +-------------------------+--------------------------------------------------------------+---------+-------------+------------------------+
    | cepa(2)                 | coupled electron pair approximation, variant 2               |    Y    |     N       |          N             |
    +-------------------------+--------------------------------------------------------------+---------+-------------+------------------------+
    | cepa(3)                 | coupled electron pair approximation, variant 3               |    Y    |     N       |          N             |
    +-------------------------+--------------------------------------------------------------+---------+-------------+------------------------+
    | acpf                    | averaged coupled-pair functional                             |    Y    |     N       |          Y             |
    +-------------------------+--------------------------------------------------------------+---------+-------------+------------------------+
    | aqcc                    | averaged quadratic coupled-cluster                           |    Y    |     N       |          Y             |
    +-------------------------+--------------------------------------------------------------+---------+-------------+------------------------+
    | sdci                    | configuration interaction with single and double excitations |    Y    |     N       |          Y             |
    +-------------------------+--------------------------------------------------------------+---------+-------------+------------------------+
    | dci                     | configuration interaction with double excitations            |    Y    |     N       |          Y             |
    +-------------------------+--------------------------------------------------------------+---------+-------------+------------------------+

Basic Keywords
~~~~~~~~~~~~~~~~~~

.. include:: /autodir_options_c/mints__basis.rst
.. include:: /autodir_options_c/globals__freeze_core.rst
.. include:: /autodir_plugins/plugin_cepa__r_convergence.rst
.. include:: /autodir_plugins/plugin_cepa__maxiter.rst
.. include:: /autodir_plugins/plugin_cepa__diis_max_vecs.rst
.. include:: /autodir_plugins/plugin_cepa__mp2_scale_os.rst
.. include:: /autodir_plugins/plugin_cepa__mp2_scale_ss.rst
.. include:: /autodir_plugins/plugin_cepa__dipmom.rst

Advanced Keywords
~~~~~~~~~~~~~~~~~~~~~
.. include:: /autodir_plugins/plugin_cepa__cepa_level.rst
.. include:: /autodir_plugins/plugin_cepa__scs_cepa.rst
.. include:: /autodir_plugins/plugin_cepa__cepa_scale_os.rst
.. include:: /autodir_plugins/plugin_cepa__cepa_scale_ss.rst
.. include:: /autodir_plugins/plugin_cepa__cepa_no_singles.rst


Bibliography
============

.. [Wennmohs:2008]
   F. Wennmohs and F. Neese, 
   *Chem. Phys. Lett.* **343**, 217-230 (2008).
