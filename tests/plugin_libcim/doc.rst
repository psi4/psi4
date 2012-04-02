
Theory, Usage, and Notes
------------------------

Cluster-in-molecule local correlation calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. codeauthor:: A. Eugene DePrince
.. sectionauthor:: A. Eugene DePrince

The cluster-in-molecule (CIM) local correlation framework decomposes the correlation 
energy into contributions from individual occupied orbitals.  Within a local molecular
orbital (LMO) basis, the contribution to the correlation energy from LMO, i, can
usually be determined quite accurately even while ignoring a large subset of occupied
LMOs that interact only weakly with LMO, i.  A simplistic CIM procedure would be one 
in which a correlated calculation was performed to obtain a correlation energy for 
each occupied LMO, i, where only those occupied orbitals, j, that interact strongly
with i are active in the calculation.  In general, there are far fewer occupied 
orbitals in each subcalculation than in the original molecule, leading to a drastic
reduction in computational cost.

Example input
^^^^^^^^^^^^^^^^

The following is a simple input that will perform a CIM-CCSD(T) computation on a 
cluster of 8 water molecules. ::

	sys.path.insert(0, '/Users/deprince/psi4/tests')
	import plugin_libcim
	molecule water_cluster {
		0 1
		O      -1.0651248364       -3.8315679992        0.7532888494
		O      -2.9143323066        0.0506352184        2.6495240589
		O       3.1809638699       -1.9024798828        1.3330724074
		O       1.3354866061        1.9547090768        3.2830335354
		O      -3.2861333963        1.9481142161       -1.3376469592
		O       1.0275334495        3.7242249597       -0.7678751934
		O       3.0157785303       -0.0712267216       -2.6986679801
		O      -1.2941710638       -1.8723989728       -3.2147288750
		H      -1.9357609046       -2.9760352188       -4.5487433902
		H      -1.2742150895       -2.9485799993       -1.5252679877
		H      -2.0054584299       -2.5457669126        1.7642905783
		H       0.7277762380       -3.4215608192        1.1740336558
		H       4.7885807922       -2.7315714564        1.7037902075
		H       3.3238228662       -1.1521031580       -0.5188299538
		H       2.4484094673        1.6594111405       -2.2060663486
		H       1.4072791668       -0.8502101370       -3.3032967324
		H      -4.4272510754        0.2214890713        3.6940950394
		H      -3.2975985621        0.9695995084        0.9113139951
		H      -2.7572072127        0.4842512058       -2.4039146177
		H      -1.7160495127        2.9938333660       -1.3066841156
		H       1.5744710900        5.4858950725       -0.8492536409
		H       1.2479856578        3.1310886550        1.1327858415
		H       2.3141902547        0.4021201957        2.8456440982
		H      -0.4189882778        1.2779824493        3.4361058526
		units bohr
	}
	memory 2000 mb
	set {
		basis 6-31g
		freeze_core true
	}
	energy('cim-ccsd(t)')

Note that we have included the path to the plugin directory (here, /Users/deprince/psi4/tests/) and imported the plugin. These commands are necesarry to call the CIM procedure via the :py:func:`~driver.energy`. Publications from the use of the CIM code should cite my non-existent paper: [DePrince!].

CIM
^^^^^^^^^^^^^^^^

In principle, the CIM framework can be
applied to any post-Hartree-Fock method. Current capabilities include
CCSD, CCSD(T), and a family of coupled-pair methods:

* cim-ccsd, coupled cluster with single and double excitations (CCSD)
* cim-ccsd(t), CCSD with perturbative triples
* cim-cepa(n), coupled electron pair approximation (CEPA), variant n
* cim-acpf, averaged coupled pair functional (ACPF)
* cim-aqcc, averaged quadratic coupled cluster (AQCC)
* cim-cisd, configuration interaction with single and double excitations (CISD)

Basic CIM Keywords
~~~~~~~~~~~~~~~~~~

.. include:: /autodir_options_c/mints__basis.rst
.. include:: /autodir_options_c/dfmp2__df_basis_mp2.rst
.. include:: /autodir_options_c/globals__freeze_core.rst
.. include:: /autodir_plugins/plugin_libcim__cim_domain_type.rst
.. include:: /autodir_plugins/plugin_libcim__cim_se_tolerance.rst
.. include:: /autodir_plugins/plugin_libcim__cim_de_tolerance1.rst
.. include:: /autodir_plugins/plugin_libcim__cim_de_tolerance2.rst
.. include:: /autodir_plugins/plugin_libcim__occ_tolerance.rst

Advanced CIM Keywords
~~~~~~~~~~~~~~~~~~~~~

.. include:: /autodir_plugins/plugin_libcim__denominator_delta.rst
.. include:: /autodir_plugins/plugin_libcim__boys_convergence.rst
.. include:: /autodir_plugins/plugin_libcim__boys_maxiter.rst
.. include:: /autodir_plugins/plugin_libcim__cim_initialize.rst
.. include:: /autodir_plugins/plugin_libcim__cim_cluster_num.rst

Running CIM in parallel
^^^^^^^^^^^^^^^^^^^^^^^

One desirable feature of the CIM framework is the embarassingly parallel structure of the energy decomposition.  In |PSIfour|, one can perform the CIM procedure in parallel (with mpi4py) on N nodes by ``mpirun -n N psi4 input.dat``.  Note that threading within the underlying correlated calculations is independent of the MPI parallelism.  The command ``mpirun -n 8 psi4 input.dat`` with :term:`num_threads` set to 8 in input.dat will spawn 8 mpi ranks, and each rank will use 8 OpenMP threads within the correlated calculation.
