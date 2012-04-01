
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
LMOs that interact only weakly with LMO, i.  In principle, the CIM framework can be
applied to any post-Hartree-Fock method, but this implementation is currently 
limited to the CCSD and CCSD(T) methods.

A First Example
^^^^^^^^^^^^^^^

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
		df_basis_scf cc-pvdz
		freeze_core true
	}
	energy('cim-ccsd(t)')

Note that we have included the path to the plugin directory (here, /Users/deprince/psi4/tests/)
and imported the plugin.  These commands are necesarry to call the CIM procedure via the 
:py:func:`~driver.energy` function.  Note also that we have specified a df basis even though this is not
a df computation.  The construction of the virtual orbitals in this CIM implementation uses
scaled opposite-spin MP2 natural orbitals, the construction of which require df technology.
Ideally, one should select a df basis that corresponds to the full basis, but, when no
corresponding basis exists, I just use the smallest basis available in |PSIfour|, cc-pVDZ.
Publications from the use of the CIM code should cite my non-existent paper: [DePrince!]_.

Basic CIM Keywords
~~~~~~~~~~~~~~~~~~

:term:`BASIS`
"""""""""""""

      Primary basis set, describes the molecular orbitals

      * **Type**: string
      * **Possible Values**: :ref:`basis string <apdx:basisElement>`
      * **Default**: No Default

:term:`DF_BASIS_SCF`
""""""""""""""""""""

      Auxilliary basis set, used to appriximate MP2 natural orbitals which are 
      used to define the virtual space for each cluster

      * **Type**: string
      * **Possible Values**: :ref:`basis string <apdx:basisElement>`
      * **Default**: No Default

:term:`FREEZE_CORE`
"""""""""""""""""""

      The scope of core orbitals to freeze in occupied orbital localization and 
      determination of occupied domains for the CIM procedure. 
      Recommended true for all CIM computations

      * **Type**: string
      * **Possible Values**: FALSE, TRUE, SMALL, LARGE
      * **Default**: FALSE

:term:`CIM_DOMAIN_TYPE`
"""""""""""""""""""""""

      CIM central domain type (dual- or single-environment CIM).  Recommeneded SECIM for all CIM computations.

      * **Type**: string
      * **Possible Values**: SECIM, DECIM
      * **Default**: SECIM

:term:`CIM_SE_TOLERANCE`
""""""""""""""""""""""""

      For a given occupied LMO, i, the minimum absolute value of the Fock matrix element, :math:`F_{ij}`, for occupied LMO, j, to be in i's cluster.  Only applies if CIM_DOMAIN_TYPE is SECIM

      * **Type**: double
      * **Default**: 0.001

:term:`CIM_DE_TOLERANCE1`
"""""""""""""""""""""""""

      For a given occupied LMO, i, the minimum absolute value of the Fock matrix element, :math:`F_{ij}`, for occupied LMO, j, to be included in the MO domain of LMO, i.
      Only applies if CIM_DOMAIN_TYPE is DECIM

      * **Type**: double
      * **Default**: 0.01

:term:`CIM_DE_TOLERANCE2`
"""""""""""""""""""""""""

      For a given occupied LMO, i, the minimum absolute value of the Fock matrix element, :math:`F_{ij}`, for occupied LMO, j, to be included in the environmental domain of LMO, i.
      Only applies if CIM_DOMAIN_TYPE is DECIM

      * **Type**: double
      * **Default**: 0.05

:term:`OCC_TOLERANCE`
"""""""""""""""""""""

      Minimum occupation (eigenvalues of the MP2 OPDM) below which virtual natural orbitals are discarded for for a given cluster

      * **Type**: double
      * **Default**: 5.0e-5


Advanced CIM Keywords
~~~~~~~~~~~~~~~~~~~~~

:term:`DENOMINATOR_DELTA`
"""""""""""""""""""""""""

      Maximum error allowed (Max error norm in Delta tensor) in the approximate energy denominators employed in evaluating the scaled opposite-spin MP2 OPDM used in defining the virtual orbitals for each CIM cluster.  The default may be more conservative than is necessary in practice.

      * **Type**: double
      * **Default**: 1.0e-6

:term:`BOYS_CONVERGENCE`
""""""""""""""""""""""""

      Convergence threshold for the localization procedure

      * **Type**: double
      * **Default**:  1.0e-6

:term:`BOYS_MAXITER`
""""""""""""""""""""

      Maximum number of iterations to converge the orbital localization procedure

      * **Type**: integer
      * **Default**: 100

:term:`CIM_INITIALIZE`
""""""""""""""""""""""

      Should the CIM procedure return after the occupied domains are determined?  This parameter is used internally by the python driver if the calculation is going to be run in parallel.  Changing this won't have any effect on the procedure

      * **Type**: bool
      * **Default**: False

:term:`CIM_CLUSTER_NUM`
"""""""""""""""""""""""

      For which cluster number should |PSIfour| evaluate the correlation energy?  This parameter is used internally by the python driver if the calculation is going to be run in parallel.  Changing this won't have any effect on the procedure

      * **Type**: integer
      * **Default**: 0


