
Theory, Usage, and Notes
------------------------

Theories of the coupled-pair type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. codeauthor:: A. Eugene DePrince
.. sectionauthor:: A. Eugene DePrince

Coupled-pair theories.

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

Note that we have included the path to the plugin directory (here, /Users/deprince/psi4/tests/) and imported the plugin. These commands are necesarry to call the CEPA(1) procedure via :py:func:`~driver.energy`. 

Coupled-pair theories
^^^^^^^^^^^^^^^^^^^^^

With minor modifications to the working equations, a family of coupled-pair
methods may be implemented.  The current capabilities include:

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

