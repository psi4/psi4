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

.. index:: OOO, OpenOrbitalOptimizer, integrals

.. _`sec:ooo`:

Interface to OpenOrbitalOptimizer by S. Lehtola
===============================================

.. codeauthor:: Susi Lehtola
.. sectionauthor:: Lori A. Burns

.. *Module:* :ref:`Keywords <apdx:efp>`, :ref:`PSI Variables <apdx:efp_psivar>`, :source:`LIBEFP <src/lib/libefp_solver>`

.. image:: https://img.shields.io/badge/home-OpenOrbitalOptimizer-5077AB.svg
   :target: https://github.com/SusiLehtola/OpenOrbitalOptimizer

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: https://github.com/SusiLehtola/OpenOrbitalOptimizer

This is a collection of SCF optimization techniques (DIIS, EDIIS, ADIIS, ODA)
developed by S. Lehtola and interfaced to |PSIfour| .
Enabling OOO and adding ``set
orbital_optimizer_package openorbitaloptimizer`` or ``set
orbital_optimizer_package ooo``
runs the SCF iterations steps from OOO. Set to ``internal`` to revoke.

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/SusiLehtola/OpenOrbitalOptimizer/badges/version.svg
     :target: https://anaconda.org/SusiLehtola/OpenOrbitalOptimizer

* OpenOrbitalOptimizer is available as a conda package for Linux, macOS, and Windows.

* The conda package is (for the C++ interface) a header-only library.

* If using the |PSIfour| binary, OpenOrbitalOptimizer has already been compiled in.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  OpenOrbitalOptimizer can be obtained through ``conda install openorbitaloptimizer -c conda-forge``.
  Then enable it as a feature with :makevar:`ENABLE_OpenOrbitalOptimizer`,
  hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect OpenOrbitalOptimizer and activate dependent code.

* To remove a conda installation, ``conda remove openorbitaloptimizer``.

**Source**

* .. image:: https://img.shields.io/github/tag/SusiLehtola/OpenOrbitalOptimizer.svg?maxAge=2592000
     :target: https://github.com/SusiLehtola/OpenOrbitalOptimizer

* If using |PSIfour| built from source and you want OpenOrbitalOptimizer built from
  from source also,
  enable it as a feature with :makevar:`ENABLE_OpenOrbitalOptimizer`,
  and let the build system fetch and build it and activate dependent code.


.. _`options:ooo`:

OpenOrbitalOptimizer Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OpenOrbitalOptimizer is ready to use for RHF and UHF SCF
iterations. (Other references will fall back on the internal
optimizer.) Note that you may need to back off on convergence
criteria. OpenOrbitalOptimizer for the SAD guess is available for
experimentation only.

.. include:: autodir_options_c/globals__orbital_optimizer_package.rst
.. include:: autodir_options_c/scf__maxiter.rst
.. include:: autodir_options_c/scf__ooo_print.rst
.. include:: autodir_options_c/scf__diis_max_vecs
.. include:: autodir_options_c/scf__scf_initial_start_diis_transition
.. include:: autodir_options_c/scf__scf_initial_finish_diis_transition
.. include:: autodir_options_c/scf__ooo_diis_restart_factor
.. include:: autodir_options_c/scf__ooo_optimal_damping_threshold
.. include:: autodir_options_c/scf__sad_orbital_optimizer_package


.. _`cmake:ooo`:

How to configure OpenOrbitalOptimizer for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, OpenOrbitalOptimizer is a library that provides alternate
  orbital optimization.

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) OpenOrbitalOptimizer

* Upstream Dependencies |w---w| OpenOrbitalOptimizer |dr| Armadillo

**CMake Variables**

* :makevar:`ENABLE_OpenOrbitalOptimizer` |w---w| CMake variable toggling whether |PSIfour| builds with OpenOrbitalOptimizer
* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For OOO, set to an installation directory containing ``include/openorbitaloptimizer/scfsolver.hpp``
* :makevar:`OpenOrbitalOptimizer_DIR` |w---w| CMake variable to specify where pre-built OpenOrbitalOptimizer can be found. Set to installation directory containing ``share/cmake/OpenOrbitalOptimizer/OpenOrbitalOptimizer/OpenOrbitalOptimizerConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_OpenOrbitalOptimizer` |w---w| CMake variable to force internal build of OpenOrbitalOptimizer instead of detecting pre-built
* :makevar:`CMAKE_INSIST_FIND_PACKAGE_OpenOrbitalOptimizer` |w---w| CMake variable to force detecting pre-built OpenOrbitalOptimizerand not falling back on internal build

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake -DENABLE_OpenOrbitalOptimizer=ON

B. Build *without* OpenOrbitalOptimizer

  .. code-block:: bash

    >>> cmake

C. Link against pre-built

  .. code-block:: bash

    >>> cmake -DENABLE_OpenOrbitalOptimizer=ON -DCMAKE_PREFIX_PATH=/path/to/OpenOrbitalOptimizer/root

  .. code-block:: bash

    >>> cmake -DENABLE_OpenOrbitalOptimizer=ON -DOpenOrbitalOptimizer_DIR=/path/to/ooo/configdir

D. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DENABLE_OpenOrbitalOptimizer=ON -DCMAKE_PREFIX_PATH=/path/to/unwanted/OpenOrbitalOptimizer/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_OpenOrbitalOptimizer=ON


