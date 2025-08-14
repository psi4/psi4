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

.. image:: https://img.shields.io/badge/home-simint-5077AB.svg
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

* .. image:: https://anaconda.org/psi4/simint/badges/version.svg
     :target: https://anaconda.org/psi4/simint

* SIMINT is available as a conda package for Linux and macOS (and Windows, through the Ubuntu shell).

* The conda package is compiled to least-common-denominator, namely SSE instruction set.

* If using the |PSIfour| binary, simint has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  simint can be obtained through ``conda install simint -c psi4``.
  Then enable it as a feature with :makevar:`ENABLE_simint`,
  hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect simint and activate dependent code.

.. * Previous bullet had details. To build |PSIfour| from source and use 
..   simint from conda without thinking, consult.

* To remove a conda installation, ``conda remove simint``.

**Source**

* .. image:: https://img.shields.io/github/tag/psi4/simint.svg?maxAge=2592000

..     :target: https://github.com/psi4/simint TODO BPP

* If using |PSIfour| built from source and you want simint built from
  from source also,
  enable it as a feature with :makevar:`ENABLE_simint`,
  and let the build system fetch and build it and activate dependent code.


.. _`cmake:simint`:

How to configure OpenOrbitalOptimizer for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, OpenOrbitalOptimizer is a library that provides alternate
  orbital optimization.

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) OpenOrbitalOptimizer

* Upstream Dependencies |w---w| OpenOrbitalOptimizer |dr| Armadillo

**CMake Variables**

* :makevar:`ENABLE_simint` |w---w| CMake variable toggling whether |PSIfour| builds with simint
* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For simint, set to an installation directory containing ``include/simint/simint.h``
* :makevar:`simint_DIR` |w---w| CMake variable to specify where pre-built simint can be found. Set to installation directory containing ``share/cmake/simint/simintConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_simint` |w---w| CMake variable to force internal build of simint instead of detecting pre-built
* :makevar:`CMAKE_INSIST_FIND_PACKAGE_simint` |w---w| CMake variable to force detecting pre-built simint and not falling back on internal build
* :makevar:`SIMINT_VECTOR` |w---w| CMake variable for simint vectorization (i.e., scalar sse avx avxfma micavx512). Default is ``avx``, *not* detected, so ``sse`` may be required for older chipsets. See http://www.bennyp.org/research/simint/README.txt for details.

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake -DENABLE_simint=ON

B. Build *without* simint

  .. code-block:: bash

    >>> cmake

C. Link against pre-built

  .. code-block:: bash

    >>> cmake -DENABLE_simint=ON -DCMAKE_PREFIX_PATH=/path/to/simint/root

  .. code-block:: bash

    >>> cmake -DENABLE_simint=ON -Dsimint_DIR=/path/to/simint/configdir

D. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DENABLE_simint=ON -DCMAKE_PREFIX_PATH=/path/to/unwanted/simint/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_simint=ON

