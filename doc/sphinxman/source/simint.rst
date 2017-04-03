.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2017 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This program is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU General Public License as published by
.. # the Free Software Foundation; either version 2 of the License, or
.. # (at your option) any later version.
.. #
.. # This program is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU General Public License for more details.
.. #
.. # You should have received a copy of the GNU General Public License along
.. # with this program; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #

.. include:: autodoc_abbr_options_c.rst

.. index:: SIMINT, integrals

.. _`sec:simint`:

Interface to SIMINT by B. Pritchard
===================================

.. codeauthor:: Benjamin P. Pritchard
.. sectionauthor:: Lori A. Burns

.. *Module:* :ref:`Keywords <apdx:efp>`, :ref:`PSI Variables <apdx:efp_psivar>`, :source:`LIBEFP <src/lib/libefp_solver>`

.. image:: https://img.shields.io/badge/home-simint-5077AB.svg
   :target: http://www.bennyp.org/research/simint/

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: http://www.bennyp.org/research/simint/

These are the vectorized implementation of the Obara-Saika (OS) method of
calculating electron repulsion integrals developed by B. Pritchard and
interfaced into libmints. Enabling simint and adding ``set
integral_package simint`` (do this in ``~/.psi4rc`` for universal effect)
runs libderiv from libint for derivative integrals and simint for
non-derivative integrals. Note that present AM maximum is ``$$(gg|gg)$$``

Installation
~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/simint/badges/version.svg
     :target: https://anaconda.org/psi4/simint

* SIMINT is available as a conda package for Linux and macOS.

* The conda package is compiled to least-common-denominator, namely SSE instruction set.

.. * If using the |PSIfour| binary, simint has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  simint can be obtained through ``conda install simint``.
  Then enable it as a feature with :makevar:`ENABLE_simint`,
  hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect simint and activate dependent code.

* To remove a conda installation, ``conda remove simint``.

**Source**

* .. image:: https://img.shields.io/github/tag/psi4/simint.svg?maxAge=2592000
..     :target: https://github.com/psi4/simint TODO BPP

* If using |PSIfour| built from source and you want simint built from
  from source also,
  enable it as a feature with :makevar:`ENABLE_simint`,
  and let the build system fetch and build it and activate dependent code.


.. _`cmake:simint`:

How to configure simint for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, simint is a library that provides alternate
  integrals.

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) simint

* Upstream Dependencies |w---w| simint |dr| None

**CMake Variables**

* :makevar:`ENABLE_simint` |w---w| CMake variable toggling whether |PSIfour| builds with simint
* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For simint, set to an installation directory containing ``include/simint/simint.h``
* :makevar:`simint_DIR` |w---w| CMake variable to specify where pre-built simint can be found. Set to installation directory containing ``share/cmake/simint/simintConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_simint` |w---w| CMake variable to force internal build of simint instead of detecting pre-built
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

