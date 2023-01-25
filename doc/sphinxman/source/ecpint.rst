.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2023 The Psi4 Developers.
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

.. index:: LibECPInt, ecpint
.. _`sec:ecpint`:

Interface to LibECPInt by R. Shaw
=================================

.. codeauthor:: Andrew C. Simmonett
.. sectionauthor:: Lori A. Burns

.. image:: https://img.shields.io/badge/home-LibECPInt-5077AB.svg
   :target: https://github.com/robashaw/libecpint

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: https://libecpint.readthedocs.io/en/latest/index.html


Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/libecpint/badges/version.svg
     :target: https://anaconda.org/psi4/libecpint

* .. image:: https://anaconda.org/conda-forge/libecpint/badges/version.svg
     :target: https://anaconda.org/conda-forge/libecpint

* LibECPInt is available as a conda package for Linux and macOS.

* If using the |PSIfour| binary, LibECPInt has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  LibECPInt can be obtained through ``conda install libecpint``.
  Then enable it as a feature with :makevar:`ENABLE_ecpint`,
  hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect LibECPInt and activate dependent code.

* To remove a conda installation, ``conda remove libecpint``.

**Source**

* .. image:: https://img.shields.io/github/tag/robashaw/libecpint.svg?maxAge=2592000
     :target: https://github.com/robashaw/libecpint

* If using |PSIfour| built from source and you want LibECPInt built from
  from source also,
  enable it as a feature with :makevar:`ENABLE_ecpint`,
  and let the build system fetch and build it and activate dependent code.


.. _`cmake:ecpint`:

How to configure LibECPInt for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, LibECPInt is a library that provides additional
  quantum chemical capabilities (ECP integrals).

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) LibECPInt

* Upstream Dependencies |w---w| LibECPInt |dr| None

**CMake Variables**

* :makevar:`ENABLE_ecpint` |w---w| CMake variable toggling whether Psi4 builds with LibECPInt
* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For LibECPInt, set to an installation directory containing ``include/libecpint.hpp``
* :makevar:`ecpint_DIR` |w---w| CMake variable to specify where pre-built LibECPInt can be found. Set to installation directory containing ``lib/cmake/ecpint/ecpint-config.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_ecpint` |w---w| CMake variable to force internal build of ecpint instead of detecting pre-built
* :makevar:`CMAKE_INSIST_FIND_PACKAGE_ecpint` |w---w| CMake variable to force detecting pre-built LibECPInt and not falling back on internal build

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake -DENABLE_ecpint=ON

B. Build *without* LibECPInt

  .. code-block:: bash

    >>> cmake

C. Link against pre-built

  .. code-block:: bash

    >>> cmake -DENABLE_ecpint=ON -DCMAKE_PREFIX_PATH=/path/to/ecpint/root

  .. code-block:: bash

    >>> cmake -DENABLE_ecpint=ON -Decpint_DIR=/path/to/ecpint/configdir

D. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DENABLE_ecpint=ON -DCMAKE_PREFIX_PATH=/path/to/unwanted/ecpint/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_ecpint=ON

