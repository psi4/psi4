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

.. index:: DKH
.. _`sec:DKH`:

Interface to DKH by A. Wolf, M. Reiher, and B. A. Hess
======================================================

.. codeauthor:: Justin M. Turney
.. sectionauthor:: Justin M. Turney

*Module:* :ref:`Keywords <apdx:dkh>`, :ref:`Samples <apdx:testSuitedkh>`

.. image:: https://img.shields.io/badge/home-dkh-5077AB.svg
   :target: https://github.com/psi4/dkh

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: http://www.reiher.ethz.ch/software/dkh-x2c.html

.. _`sec:dkhinstall`:

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/dkh/badges/version.svg
     :target: https://anaconda.org/psi4/dkh

* DKH is available as a conda package for Linux and macOS.

* If using the |PSIfour| binary, dkh has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  dkh can be obtained through ``conda install dkh``.
  Then enable it as a feature with :makevar:`ENABLE_dkh`,
  hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect dkh and activate dependent code.

* To remove a conda installation, ``conda remove dkh``.

**Source**

* .. image:: https://img.shields.io/github/tag/psi4/dkh.svg?maxAge=2592000
     :target: https://github.com/psi4/dkh

* If using |PSIfour| built from source and you want dkh built from
  from source also,
  enable it as a feature with :makevar:`ENABLE_dkh`,
  and let the build system fetch and build it and activate dependent code.

.. _`sec:dkhinput`:

Input
~~~~~

For all electron calculations one can use the Douglas-Kroll-Hess (DKH)
Hamiltonian to take into account scalar relativistic effects.

Minimal input for DKH single-point computation looks like this::

    molecule {
    Mg
    }

    set basis aug-cc-pvdz-dk
    set relativistic dkh

    energy('scf')

By default a 2nd-order DKH calculation is performed. To change the default
order use the |globals__dkh_order| option. The version of the code found in
|Psifour| is capable of up to 4th-order DKH calculations.

Keywords
~~~~~~~~

.. include:: autodir_options_c/globals__relativistic.rst
.. include:: autodir_options_c/globals__dkh_order.rst

Reference
~~~~~~~~~

When using this code please make reference to the appropriate following paper:

* "The Generalized Douglas-Kroll Transformation," A. Wolf,
  M. Reiher, and B. A. Hess, *J. Chem. Phys.* **117**, 9215 (2002).
  (doi: `10.1063/1.1515314 <http://dx.doi.org/10.1063/1.1515314>`_)

.. _`cmake:dkh`:

How to configure dkh for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, DKH is a library that provides additional
  quantum chemical capabilities (relativistic effects).

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) dkh

* Upstream Dependencies |w---w| dkh |dr| Fortran

**CMake Variables**

* :makevar:`ENABLE_dkh` |w---w| CMake variable toggling whether Psi4 builds with dkh
* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For dkh, set to an installation directory containing ``include/DKH/DKH_MANGLE.h``
* :makevar:`dkh_DIR` |w---w| CMake variable to specify where pre-built dkh can be found. Set to installation directory containing ``share/cmake/dkh/dkhConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_dkh` |w---w| CMake variable to force internal build of dkh instead of detecting pre-built

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake -DENABLE_dkh=ON

B. Build *without* dkh

  .. code-block:: bash

    >>> cmake

C. Link against pre-built

  .. code-block:: bash

    >>> cmake -DENABLE_dkh=ON -DCMAKE_PREFIX_PATH=/path/to/dkh/root

  .. code-block:: bash

    >>> cmake -DENABLE_dkh=ON -Ddkh_DIR=/path/to/dkh/configdir

D. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DENABLE_dkh=ON -DCMAKE_PREFIX_PATH=/path/to/unwanted/dkh/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_dkh=ON

