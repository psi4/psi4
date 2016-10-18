.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2016 The Psi4 Developers.
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

.. :target: http://sebwouters.github.io/CheMPS2/index.html

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

