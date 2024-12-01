.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2024 The Psi4 Developers.
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

.. index:: CPPE, PE

.. _`sec:cppe`:

Interface to CPPE by M. Scheurer
=======================================

.. codeauthor:: Maximilian Scheurer
.. sectionauthor:: Maximilian Scheurer

*Module:* :ref:`Keywords <apdx:pe>`, :ref:`PSI Variables <apdx:pe_psivar>`

.. image:: https://img.shields.io/badge/home-cppe-informational.svg
   :target: https://github.com/maxscheurer/cppe

.. .. raw:: html
.. 
..    <br>
.. 
.. .. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
..    :target: http://pcmsolver.readthedocs.io/en/latest/

|PSIfour| contains code to interface to the CPPE library developed
by M. Scheurer.
The CPPE library requires no additional licence, downloads, or
configuration. The library allows for calculations in solution with the
polarizable embedding model (PE), an explicit, fragment-based solvent model [Olsen:2010:3721]_.

For a general tutorial on how to prepare/perform PE calculations, read the
`tutorial review <https://onlinelibrary.wiley.com/doi/full/10.1002/qua.25717>`_.

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/cppe/badges/version.svg
     :target: https://anaconda.org/psi4/cppe

* CPPE is available as a conda package for Linux and macOS (and Windows, through the Ubuntu shell).

* If using the |PSIfour| binary, CPPE has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  CPPE can be obtained through ``conda install cppe -c psi4``.
  Then enable it as a feature with :makevar:`ENABLE_cppe`,
  hint its location with :makevar:`cppe_DIR`,
  and rebuild |PSIfour| to detect CPPE and activate dependent code.

.. * Previous bullet had details. To build |PSIfour| from source and use
..   cppe from conda without thinking, consult.

* To remove a conda installation, ``conda remove cppe``.

**Source**

* .. image:: https://img.shields.io/github/tag-date/maxscheurer/cppe.svg?maxAge=2592000
     :target: https://github.com/maxscheurer/cppe

* If using |PSIfour| built from source and you want CPPE built from
  from source also,
  enable it as a feature with :makevar:`ENABLE_cppe`,
  and let the build system fetch and build it and activate dependent code.

.. index:: PE; Using PE

.. _`sec:usingPE`:

Using the polarizable embedding model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The inclusion of a PE description of the solvent into your calculation
is achieved by setting |globals__pe| ``true`` in your input file.

.. note:: At present, PE can only be used for energy calculations with SCF
          wavefunctions and CC wavefunctions in the PTE approximation [Cammi:2009:164104]_.
          All ERI algorithms (``PK``, ``OUT_OF_CORE``, ``DIRECT``, ``DF``, ``CD``) are supported.

.. note:: linear response calculations (static polarisabilities, TD-SCF) are supported for RHF/UHF if available.

.. warning:: The CPPE library **cannot** exploit molecular point group symmetry.

.. .. warning:: Analytic gradients and Hessians **are not** available with PE. Finite differences will be used
..              regardless of the ``dertype`` passed to the ``optimize`` function.
..              See :srcsample:`pcmsolver/opt-fd` for a sample input.

.. The PCM model and molecular cavity are specified in a ``pcm`` section that has
.. to be explicitly typed in by the user. This additional section follows a syntax
.. that is slightly different from that of |Psifour| and is fully documented
.. `here <http://pcmsolver.readthedocs.io/en/latest/users/input.html>`_

A typical input for a Hartree--Fock calculation with PE would look like the following: ::

    molecule pna {
        C          8.64800        1.07500       -1.71100
        C          9.48200        0.43000       -0.80800
        C          9.39600        0.75000        0.53800
        C          8.48200        1.71200        0.99500
        C          7.65300        2.34500        0.05500
        C          7.73200        2.03100       -1.29200
        H         10.18300       -0.30900       -1.16400
        H         10.04400        0.25200        1.24700
        H          6.94200        3.08900        0.38900
        H          7.09700        2.51500       -2.01800
        N          8.40100        2.02500        2.32500
        N          8.73400        0.74100       -3.12900
        O          7.98000        1.33100       -3.90100
        O          9.55600       -0.11000       -3.46600
        H          7.74900        2.71100        2.65200
        H          8.99100        1.57500        2.99500
        symmetry c1
        no_reorient
        no_com
    }

    set {
     basis sto-3g
     pe true
     e_convergence 10
     d_convergence 10
     scf_type pk
    }

    set pe {
     potfile pna_6w.pot
    }

    scf_energy, wfn = energy('scf', return_wfn=True)


The corresponding potential file `pna_6w.pot` can be downloaded
`here <https://raw.githubusercontent.com/maxscheurer/cppe/master/tests/potfiles/pna_6w.pot>`_.

Keywords for CPPE
~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/globals__pe.rst
.. include:: autodir_options_c/pe__potfile.rst
.. include:: autodir_options_c/pe__isotropic_pol.rst
.. include:: autodir_options_c/pe__induced_convergence.rst
.. include:: autodir_options_c/pe__maxiter.rst
.. include:: autodir_options_c/pe__border.rst
.. include:: autodir_options_c/pe__border_type.rst
.. include:: autodir_options_c/pe__border_n_redist.rst
.. include:: autodir_options_c/pe__border_redist_order.rst
.. include:: autodir_options_c/pe__border_rmin.rst
.. include:: autodir_options_c/pe__border_rmin_unit.rst
.. include:: autodir_options_c/pe__border_redist_pol.rst


.. _`cmake:cppe`:

How to configure CPPE for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, CPPE is a library that provides additional
  quantum chemical capabilities (explicit solvation modeling).

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) CPPE

* Upstream Dependencies |w---w| CPPE

**CMake Variables**

* :makevar:`ENABLE_cppe` |w---w| CMake variable toggling whether Psi4 builds with CPPE
* :makevar:`cppe_DIR` |w---w| CMake variable to specify where pre-built CPPE can be found. Set to installation directory containing ``share/cmake/cppe/cppeConfig.cmake``

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake -DENABLE_cppe=ON

B. Build *without* CPPE

  .. code-block:: bash

    >>> cmake

