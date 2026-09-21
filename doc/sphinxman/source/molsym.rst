.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2026 The Psi4 Developers.
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

.. index:: MolSym, SALC, symmetry

.. _`sec:molsym`:

Interface to MolSym by S. M. Goodlett and N. L. Kitzmiller
============================================================

.. codeauthor:: Stephen M. Goodlett, Nathaniel L. Kitzmiller
.. sectionauthor:: Nathaniel L. Kitzmiller

*Module:* :ref:`Keywords <apdx:findif>`

.. image:: https://img.shields.io/badge/home-molsym-informational.svg
   :target: https://github.com/NASymmetry/MolSym

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: https://molsym.readthedocs.io/en/latest/

|PSIfour| contains code to interface to MolSym, a molecular point-group detection
and symmetry-adapted linear combination (SALC) generation package developed by
S. M. Goodlett and N. L. Kitzmiller. Within |PSIfour|, MolSym is an alternative to
the internal ``CdSalc`` machinery for building the Cartesian SALCs used by the
finite-difference driver, and additionally exposes exploitation of non-Abelian
point-group degeneracy for degenerate irreps. For the underlying theory and
MolSym's standalone API (point-group detection, SALC construction for arbitrary
function sets), see the `MolSym documentation <https://molsym.readthedocs.io/en/latest/>`_.

No additional licence or configuration is required to use MolSym with |PSIfour|.

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/conda-forge/molsym/badges/version.svg
     :target: https://anaconda.org/conda-forge/molsym

* .. image:: https://img.shields.io/pypi/v/molsym
     :target: https://pypi.org/project/molsym

* MolSym is available for Linux and macOS as the ``molsym`` package on
  conda-forge and on PyPI.

* To install from conda run ``conda install molsym -c conda-forge``.

* To remove a conda installation, ``conda remove molsym``.

**Source**

* .. image:: https://img.shields.io/github/tag-date/NASymmetry/MolSym.svg?maxAge=2592000
     :target: https://github.com/NASymmetry/MolSym

* If using |PSIfour| built from source and you want MolSym installed as well,
  enable it as a feature with :makevar:`ENABLE_molsym`,
  and let the build system fetch and install it.

.. _`sec:usingMolSym`:

Using MolSym for finite-difference SALCs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, |PSIfour|'s finite-difference driver builds Cartesian SALCs using
its internal ``CdSalc`` code. Setting |findif__salc_package| to ``MOLSYM``
switches SALC generation (and the projection of translational/rotational
Eckart conditions) over to MolSym instead. This does not change the
displacements performed or how derivatives are assembled from them, only how
the symmetry-adapted displacement coordinates are constructed.

When MolSym is selected, |findif__molsym_exploit_degeneracy| additionally
controls whether |PSIfour| takes advantage of MolSym's non-Abelian point-group
handling to skip displacements along SALCs that are related by symmetry to
ones already displaced, generating the corresponding results by
transformation rather than by additional finite-difference steps.

A minimal input requesting a MolSym-based frequency analysis would look like
the following: ::

    molecule nh3 {
        N  0.          0.          0.07056746
        H  0.46649474 -0.80799259 -0.32682968
        H  0.46649474  0.80799259 -0.32682968
        H -0.93298949  0.         -0.32682968
    }

    set {
        basis sto-3g
        scf_type pk
        salc_package molsym
        molsym_exploit_degeneracy true
    }

    freqs, wfn = frequencies('scf', return_wfn=True)

.. include:: autodir_options_c/findif__salc_package.rst
.. include:: autodir_options_c/findif__molsym_exploit_degeneracy.rst


.. _`cmake:molsym`:

How to configure MolSym for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, MolSym is a library for point-group detection and
  symmetry-adapted linear combination (SALC) generation, used as an optional
  backend for the finite-difference driver.

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) MolSym

* Upstream Dependencies |w---w| MolSym |dr| NumPy

**CMake Variables**

* :makevar:`ENABLE_molsym` |w---w| CMake variable toggling whether Psi4 automatically installs MolSym

**Examples**

A. Build and install MolSym if needed

  .. code-block:: bash

    >>> cmake -DENABLE_molsym=ON

B. Build *without* MolSym

  .. code-block:: bash

    >>> cmake
