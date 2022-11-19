.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2022 The Psi4 Developers.
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

.. index:: ddx, COSMO, PCM, continuum solvation

.. _`sec:ddx`:

Interface to ddx by A. Mikhalev, A. Jha, M. Nottoli and M. F. Herbst
====================================================================

.. codeauthor:: Michael F. Herbst
.. sectionauthor:: Michael F. Herbst

*Module:* :ref:`Keywords <apdx:ddx>`, :ref:`PSI Variables <apdx:ddx_psivar>`

.. image:: https://img.shields.io/badge/home-ddx-informational.svg
   :target: https://github.com/ddsolvation/ddX

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: https://ddsolvation.github.io/ddX/

|PSIfour| contains code to interface to the ddx FORTRAN library developed
by A. Mikhalev *et. al.*. The library provides a linear-scaling implementation
of standard continuum solvation models using a domain-decomposition ansatz
[Cances:2013:054111]_ [Stamm:2016:054101]_.
Currently the conductor-like screening model (COSMO) [Barone:1998:1995]_
and the polarisable continuum model (PCM) [Tomasi:2005:2999]_
are supported.
No additional licence or configuration is required to use ddx with Psi4.


Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/conda-forge/pyddx/badges/version.svg
     :target: https://anaconda.org/conda-forge/pyddx

* .. image:: https://img.shields.io/pypi/v/pyddx
     :target: https://pypi.org/project/pyddx

* ddx is available for Linux and macOS in form of the ``pyddx`` package
  on conda-forge and on pypi.

* To install from conda run ``conda install pyddx -c conda-forge``.

* To remove a conda installation, ``conda remove pyddx``.

**Source**

* .. image:: https://img.shields.io/github/tag-date/ddsolvation/ddx.svg?maxAge=2592000
     :target: https://github.com/ddsolvation/ddx

* If using |PSIfour| built from source and you want ddx installed as well,
  enable it as a feature with :makevar:`ENABLE_ddx`,
  and let the build system fetch and install it.

.. _`sec:usingDDX`:

Using dd-based continum solvation models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In |PSIfour| two option to enable continuum solvation models
are currently implemented using either the PCMSolver or ddx package.
PCMSolver is based on a boundary-element discretisation [Cances:1998:309]_,
while ddx is based on a domain decomposition approach
[Cances:2013:054111]_ making it linear scaling.
For more details about PCMSolver see the :ref:`section on PCMsolver <sec:pcmsolver>`.

The usage of ddx-based solvation models is enabled
by specifying |globals__ddx| ``true`` in your input file.
The solvation model itself is selected using the |ddx__model| parameter.

|Psifour| understands a number of other options to configure
numerical details in the library (discretisation, iterative solver, cavity setup),
see the next section.

.. note:: At present dd-based solvation models
          can only be used for energy calculations with SCF
          wavefunctions. All ERI algorithms (``PK``, ``OUT_OF_CORE``, ``DIRECT``, ``DF``,
          ``CD``) are supported.

.. warning:: Currently the ddx interface **cannot** exploit molecular point group symmetry.

.. warning:: Analytic gradients and Hessians are currently **not available**
             with dd-based solvation models.

A minimal input for a Hartree--Fock calculation with dd-based PCM would look like
the following: ::

    import psi4
    nh3 = psi4.geometry("""
        N     -0.0000000001    -0.1040380466      0.0000000000
        H     -0.9015844116     0.4818470201     -1.5615900098
        H     -0.9015844116     0.4818470201      1.5615900098
        H      1.8031688251     0.4818470204      0.0000000000
        symmetry c1
        no_reorient
        no_com
        units bohr
        """)

    psi4.set_options({
        "basis": "sto-3g",
        "scf_type": "pk",
        "ddx": True,
    })

    psi4.set_module_options("ddx", {
        "model":     "pcm",
        "solvent":   "water",
        "radii_set": "uff",
    })

    scf_e = psi4.energy('SCF')


Cavity and solvent setup
~~~~~~~~~~~~~~~~~~~~~~~~

TODO

.. include:: autodir_options_c/globals__ddx.rst
.. include:: autodir_options_c/ddx__model.rst
.. include:: autodir_options_c/ddx__radii.rst
.. include:: autodir_options_c/ddx__radii_scaling.rst
.. include:: autodir_options_c/ddx__radii_set.rst
.. include:: autodir_options_c/ddx__solvent_epsilon.rst
.. include:: autodir_options_c/ddx__solvent.rst

Numerical integration and discretisation parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO

.. include:: autodir_options_c/ddx__dft_radial_points.rst
.. include:: autodir_options_c/ddx__dft_spherical_points.rst
.. include:: autodir_options_c/ddx__lmax.rst
.. include:: autodir_options_c/ddx__n_lebedev.rst

Iterative solver parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

These parameters determine how the forward and adjoint linear systems
of the solvation model are solved. Usually these parameters do not need
to be changed. Occasionally |ddx__solvation_convergence| might need to be adapted,
e.g. if only a very crude or a highly accurate SCF solution is targeted.

.. include:: autodir_options_c/ddx__diis_max_vecs.rst
.. include:: autodir_options_c/ddx__maxiter.rst
.. include:: autodir_options_c/ddx__solvation_convergence.rst

Further keywords for ddx
~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/ddx__eta.rst
.. include:: autodir_options_c/ddx__fmm_local_lmax.rst
.. include:: autodir_options_c/ddx__fmm_multipole_lmax.rst
.. include:: autodir_options_c/ddx__fmm.rst
.. include:: autodir_options_c/ddx__incore.rst
.. include:: autodir_options_c/ddx__logfile.rst
.. include:: autodir_options_c/ddx__shift.rst


.. _`cmake:ddx`:

How to configure ddx for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, ddx is a library for providing fast continuum
  solvation models.

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) ddx

* Upstream Dependencies |w---w| ddx |dr| Fortran

**CMake Variables**

* :makevar:`ENABLE_ddx` |w---w| CMake variable toggling whether Psi4 automatically installs ddx

**Examples**

A. Build and install ddx if needed

  .. code-block:: bash

    >>> cmake -DENABLE_ddx=ON

B. Build *without* ddx

  .. code-block:: bash

    >>> cmake

