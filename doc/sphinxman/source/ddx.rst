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
For a concise introduction to the `theory behind ddx <https://ddsolvation.github.io/ddX/md_docs_theory.html>`
or further `literature references <https://ddsolvation.github.io/ddX/label_references.html>`
see the ddx documentation.

The usage of ddx-based solvation models is enabled
by specifying |globals__ddx| ``true`` in your input file.
The solvation model itself is selected using the |ddx__model| parameter.
Additionally the definition of the solvent and solute cavity is required
and further parameters allow to influence details of discretisation,
numerical integration and iterative solvers,
see the next sections for details.

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
        "ddx__model":     "pcm",
        "ddx__solvent":   "water",
        "ddx__radii_set": "uff",
    })

    scf_e = psi4.energy('SCF')

Solvent model and solvent cavity definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Beyond setting |globals__ddx| to ``true`` and selecting
a solvent model using |ddx__model|,
the definition of the solvent is mandatory.
Regularly one might want to influence the setup of the solvent
cavity as well.

The **solvent** can be defined either by directly providing a dielectric
constant using |ddx__solvent_epsilon| or by looking up the dielectric
constant in from a solvent trivial name provided by |ddx__solvent|
(e.g. ``water``, ``ethanol``, ``cis-1,2-dimethylcyclohexane``).
By convention solvent names are all lowercase and use dashes (``-``) to separate
quantifiers like ``o``, ``n`` etc.
The full list understood by ddx can be obtained using ::

    import pyddx
    print(pyddx.data.solvent_epsilon.keys())

The **cavity** in ddx is defined as a union of spheres around each atom.
Usually the spehere radii for each atom are determined using a standard
set of tabulated radii per atomic species, determined by the |ddx__radii_set| parameter.
Currently ``bondi`` and ``uff`` are supported for |ddx__radii_set|
with ``uff`` selected by default.
These radius values are conventionally scaled by an additional factor before use,
conventionally 1.1 for ``uff`` and 1.2 for ``bondi``. Customisation of the scaling
is possible using the |ddx__radii_scaling| parameter.
A more fine-grained control over the sphere radii is available by explicitly providing
a list of radii (one per atom, exactly in the order of the input geometry)
using the |ddx__radii| parameter.

.. include:: autodir_options_c/globals__ddx.rst
.. include:: autodir_options_c/ddx__model.rst
.. include:: autodir_options_c/ddx__radii.rst
.. include:: autodir_options_c/ddx__radii_scaling.rst
.. include:: autodir_options_c/ddx__radii_set.rst
.. include:: autodir_options_c/ddx__solvent_epsilon.rst
.. include:: autodir_options_c/ddx__solvent.rst

Numerical integration and discretisation parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These parameters can be altered to balance the cost and accuracy
of the implict description of the solvation.

|ddx__dft_radial_points| and |ddx__dft_spherical_points| influence
the accuracy of the numerical grid used to obtain the representation
of the electric field of the solute density
(for which a standard DFT integration grid is used).
Unlike the case for DFT grids a finer grid than the employed default is
rarely needed. In fact often even a coarser grid than the default
can be used without trading too much accuracy
(e.g. 50 radial and 230 spherical points).

|ddx__lmax| and |ddx__n_lebedev| determine the accuracy of the discretisation
on the boundary of the spheres around each atom. The defaults are usually good.

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

These parameter should rarely require changes.
In particular |ddx__eta|, |ddx__shift| and |ddx__logfile|
are expert parameters and should not be altered beyond debugging.

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

