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
Currently the conductor-like screening model (COSMO) [Klamt:1993:799]_ [Lipparini:2014:184108]_,
the polarisable continuum model (PCM) [Tomasi:2005:2999]_ [Nottoli:2019:6061]_
and the linearized poisson-boltzmann model (LPB) [Lu:2008:973]_ [Jha:2023:104105]_ are supported.
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
For a concise introduction to the
`theory behind ddx <https://ddsolvation.github.io/ddX/md_docs_theory.html>`_
or further `literature references <https://ddsolvation.github.io/ddX/label_references.html>`_
see the ddx documentation.

The usage of ddx-based solvation models is enabled
by specifying |globals__ddx| ``true`` in your input file.
The solvation model itself is selected using the |ddx__ddx_model| parameter.
Additionally the definition of the solvent and solute cavity is required
and further parameters allow to influence details of discretisation,
numerical integration and iterative solvers,
see the next sections for details.

.. note:: At present PCM can only be used for energy calculations with SCF
          wavefunctions in the PTE approximation [Cammi:2009:164104]_.
          All ERI algorithms (``PK``, ``OUT_OF_CORE``, ``DIRECT``, ``DF``, ``CD``) are supported.

.. note:: linear response calculations (static polarisabilities, TD-SCF) are supported for RHF/UHF if available.

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
        "ddx_model":     "pcm",
        "ddx_solvent":   "water",
        "ddx_radii_set": "uff",
    })

    scf_e = psi4.energy('SCF')

Solvent model and solvent cavity definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Beyond setting |globals__ddx| to ``true`` and selecting
a solvent model using |ddx__ddx_model|,
the definition of the solvent is mandatory.
Regularly one might want to influence the setup of the solvent
cavity as well.

The **solvent** can be defined either by directly providing a dielectric
constant using |ddx__ddx_solvent_epsilon| or by looking up the dielectric
constant in from a solvent trivial name provided by |ddx__ddx_solvent|
(e.g. ``water``, ``ethanol``, ``cis-1,2-dimethylcyclohexane``).
By convention solvent names are all lowercase and use dashes (``-``) to separate
quantifiers like ``o``, ``n`` etc.
The full list understood by ddx can be obtained using ::

    import pyddx
    print(pyddx.data.solvent_epsilon.keys())

For when an LPB solvent model is selected (|ddx__ddx_model| is ``LPB``)
additionally the **Debye-Hückel parameter** |ddx__ddx_solvent_kappa| needs to be provided
(in units of inverse Bohr or inverse Angström, depending on the unit used to define
the molecular geometry). ``pyddx`` provides a handy utility function to compute
the Debye-Hückel parameter. For example ::

    import pyddx
    from qcelemental import constants

    list_of_ions = [(+1, 0.1), (-1, 0.1)]
    dielectric_constant = pyddx.data.solvent_epsilon["water"]
    temperature = 298.15  # Kelvin
    kappa_invbohr = pyddx.solvent_kappa(list_of_ions, dielectric_constant, temperature)
    kappa_invang = kappa_invbohr / constants.bohr2angstroms

computes the parameter (in inverse Angström) for a 0.1 mol/l solution of sodium
chloride in water, thus a solution woith 0.1 mol/l of a ``+1``-charged ion
and 0.1 mol/l of a ``-1``-charged ion.

The **cavity** in ddx is defined as a union of spheres around each atom.
Usually the spehere radii for each atom are determined using a standard
set of tabulated radii per atomic species, determined by the |ddx__ddx_radii_set| parameter.
Currently ``bondi`` [Bondi:1964:441]_ and ``uff`` [Rappe:1992:114]_
are supported for |ddx__ddx_radii_set| with ``uff`` selected by default.
These radius values are conventionally scaled by an additional factor before use,
conventionally 1.1 for ``uff`` and 1.2 for ``bondi``. Customisation of the scaling
is possible using the |ddx__ddx_radii_scaling| parameter.
A more fine-grained control over the sphere radii is available by explicitly providing
a list of radii (one per atom, exactly in the order of the input geometry)
using the |ddx__ddx_radii| parameter. Note that the same unit as for the molecular
input is expected for the radii.

.. include:: autodir_options_c/globals__ddx.rst
.. include:: autodir_options_c/ddx__ddx_model.rst
.. include:: autodir_options_c/ddx__ddx_radii.rst
.. include:: autodir_options_c/ddx__ddx_radii_scaling.rst
.. include:: autodir_options_c/ddx__ddx_radii_set.rst
.. include:: autodir_options_c/ddx__ddx_solvent_epsilon.rst
.. include:: autodir_options_c/ddx__ddx_solvent.rst
.. include:: autodir_options_c/ddx__ddx_solvent_kappa.rst

Numerical integration and discretisation parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These parameters can be altered to balance the cost and accuracy
of the implict description of the solvation.

|ddx__ddx_solute_radial_points| and |ddx__ddx_solute_spherical_points| influence
the accuracy of the numerical grid used to obtain the representation
of the electric potential / field of the solute density,
since a standard DFT integration grid is used to obtain these quantities.
In contrast to the integration of DFT quantities much lower accuracy
is required, such that for this step considerably smaller grids are employed.
If extremely high accuracy reference solutions are required, the DDX
DFT integration parameters might need to be increased, but this is rarely needed.

|ddx__ddx_lmax| and |ddx__ddx_n_lebedev| determine the accuracy of the computations
on the boundary of the spheres around each atom performed by DDX. |ddx__ddx_lmax|
determines the largest angular momentum of the spherical harmonics basis used
to discretise quantities on the atomic spheres and |ddx__ddx_n_lebedev| determines the
number of points of the Lebedev angular grid used for integration on the spheres.
|ddx__ddx_n_lebedev| should be chosen higher than |ddx__ddx_solute_spherical_points|
and the defaults are usually good.

.. include:: autodir_options_c/ddx__ddx_solute_radial_points.rst
.. include:: autodir_options_c/ddx__ddx_solute_spherical_points.rst
.. include:: autodir_options_c/ddx__ddx_lmax.rst
.. include:: autodir_options_c/ddx__ddx_n_lebedev.rst

Iterative solver parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

These parameters determine how the forward and adjoint linear systems
of the solvation model are solved. Usually these parameters do not need
to be changed. Occasionally |ddx__ddx_solvation_convergence| might need to be adapted,
e.g. if only a very crude or a highly accurate SCF solution is targeted.

.. include:: autodir_options_c/ddx__ddx_diis_max_vecs.rst
.. include:: autodir_options_c/ddx__ddx_maxiter.rst
.. include:: autodir_options_c/ddx__ddx_solvation_convergence.rst

Further keywords for ddx
~~~~~~~~~~~~~~~~~~~~~~~~

These parameter should rarely require changes.
In particular |ddx__ddx_eta|, |ddx__ddx_shift| and |ddx__ddx_logfile|
are expert parameters and should not be altered beyond debugging.

.. include:: autodir_options_c/ddx__ddx_eta.rst
.. include:: autodir_options_c/ddx__ddx_fmm_local_lmax.rst
.. include:: autodir_options_c/ddx__ddx_fmm_multipole_lmax.rst
.. include:: autodir_options_c/ddx__ddx_fmm.rst
.. include:: autodir_options_c/ddx__ddx_incore.rst
.. include:: autodir_options_c/ddx__ddx_logfile.rst
.. include:: autodir_options_c/ddx__ddx_shift.rst


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

