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

.. index:: OTR, OpenTrustRegion

.. _`sec:otr`:

Interface to OpenTrustRegion by J. Greiner
==========================================

.. codeauthor:: Jonas Greiner
.. sectionauthor:: Lori A. Burns

.. image:: https://img.shields.io/badge/home-OpenTrustRegion-5077AB.svg
   :target: https://github.com/eriksen-lab/OpenTrustRegion

.. raw:: html

   <br>

OpenTrustRegion is a black-box second-order orbital optimizer developed by
J. Greiner in the Eriksen lab and interfaced to |PSIfour|. Rather than iterating the Fock matrix
to self-consistency, it minimizes the SCF energy directly with respect to the
orbital rotation parameters, using a trust-region method with the orbital
Hessian applied on the fly.

OpenTrustRegion serves as |PSIfour|'s *second-order* SCF optimizer, an alternative
to the internal SOSCF code. It is engaged the same way SOSCF always has been --
by setting |scf__soscf| -- and additionally
``set second_order_orbital_optimizer_package otr``. First-order DIIS iterations
run as usual until the orbital gradient falls below |scf__soscf_start_convergence|,
at which point OpenTrustRegion takes over and carries the calculation to
convergence. Set |globals__second_order_orbital_optimizer_package| to ``internal``
to revoke.

The package driving the *first-order* iterations is chosen separately with
|globals__orbital_optimizer_package| (``internal`` or ``ooo``), so the two
keywords combine freely.

Installation
~~~~~~~~~~~~

**Binary**

* A conda package for OpenTrustRegion is available. Obtain it through
  ``conda install opentrustregion -c conda-forge``, then enable it as a feature
  with :makevar:`ENABLE_OpenTrustRegion`, hint its location with
  :makevar:`CMAKE_PREFIX_PATH`, and rebuild |PSIfour| to detect OpenTrustRegion
  and activate dependent code.

**Source**

* .. image:: https://img.shields.io/github/tag/eriksen-lab/OpenTrustRegion.svg?maxAge=2592000
     :target: https://github.com/eriksen-lab/OpenTrustRegion

* If using |PSIfour| built from source and you want OpenTrustRegion built from
  source also,
  enable it as a feature with :makevar:`ENABLE_OpenTrustRegion=ON`,
  and let the build system fetch and build it and activate dependent code.
  Note that OpenTrustRegion is written in Fortran, so a Fortran compiler is
  required for the source build.


.. _`options:otr`:

OpenTrustRegion Options
~~~~~~~~~~~~~~~~~~~~~~~

OpenTrustRegion is ready to use for RHF, UHF, ROHF and the corresponding KS
references. Where one of its settings has a genuine |PSIfour| counterpart, that
existing keyword is used:

* |scf__soscf| turns second-order iterations on, and
  |scf__soscf_start_convergence| sets the orbital-gradient RMS at which the
  handoff from DIIS occurs.
* the convergence threshold on the RMS orbital gradient is taken as the tighter
  of |scf__e_convergence| and |scf__d_convergence|, then tightened by a further
  two digits. OpenTrustRegion stops as soon as its own gradient test is met,
  whereas DIIS tends to sail well past the requested threshold; the margin keeps
  response properties and analytic Hessians, which rely on that incidental
  tightness, in agreement with the internal optimizer.
* |scf__maxiter| caps the number of macro-iterations. The micro-iteration cap is
  deliberately *not* taken from |scf__soscf_max_iter|, whose default of 5 is far
  below the 50 OpenTrustRegion expects.
* |scf__soscf_print| turns on the micro-iteration detail, as it does for the
  internal code. |scf__otr_print| reaches the levels a boolean cannot -- silencing
  the solver entirely, or adding its linear-solver trace -- and takes precedence
  wherever it has been set.

The remaining settings have no |PSIfour| counterpart and so are exposed as
dedicated keywords, each defaulting to OpenTrustRegion's own default:
|scf__otr_subsystem_solver|, |scf__otr_n_random_trial_vectors|,
|scf__otr_jacobi_davidson_start|, |scf__otr_line_search|,
|scf__otr_start_trust_radius|, |scf__otr_global_red_factor|,
|scf__otr_local_red_factor|, |scf__otr_seed| and |scf__otr_print|.

Because OpenTrustRegion owns the iteration loop from the handoff onward,
|PSIfour| falls back on the internal second-order code whenever the computation
needs something applied per-iteration by the driver, or something
OpenTrustRegion cannot express:

* a ``CUHF`` reference, for which no orbital Hessian is implemented,
* meta-GGA and VV10 functionals,
* a GRAC-shifted potential, which is spliced into V_xc outside the kernel the
  Hessian differentiates,
* MOM (|scf__mom_start|) and fractional occupation (|scf__frac_start|), which
  change the occupation during the SCF,
* EFP, PCM, DDX and PE embedding, whose contributions are added to the Fock
  matrix by the Python driver on each iteration.

A note printed to the output file names whichever condition applied. Separately,
|PSIfour| declines SOSCF altogether for the semi-numerical exchange builds
(``DFDIRJ+COSX``, ``DFDIRJ+LINK``, ``DFDIRJ+SNLINK``), which cannot supply the
non-symmetric K matrices any second-order method needs; that check is not
specific to OpenTrustRegion.

A line naming the optimizers actually in force is printed above the iteration
table, for instance::

  The orbital optimizer module is Internal, second-order OpenTrustRegion

The same information is recorded on the wavefunction, keyed by role, so it can be
inspected without scraping the output file. It reports whichever package actually
ran, not the one requested, so a fallback is visible::

  >>> e, wfn = energy("scf", return_wfn=True)
  >>> wfn.module_roles()
  {'orbital_optimizer': 'internal', 'second_order_orbital_optimizer': 'opentrustregion'}

This is distinct from :py:meth:`~psi4.core.Wavefunction.module`, which names the
level-of-theory module responsible for the current energy.

Some further differences from the internal SOSCF code are worth knowing about:

* OpenTrustRegion optimizes the orbitals at a *fixed* occupation, whereas the
  internal optimizer re-runs the aufbau assignment every iteration. Coming in
  after DIIS has reached |scf__soscf_start_convergence| the occupation is
  usually settled, but if canonicalizing the converged orbitals changes it,
  |PSIfour| reconverges from the new occupation.
* Being a genuine minimizer rather than a stationary-point finder, it can walk
  off a saddle that DIIS would have stopped on. This is usually welcome, but it
  means a deliberately prepared excited or symmetry-broken state may collapse to
  the ground state unless the occupation is pinned.
* Stability analysis (|scf__stability_analysis|) runs through |PSIfour|'s own
  code as usual; OpenTrustRegion's internal stability check is left off so that
  both optimizers land on the same solution.

.. include:: autodir_options_c/globals__second_order_orbital_optimizer_package.rst
.. include:: autodir_options_c/globals__orbital_optimizer_package.rst
.. include:: autodir_options_c/scf__soscf.rst
.. include:: autodir_options_c/scf__soscf_start_convergence.rst
.. include:: autodir_options_c/scf__e_convergence.rst
.. include:: autodir_options_c/scf__d_convergence.rst
.. include:: autodir_options_c/scf__maxiter.rst
.. include:: autodir_options_c/scf__soscf_print.rst
.. include:: autodir_options_c/scf__otr_subsystem_solver.rst
.. include:: autodir_options_c/scf__otr_n_random_trial_vectors.rst
.. include:: autodir_options_c/scf__otr_jacobi_davidson_start.rst
.. include:: autodir_options_c/scf__otr_line_search.rst
.. include:: autodir_options_c/scf__otr_start_trust_radius.rst
.. include:: autodir_options_c/scf__otr_global_red_factor.rst
.. include:: autodir_options_c/scf__otr_local_red_factor.rst
.. include:: autodir_options_c/scf__otr_seed.rst
.. include:: autodir_options_c/scf__otr_print.rst


.. _`cmake:otr`:

How to configure OpenTrustRegion for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, OpenTrustRegion is a library that provides alternate
  orbital optimization.

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) OpenTrustRegion

* Upstream Dependencies |w---w| OpenTrustRegion |dr| LAPACK

**CMake Variables**

* :makevar:`ENABLE_OpenTrustRegion` |w---w| CMake variable toggling whether |PSIfour| builds with OpenTrustRegion
* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For OTR, set to an installation directory containing ``include/opentrustregion.h``
* :makevar:`OpenTrustRegion_DIR` |w---w| CMake variable to specify where pre-built OpenTrustRegion can be found. Set to installation directory containing ``lib/cmake/OpenTrustRegion/OpenTrustRegionConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_OpenTrustRegion` |w---w| CMake variable to force internal build of OpenTrustRegion instead of detecting pre-built
* :makevar:`CMAKE_INSIST_FIND_PACKAGE_OpenTrustRegion` |w---w| CMake variable to force detecting pre-built OpenTrustRegion and not falling back on internal build

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake -DENABLE_OpenTrustRegion=ON

B. Build *without* OpenTrustRegion

  .. code-block:: bash

    >>> cmake

C. Link against pre-built

  .. code-block:: bash

    >>> cmake -DENABLE_OpenTrustRegion=ON -DCMAKE_PREFIX_PATH=/path/to/OpenTrustRegion/root

  .. code-block:: bash

    >>> cmake -DENABLE_OpenTrustRegion=ON -DOpenTrustRegion_DIR=/path/to/otr/cmakeconfigdir

D. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DENABLE_OpenTrustRegion=ON -DCMAKE_PREFIX_PATH=/path/to/unwanted/OpenTrustRegion/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_OpenTrustRegion=ON
