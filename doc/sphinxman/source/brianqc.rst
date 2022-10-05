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

.. index:: BrianQC
.. _`sec:brianqc`:

Interface to the BrianQC GPU module by the BrianQC team
=======================================================

.. codeauthor:: Gergely Kis
.. sectionauthor:: Gergely Kis

|PSIfour| contains code to interface to the BrianQC GPU module developed
by the `BrianQC team <https://www.brianqc.com/team>`_, which is available after a license agreement from
`https://brianqc.com/ <https://brianqc.com/>`_.

Installing BrianQC
~~~~~~~~~~~~~~~~~~

Please contact BrianQC at `https://brianqc.com/ <https://brianqc.com/>`_
to download the BrianQC GPU module and obtain a license.

Note that there are several prerequisites for using BrianQC, including
having a supported GPU available in the computing node and having the
proper GPU drivers installed. Please refer to the `BrianQC manual <https://www.brianqc.com/download/>`_
for a full list of prerequisites.

When installing BrianQC, choose the SDK installation by setting the
:envvar:`BRIANQC_SDK_INSTALL` envoronment variable to `1`.

Building BrianQC's user-built components
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After the installation, build BrianQC's example and sample programs,
which also creates the wrapper library and CMake configuration files
required to build |PSIfour| with BrianQC.

1. Create a build directory to keep the source tree clean.

.. code-block:: bash

    cd <brianqc_install_path>
    mkdir build
    cd build

2. Configure project and generate makefiles with CMake.
   You will require `Eigen <http://eigen.tuxfamily.org>`_ (tested with version 3.1.2)
   and `boost <https://www.boost.org/>`_ (tested with version 1.62).

.. code-block:: bash

    cmake ..

3. Build the examples and samples.

.. code-block:: bash

    make

4. Test the installation by starting a small calculation.
   Make sure to set the :envvar:`BRIANQC_INSTALL_PATH` environment variable to `<brianqc_install_path>`!

.. code-block:: bash

    export BRIANQC_INSTALL_PATH=<brianqc_install_path>
    bin/sample_hf_and_dft --molecule ../share/qc_molecules/cis-decalin.raw --basis ../share/basis_sets/cc-pvdz

Building |PSIfour| with BrianQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When CMake-ing |PSIfour|, set the :makevar:`ENABLE_BrianQC` CMake variable to `1`
and set the :makevar:`BrianQC_DIR` CMake variable to the path where BrianQC's
components have been built (usually `<brianqc_install_path>/build`), then build |PSIfour| normally.

Using BrianQC from |PSIfour|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use BrianQC, the :envvar:`BRIANQC_INSTALL_PATH` environment variable
must be set to the full path of your BrianQC installation.

There are two ways to enable BrianQC for a specific calculation:

* setting the |globals__brianqc_enable| input option to `True`;

* setting the :envvar:`BRIANQC_ENABLE` environment variable to `1`.

If both the input option and the environment variable are specified,
then the environment variable takes precedence.

When to enable BrianQC
""""""""""""""""""""""

The operation of BrianQC is transparent from a user's point of view; every
internal computation is either performed by |PSIfour| (if unsupported by BrianQC),
or taken over by BrianQC, yielding the same result to within the required
precision. In particular, BrianQC will respect all of |PSIfour|'s usual precision
parameters, like |scf__ints_tolerance| and |scf__e_convergence|.

BrianQC can speed up a number of internal computations, including Fock and
gradient computation. Thus, BrianQC will speed up any calculation involving
those terms, such as

* HF and DFT single point energies

* HF and DFT geometry optimizations

* HF and DFT frequency analysis

Note that not every term of every calculation can be handled by BrianQC, thus,
the actual speedup depends on the specifics of the calculation.

Necessary and recommended option values
"""""""""""""""""""""""""""""""""""""""

To ensure that a calculation works with BrianQC, the following options need to
be set when enabling BrianQC:

* BrianQC currently only handles the C1 molecular symmetry point group. Thus, if
  the molecule has any symmetries (which |PSIfour| would detect by default), the
  input geometry must contain the line ``symmetry c1`` to force |PSIfour| to
  disregard the symmetry.

To achieve peak performance, the following settings are recommended when enabling
BrianQC:

* By default, |PSIfour| uses density-fitted preiterations for SCF, which BrianQC
  doesn't handle yet. To ensure that no non-BrianQC-accelerated iterations are
  performed, disable the preiterations by setting |scf__df_scf_guess| to `False`.
* By default, |PSIfour| uses a disk-based Fock building, but BrianQC currently only
  accelerates direct Fock builds. To ensure that Fock building is accelerated by
  BrianQC, set |globals__scf_type| to `direct`.

Keywords
~~~~~~~~

.. include:: autodir_options_c/globals__brianqc_enable.rst

Environment Variables
~~~~~~~~~~~~~~~~~~~~~

.. envvar:: BRIANQC_SDK_INSTALL

   Toggle for type of BrianQC installation.
   See https://s3.amazonaws.com/streamnovation-brianqc-manual/BrianQC_Manual_v1_1_0.pdf for details.

.. envvar:: BRIANQC_INSTALL_PATH

   Location of BrianQC installation.
   See https://s3.amazonaws.com/streamnovation-brianqc-manual/BrianQC_Manual_v1_1_0.pdf for details.

.. envvar:: BRIANQC_ENABLE

   Runtime toggle of whether to use BrianQC.

