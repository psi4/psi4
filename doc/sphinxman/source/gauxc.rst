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

.. index:: GauXC
.. _`sec:gauxc`:

Interface to GauXC
==================

.. codeauthor:: Jonathon Misiewicz
.. sectionauthor:: Jonathon Misiewicz1

|PSIfour| contains code to interface to the open-source GauXC numerical integration module developed by Microsoft, which is available on conda.

Installing GauXC
~~~~~~~~~~~~~~~~~~

For users who obtain |PSIfour| from a package manager such as conda, use the same package manager to install GauXC, and then re-install |PSIfour|. For users installing |PSIfour| themselves, ensure that GauXC is available in your installation environment.

To control compilation and linking of the optional GauXC dependency required for the sn-LinK algorithm, 
here are the list of compile-time options provided.
  
* :makevar:`ENABLE_gauxc`: Compile Psi4 with support for GauXC.

* :makevar:`gauxc_DIR`: Location of the external GauXC install to compile Psi4 with, if using an external GauXC instance.

* :makevar:`gauxc_ENABLE_GPU`: Enable GPU support for the Psi4-GauXC interface class. When building GauXC internally within Psi4, this keyword controls whether to enable GPU support on the internally-built GauXC instance. When using an external GauXC build, this keyword must match the install: ``ON`` requires a GPU-enabled GauXC and ``OFF`` a CPU-only one, and configuration fails on a mismatch. With ``ON``, Psi4 also locates the CUDA Toolkit that GauXC links against (no CUDA compiler is needed); hint it with :makevar:`CUDAToolkit_ROOT` if necessary. Whether the compiled Psi4 has GPU-capable GauXC can be checked at runtime with ``psi4.addons("gauxc_gpu")``.

Seminumerical Linear Exchange
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GauXC can be used to compute exchange terms in an exact exchange self-consistent field computation. Literature references, keywords, and an overview of the algorithm may be found in :ref:`sec:scfsnlink`.

Kohn-Sham DFT
~~~~~~~~~~~~~

GauXC may be used as a replacement for |PSIfour|'s own KS-DFT engine. At present, we support outsourcing to GauXC for energies and analytic gradients and electronic hessians, as needed in stability analysis and TDDFT. GauXC does not provide support for analytic hessians. GauXC may be enabled by setting |scf__dft_xc_integrator_packages|, which expects an array of the available integrators, to be tried in order. GauXC may be sent by including "gauxc" in that list. See :ref:`notes <sec:ks-integrators>` on how to enable and disable other integrators.

The :ref:`keywords <apdx:scfgauxcdft>` are provided here. The most common keywords that will need editing are the number of radial integration grid points |scf__GAUXC_RADIAL_POINTS|, the number of angular integration grid points |scf__GAUXC_SPHERICAL_POINTS|, and whether or not to use GPUs |scf__GAUXC_USE_GPU|.

Programs using GauXC to compute DFT on KS-CPUs should cite [Petrone:2018:169]_, and programs using GauXC to compute KS-DFT shoudl cite [Williams-Young:2020:581058]_ (CUDA and distributed memory) and [Williams-Young:2021:102829]_ (HIP/SYCL). All programs should cite [Lehtola:2018:1]_.
