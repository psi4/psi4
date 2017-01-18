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

.. _`apdx:testSuite`:

Test Suite and Sample Inputs
============================

|PSIfour| is distributed with an extensive test suite, which can
be found in :source:`tests`. After building the source code, these
can automatically be run by running ``ctest`` in the compilation
directory. More info on ``ctest`` options can be found
:ref:`here <faq:subsettests>`. Sample input files
can be found in the :source:`samples` subdirectory of the top-level Psi
directory. The samples and a brief description are provided below.

Sample inputs accessible through :ref:`interfaced executables
<sec:interfacing>` are bulleted below.

.. toctree::

   autodoc_testsuite_cfour
   autodoc_testsuite_chemps2
   autodoc_testsuite_dftd3
   autodoc_testsuite_dkh
   autodoc_testsuite_libefp
   autodoc_testsuite_erd
   autodoc_testsuite_gdma
   autodoc_testsuite_mrcc
   autodoc_testsuite_pcmsolver

Sample inputs for |PSIfour| as distributed are below.

.. comment This toctree directive only here to suppress warning at build time.
   include line below is doing the work.

.. toctree::
   :hidden:

   autodoc_testsuite_corepsi4

.. include:: autodoc_testsuite_corepsi4.rst

