
.. _`apdx:testSuite`:

Test Suite and Sample Inputs
============================

|PSIfour| is distributed with an extensive test suite, which can
be found in :source:`tests`. After building the source code, these
can automatically be run by running ``ctest`` in the compilation
directory. More info on ``ctest`` options can be found on the
`Wiki <https://github.com/psi4/psi4/wiki/4_Testing>`_. Sample input files
can be found in the the :source:`samples` subdirectory of the top-level Psi
directory. The samples and a brief description are provided below.

Sample inputs accessible through :ref:`interfaced executables
<sec:interfacing>` are bulleted below.

.. toctree::

   autodoc_testsuite_dftd3
   autodoc_testsuite_mrcc
   autodoc_testsuite_cfour
   autodoc_testsuite_libefp
   autodoc_testsuite_pcmsolver
   autodoc_testsuite_dmrg

Sample inputs for |PSIfour| as distributed are below.

.. comment This toctree directive only here to suppress warning at build time.
   include line below is doing the work.

.. toctree::
   :hidden:

   autodoc_testsuite_corepsi4

.. include:: autodoc_testsuite_corepsi4.rst

