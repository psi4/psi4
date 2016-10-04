.. include:: autodoc_abbr_options_c.rst

.. index::
   :pair: plugin; v2rdm_casscf

.. _`sec:v2rdm_casscf`:

Plugin v2rdm_casscf by A. E. DePrince
=====================================

.. codeauthor:: A. E. DePrince
.. sectionauthor:: Lori A. Burns

.. *Module:* :ref:`Keywords <apdx:pcm>`, :ref:`PSI Variables <apdx:pcm_psivar>`, :source:`PCMSolver <src/lib/libpsipcm>`

.. image:: https://img.shields.io/badge/home-v2rdm_casscf-5077AB.svg
   :target: https://github.com/edeprince3/v2rdm_casscf

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: https://github.com/edeprince3/v2rdm_casscf/blob/master/README.md

A variational 2-RDM-driven CASSCF plugin to |PSIfour|

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/v2rdm_casscf/badges/version.svg
     :target: https://anaconda.org/psi4/v2rdm_casscf

* v2rdm_casscf is available as a conda package for Linux and macOS.

* If using the |PSIfour| binary, v2rdm_casscf has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  v2rdm_casscf can be obtained through ``conda install v2rdm_casscf``.
  Then, hint its location with :envvar:`PYTHONPATH`.

* To remove a conda installation, ``conda remove v2rdm_casscf``.

**Source**

* .. image:: https://img.shields.io/github/tag/edeprince3/v2rdm_casscf.svg?maxAge=2592000
     :target: https://github.com/edeprince3/v2rdm_casscf

* If using |PSIfour| built from source and you want v2rdm_casscf built from
  from source also,
  built it, then hint its location with :envvar:`PYTHONPATH`.

