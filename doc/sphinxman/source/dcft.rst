
.. include:: autodoc_abbr_options_c.rst

.. index::
   single: DCFT
   pair: DCFT; theory

.. _`sec:dcft`:

DCFT: Density Cumulant Functional Theory
========================================

.. codeauthor:: Alexander Yu. Sokolov and Andrew C. Simmonett
.. sectionauthor:: Alexander Yu. Sokolov

*Module:* :ref:`Keywords <apdx:dcft>`, :ref:`PSI Variables <apdx:dcft_psivar>`, :source:`DCFT <src/bin/dcft>`

Introduction
~~~~~~~~~~~~

Density cumulant functional theory (DCFT) is a density-based *ab initio* theory
that can compute electronic energies without the use of the wavefunction. The
theory starts by writing the exact energy expression in terms of the one- and
two-particle density matrices (OPDM and TPDM):

.. math:: 

    E_{\mathrm{DCFT}}  
    &= h_p^q \gamma_q^p + \frac{1}{2} g_{pq}^{rs} \gamma_{rs}^{pq}

Here we used Einstein convention for the summation over the repeated indices,
`h_p^q` and `g_{pq}^{rs}` are the standard one- and two-electron integrals,
`\gamma_p^q` and `\gamma_{pq}^{rs}` are the elements of OPDM and TPDM,
respectively. Na\"\ively, one might expect that it is possible to minimize the
energy functional in the equation above and obtain the exact energy. This is,
however, not trivial, as the density matrix elements `\gamma_p^q` and
`\gamma_{pq}^{rs}` cannot be varied arbitrarily, but must satisfy some
conditions that make sure that the density matrices are N-representable, *i.e.*
correspond to an antisymmetric N-electron wavefunction. Unfortunately, no
simple set of necessary and sufficient N-representability conditions are known,
and some of the known conditions are not easily imposed. In addition, the lack
of separability of the density matrices may result in the loss of
size-consistency and size-extensivity.


Theory
~~~~~~

Hello, World!

Minimal Input
~~~~~~~~~~~~~

Hello, World!

Recommendations for Convergence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

