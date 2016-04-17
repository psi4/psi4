
.. include:: autodoc_abbr_options_c.rst

.. index::
   single: harmonic vibrational analysis, vibrational analysis, thermochemical analysis

.. _`sec:thermo`:

Vibrational and Thermochemical Analysis
=======================================

.. codeauthor:: Rollin A. King
.. comment.. sectionauthor:: Rollin A. King and Lori A. Burns

*Module:* :ref:`Keywords <apdx:thermo`, :ref:`PSI Variables <apdx:thermo_psivar>`, :source:`THERMO <src/bin/thermo>`

Keywords
^^^^^^^^

.. include:: autodir_options_c/thermo__t.rst
.. include:: autodir_options_c/thermo__p.rst
.. include:: autodir_options_c/thermo__rotational_symmetry_number.rst

.. _`sec:thermoExamples`:

Examples
^^^^^^^^


.. index::
   pair: vibrational analysis; output

Output
^^^^^^

Text, text everywhere, and not a character on point. ::

   --------------------------------------------------------------------------------------------- ~
    Step     Total Energy     Delta E     MAX Force     RMS Force      MAX Disp      RMS Disp    ~
   --------------------------------------------------------------------------------------------- ~
     Convergence Criteria    1.00e-06 *    3.00e-04 *             o    1.20e-03 *             o  ~
   --------------------------------------------------------------------------------------------- ~
       1     -38.91591820   -3.89e+01      6.91e-02      5.72e-02 o    1.42e-01      1.19e-01 o  ~
       2     -38.92529543   -9.38e-03      6.21e-03      3.91e-03 o    2.00e-02      1.18e-02 o  ~
       3     -38.92540669   -1.11e-04      4.04e-03      2.46e-03 o    3.63e-02      2.12e-02 o  ~
       4     -38.92548668   -8.00e-05      2.30e-04 *    1.92e-04 o    1.99e-03      1.17e-03 o  ~
       5     -38.92548698   -2.98e-07 *    3.95e-05 *    3.35e-05 o    1.37e-04 *    1.05e-04 o  ~

The full list of keywords for thermo is provided in Appendix :ref:`apdx:thermo`.

Information on the Psithon function that drives frequency analyses is provided
at :py:func:`~driver.frequency`.

