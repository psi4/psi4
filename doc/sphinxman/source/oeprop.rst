
.. include:: autodoc_abbr_options_c.rst

.. index::
   single: OEProp
   pair: OEProp; theory


.. _`sec:oeprop`:

Evaluation of One-Electron Properties
=====================================

.. codeauthor:: Robert M. Parrish and Andrew C. Simmonett
.. sectionauthor:: Andrew C. Simmonett

|PSIfour| is capable of computing a number of one-electron properties
summarized in the table below. 

.. table:: Current one-electron property capabilities of |PSIfour|

   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Feature                            | Keyword               | Notes                                                                             |
   +====================================+=======================+===================================================================================+
   | Electric dipole moment             | DIPOLE                |                                                                                   |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Electric quadrupole moment         | QUADRUPOLE            | Raw (traced) moments and traceless multipoles                                     |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | All moments up order N             | MULTIPOLE(N)          | Only raw (traced) moments. Sets global variables e.g. "DIPOLE X", "32-POLE XYYZZ" |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Transition dipole moment           | TRANSITION_DIPOLE     |                                                                                   |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Transition quadrupole moment       | TRANSITION_QUADRUPOLE |                                                                                   |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Electrostatic potential, at nuclei | ESP_AT_NUCLEI         | Sets global variables "ESP AT CENTER n", n = 1 to natoms                          |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Molecular orbital extents          | MO_EXTENTS            |                                                                                   |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Mulliken atomic charges            | MULLIKEN_CHARGES      |                                                                                   |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Löwdin atomic charges              | LOWDIN_CHARGES        |                                                                                   |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Wiberg bond indices                | WIBERG_LOWDIN_INDICES | Uses (Löwdin) symmetrically orthogonalized orbitals                               |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Mayer bond indices                 | MAYER_INDICES         |                                                                                   |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Natural orbital occupations        | NO_OCCUPATIONS        |                                                                                   |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+

There are two ways the computation of one-electron properties can be requested. 
Firstly, the properties can be evaluated from the last
computed one-particle density, using the following syntax::

  oeprop("MO_EXTENTS", "MULTIPOLE(4)", title = "hello!")

Note that it is the user's responsibility to ensure that the relaxed density
matrix is computed using the method of interest, which may require setting
additional keywords (see the method's manual section for details). The named
argument, *title*, is completely optional and is prepended to any
globals variables set during the computation.  The unnamed arguments are the
properties to be computed.  These can appear in any order, and multiple
properties may be requested, as in the example above.  Note that, due to Python
syntax restrictions, the title argument must appear after the list of
properties to compute.  The available properties are shown in the table above.

The syntax above works well for computing properties using the SCF
wavefunction, however, may be difficult (or impossible) to use for some of the
correlated levels of theory. Alternatively, one-electron properties can be
computed using the built-in property() function, e.g.::

  property('ccsd', properties=['dipole'])

The :py:func:`~driver.property` function provides limited functionality, but is a lot easier to
use for correlated methods. For capabilities of :py:func:`~driver.property` see the
corresponding section of the manual.


Basic Keywords
^^^^^^^^^^^^^^

Multipole moments may be computed at any origin, which is controlled by the
global |globals__properties_origin| keyword.  The keyword takes an array with
the following possible values:

.. table:: Allowed origin specifications

   +-------------------------------+-------------------------------------------------------------------------------+
   | Keyword                       | Interpretation                                                                |
   +===============================+===============================================================================+
   | [x, y, z]                     | Origin is at the coordinates, in the same units as the geometry specification |
   +-------------------------------+-------------------------------------------------------------------------------+
   | ["COM"]                       | Origin is at the center of mass                                               |
   +-------------------------------+-------------------------------------------------------------------------------+
   | ["NUCLEAR_CHARGE"]            | Origin is at the center of nuclear charge                                     |
   +-------------------------------+-------------------------------------------------------------------------------+
