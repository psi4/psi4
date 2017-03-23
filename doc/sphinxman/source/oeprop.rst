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

.. include:: autodoc_abbr_options_c.rst

.. index::
   single: OEProp
   pair: OEProp; theory


.. _`sec:oeprop`:

Evaluation of One-Electron Properties |w---w| :py:func:`~psi4.oeprop`
=====================================================================

.. codeauthor:: Robert M. Parrish and Andrew C. Simmonett
.. sectionauthor:: Andrew C. Simmonett

.. autofunction:: psi4.oeprop(wfn, \*args[, title])

|PSIfour| is capable of computing a number of one-electron properties
summarized in the table below. 

.. _`table:oe_features`:

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
   | Electrostatic potential, at nuclei | ESP_AT_NUCLEI         | Sets global variables "ESP AT CENTER n", n = 1 to natoms                          |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Electrostatic potential, on grid   | GRID_ESP              | Generates V at each point in grid_esp.dat. See :ref:`sec:oeprop_grid`             |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Electric field, on grid            | GRID_FIELD            | Generates {Ex,Ey,Ez} at each point grid_field.dat. See :ref:`sec:oeprop_grid`     |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Molecular orbital extents          | MO_EXTENTS            |                                                                                   |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Mulliken atomic charges            | MULLIKEN_CHARGES      |                                                                                   |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | L\ |o_dots|\ wdin atomic charges   | LOWDIN_CHARGES        |                                                                                   |
   +------------------------------------+-----------------------+-----------------------------------------------------------------------------------+
   | Wiberg bond indices                | WIBERG_LOWDIN_INDICES | Uses (L\ |o_dots|\ wdin) symmetrically orthogonalized orbitals                    |
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

The :py:func:`~psi4.property` function provides limited functionality, but is a lot easier to
use for correlated methods. For capabilities of :py:func:`~psi4.property` see the
corresponding section of the manual.


Basic Keywords
^^^^^^^^^^^^^^

Multipole moments may be computed at any origin, which is controlled by the
global |globals__properties_origin| keyword.  The keyword takes an array with
the following possible values:

.. _`table:oe_origin`:

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


.. _`sec:oeprop_grid`:


Properties evaluated on a grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Certain properties may be evaluated a user-specified grid points.  The grid
points are completely arbitrary and are specified by providing a file called
grid.dat containing the x,y,z values separated with spaces for each point in order::

    x1 y1 z1
    x2 y2 z2
    ..........
    xn yn zn

The grid.dat file is completely free form; any number of spaces and/or newlines
between entries is permitted.  The units of the coordinates in grid.dat are the
same as those used to specify the molecule's geometry, and the output
quantities are always in atomic units.  The requested properties will be
written out in the same order as the grid point specification in grid.dat; see
the above table for the format and file name of the output.

The grid may be generated in the input file using standard Python loops.  By
capturing the wavefunction used to evaluate the one-electron properties, the
values at each grid point may be captured as Python arrays in the input file::

    E, wfn = prop('scf', properties=["GRID_ESP", "GRID_FIELD"], return_wfn=True)
    Vvals = wfn.oeprop.Vvals()
    Exvals = wfn.oeprop.Exvals()
    Eyvals = wfn.oeprop.Eyvals()
    Ezvals = wfn.oeprop.Ezvals()

In this example, the *Vvals* array contains the electrostatic potential at each
grid point, in the order that the grid was specified, while the *Exvals*,
*Eyvals* and *Ezvals* arrays contain the *x*, *y* and *z* components of the
electric field, respectively; all of these arrays can be iterated and
manipulated using standard Python syntax.  For a complete demonstration of this
utility, see the :srcsample:`props4` test case.
