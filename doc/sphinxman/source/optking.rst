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
   single: geometry optimization, optimization

.. _`sec:optking`:

Geometry Optimization
=====================

.. codeauthor:: Rollin A. King
.. sectionauthor:: Rollin A. King and Lori A. Burns

*Module:* :ref:`Keywords <apdx:optking>`, :ref:`PSI Variables <apdx:optking_psivar>`, :source:`OPTKING <psi4/src/psi4/optking>`

|PSIfour| carries out molecular optimizations using a module called
optking.  The optking program takes as input nuclear gradients and,
optionally, nuclear second derivatives |w---w| both in Cartesian coordinates.
The default minimization algorithm employs an empirical model Hessian,
redundant internal coordinates, an RFO step, and the BFGS Hessian update.

The principal literature references include the introduction of redundant
internal coordinates by Peng et al. [Peng:1996:49]_.
The general approach employed in this code
is similar to the "model Hessian plus RF method" described and tested by Bakken and
Helgaker [Bakken:2002:9160]_. (However, for separated
fragments, we have chosen not to employ by default their "extra-redundant"
coordinates defined by their "auxiliary interfragment" bonds.  These can be
included via the option |optking__add_auxiliary_bonds|).

The internal coordinates are generated automatically based on an assumed bond
connectivity.  The connectivity is determined by testing if the interatomic
distance is less than the sum of atomic radii times the value of
|optking__covalent_connect|. If the user finds that some
connectivity is lacking by default, then this value may be increased.
Otherwise, the internal coordinate definitions may be modified directly.  If one
desires to see or modify the internal coordinates being used, then one can set
|optking__intcos_generate_exit| to true.  The internal coordinate
definitions are provided in the file with extension ".intco".  See the :ref:`sec:optkingExamples`
section for more detail.

.. warning:: The selection of a Z-matrix input, and in particular the inclusion
   of dummy atoms, has no effect on the behavior of the optimizer, which begins
   from a Cartesian representation of the system.

The ongoing development of optking is providing for unique treatment of
coordinates which connect distinct molecular fragments.  Thus, several keywords
relate to "interfragment modes", though many of these capabilities are
still under development.  Presently by default, separate fragments are bonded by
nearest atoms, and the whole system is treated as if it were part of one
molecule.  However, with the option |optking__frag_mode|, fragments
may instead be related by a unique set of interfragment coordinates defined by
reference points within each fragment.  The reference points can be atomic
positions (current default), linear combinations of
atomic positions, or located on the principal axes (not yet working).

Basic Keywords
^^^^^^^^^^^^^^

.. include:: autodir_options_c/optking__opt_type.rst
.. include:: autodir_options_c/optking__step_type.rst
.. include:: autodir_options_c/optking__geom_maxiter.rst
.. include:: autodir_options_c/optking__g_convergence.rst
.. include:: autodir_options_c/optking__full_hess_every.rst
.. include:: autodir_options_c/optking__intcos_generate_exit.rst

.. index:: geometry optimization; minima
.. _`sec:optkingExamples`:

Optimizing Minima
^^^^^^^^^^^^^^^^^

First, define the molecule and basis in the input. ::

   molecule h2o {
     O
     H 1 1.0
     H 1 1.0 2 105.0
   }
    
   set basis dz

Then the following are examples of various types of calculations that can be completed.

* Optimize a geometry using default methods (RFO step)::

   optimize('scf')

* Optimize using Newton-Raphson steps instead of RFO steps::

   set step_type nr
   optimize('scf')

* Optimize using energy points instead of gradients::

   optimize('scf', dertype='energy')

* Optimize while limiting the initial step size to 0.1 au::

   set intrafrag_step_limit 0.1
   optimize('scf')

* Optimize while always limiting the step size to 0.1 au::

   set {
     intrafrag_step_limit     0.1
     intrafrag_step_limit_min 0.1
     intrafrag_step_limit_max 0.1
   }

   optimize('scf')

* Optimize while calculating the Hessian at every step::

   set full_hess_every 1
   optimize('scf')

Hessian
^^^^^^^

If Cartesian second derivatives are available, optking can read them
and transform them into internal coordinates to make an initial Hessian in
internal coordinates.  Otherwise, several empirical Hessians are available,
including those of Schlegel [Schlegel:1984:333]_ and Fischer and Almlof
[Fischer:1992:9770]_.
Either of these or a simple diagonal Hessian may be selected using the 
|optking__intrafrag_hess| keyword.

All the common Hessian update schemes are available.  For formulas, see
Schlegel [Schlegel:1987:AIMQC]_ and Bofill [Bofill:1994:1]_.

The Hessian may be computed during an optimization using the 
|optking__full_hess_every| keyword.

.. index:: 
   pair: geometry optimization; transition state
   pair: geometry optimization; IRC
   single: geometry optimization; constrained

Transition States, Reaction Paths, and Constrained Optimizations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Calculate a starting Hessian and optimize the "transition state" of
  linear water (note that without a reasonable starting geometry and
  Hessian, such a straightforward search often fails)::

   molecule h2o {
      O
      H 1 1.0
      H 1 1.0 2 160.0
   }

   set {
      basis dz
      full_hess_every 0
      opt_type ts
   }

   optimize('scf')

* At a transition state (planar HOOH), compute the second derivative, and
  then follow the intrinsic reaction path to the minimum::

   molecule hooh {
      symmetry c1
      H
      O 1 0.946347
      O 2 1.397780 1  107.243777
      H 3 0.946347 2  107.243777   1 0.0
   }

   set {
      basis dzp
      opt_type irc
      geom_maxiter 50
   }

   frequencies('scf')
   optimize('scf')

* Optimize a geometry (HOOH) at a frozen dihedral angle of 90 degrees. ::

   molecule {
     H
     O 1 0.90
     O 2 1.40 1 100.0
     H 3 0.90 2 100.0 1 90.0
   }

   set optking {
     frozen_dihedral = ("
       1 2 3 4
     ")
   }
   optimize('scf')

* To instead freeze the two O-H bond distances ::

   set optking {
     frozen_distance = ("
       1  2
       3  4
     ")
   }

For bends, the corresponding keyword is "frozen_bend".

* To freeze the cartesian coordinates of atom 2 ::

   freeze_list = """
     2 xyz
   """
   set optking frozen_cartesian $freeze_list

* To freeze only the y coordinates of atoms 2 and 3 ::

   freeze_list = """
     2 y
     3 y
   """
   set optking frozen_cartesian $freeze_list

* To optimize toward a value of 0.95 Angstroms for the distance between 
  atoms 1 and 3, as well as that between 2 and 4 ::

   set optking {
     fixed_distance = ("
       1  3 0.95
       2  4 0.95
     ")
   }

Note that the effect of the frozen and fixed keywords is independent of
how the geometry of the molecule was input (whether Z-matrix or cartesian, etc.)..

* To scan the potential energy surface by optimizing at several fixed values
  of the dihedral angle of HOOH. ::

   molecule hooh {
    0 1
    H  0.850718   0.772960    0.563468
    O  0.120432   0.684669   -0.035503
    O -0.120432  -0.684669   -0.035503
    H -0.850718  -0.772960    0.563468
   }
   
   set {
     basis cc-pvdz
     intrafrag_step_limit 0.1
   }
   
   dihedrals = [100,110,120,130,140,150]
   PES = []
   
   for phi in dihedrals:
     my_string = "1 2 3 4 " + str(phi)
     set optking fixed_dihedral = $my_string
     E = optimize('scf')
     PES.append((phi, E))
   
   print "\n\tcc-pVDZ SCF energy as a function of phi\n"
   for point in PES:
     print "\t%5.1f%20.10f" % (point[0], point[1])


Dealing with problematic optimizations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Although optking is continuously improved with robustness in mind, some
attempted optimizations will inevitably fail to converge to the desired minima.
For difficult cases, the following suggestions are made.

* As for any optimizer, computing the Hessian and limiting the step size will
  successfully converge a higher percentage of cases.  The default settings have
  been chosen because they perform efficiently for common, representative test sets.
  More restrictive, cautious steps are sometimes necessary.

* |optking__dynamic_level| allows optking to change the method of optimization
  toward algorithms that, while often less efficient, may help to converge difficult
  cases.  If this is initially set to 1, then optking, as poor steps are detected,
  will increase the level through several forms of more robust and cautious algorithms.
  The changes will reduce the trust radius, allow backward steps (partial line
  searching), add cartesian coordinates, switch to cartesian coordinates, and take
  steepest-descent steps.

* The developers have found the |optking__opt_coordinates| set to "BOTH" which
  includes both the redundant internal coordinate set, as well as cartesian coordinates,
  works well for systems with long 'arms' or floppy portions of a molecule poorly
  described by local internals.

Direct manipulation of the optmization coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* Generate the internal coordinates and then stop::

   set intcos_generate_exit true
   optimize('scf')

  The coordinates may then be found in the "intco" file.  In this case, the file contains::
  
     F 1 3
     R      1     2
     R      1     3
     B      2     1     3
     C      1
            1    1.000000
     C      1
            2    1.000000
     C      1
            3    1.000000

  The first line indicates a fragment containing atoms 1-3.  The following lines define
  two distance coordinates (bonds) and one bend coordinate.  This file can be modified, and if present,
  is used in subsequent optimizations.  The lines below the simple internal coordinates
  specify linear combinations of coordinates.  In the simplest default case, the lines
  above simply define combination coordinates which are identical to the simple internals.
  If |optking__opt_coordinates| specifies delocalized coordinates, then the combinations
  will be more complex.
 
  Since the multiple-fragment coordinates are still under
  development, they are not documented here.  However, if desired, one can change the value
  of |optking__frag_mode|, generate the internal coordinates, and see how multiple
  fragment systems are defined.
  
  Coordinates may be frozen by adding an asterisk after the letter of the coordinate.  The
  asterisk results in that internal coordinate being frozen at its initial value.  The
  "intco" file below for water specifies an optimization with both O-H bonds frozen.::
  
     F 1 3
     R*     1     2
     R*     1     3
     B      2     1     3

  If one instead wishes to optimize toward ("fix") a value that is not satisfied by the
  initial structure, then the value is added to the end of the line.  The following
  corresponds to an optimization that will add additional forces to move the O-H bonds
  to 1.70 au. ::

     F 1 3
     R      1     2     1.70
     R      1     3     1.70
     B      2     1     3

.. index:: 
   pair: geometry optimization; convergence criteria

Convergence Criteria
^^^^^^^^^^^^^^^^^^^^

Optking monitors five quantities to evaluate the progress of a geometry 
optimization. These are (with their keywords) the change in energy 
(|optking__max_energy_g_convergence|), the maximum element of 
the gradient (|optking__max_force_g_convergence|), the root-mean-square 
of the gradient (|optking__rms_force_g_convergence|), the maximum element
of displacement (|optking__max_disp_g_convergence|), and the 
root-mean-square of displacement (|optking__rms_disp_g_convergence|), 
all in internal coordinates and atomic units. Usually, these options will not 
be set directly. Primary control for geometry convergence lies with the keyword 
|optking__g_convergence| which sets the aforementioned in accordance 
with Table :ref:`Geometry Convergence <table:optkingconv>`.

|
|

.. _`table:optkingconv`:

.. table:: Summary of sets of geometry optimization criteria available in |PSIfour|

    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | |optking__g_convergence|    | Max Energy                 | Max Force                  | RMS Force                  | Max Disp                   | RMS Disp                   |
    +=============================+============================+============================+============================+============================+============================+
    | NWCHEM_LOOSE [#fd]_         |                            | :math:`4.5 \times 10^{-3}` | :math:`3.0 \times 10^{-3}` | :math:`5.4 \times 10^{-3}` | :math:`3.6 \times 10^{-3}` |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | GAU_LOOSE [#ff]_            |                            | :math:`2.5 \times 10^{-3}` | :math:`1.7 \times 10^{-3}` | :math:`1.0 \times 10^{-2}` | :math:`6.7 \times 10^{-3}` |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | TURBOMOLE [#fd]_            | :math:`1.0 \times 10^{-6}` | :math:`1.0 \times 10^{-3}` | :math:`5.0 \times 10^{-4}` | :math:`1.0 \times 10^{-3}` | :math:`5.0 \times 10^{-4}` |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | GAU [#fc]_ [#ff]_           |                            | :math:`4.5 \times 10^{-4}` | :math:`3.0 \times 10^{-4}` | :math:`1.8 \times 10^{-3}` | :math:`1.2 \times 10^{-3}` |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | CFOUR [#fd]_                |                            |                            | :math:`1.0 \times 10^{-4}` |                            |                            |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | QCHEM [#fa]_ [#fe]_         | :math:`1.0 \times 10^{-6}` | :math:`3.0 \times 10^{-4}` |                            | :math:`1.2 \times 10^{-3}` |                            |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | MOLPRO [#fb]_ [#fe]_        | :math:`1.0 \times 10^{-6}` | :math:`3.0 \times 10^{-4}` |                            | :math:`3.0 \times 10^{-4}` |                            |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | GAU_TIGHT [#fc]_ [#ff]_     |                            | :math:`1.5 \times 10^{-5}` | :math:`1.0 \times 10^{-5}` | :math:`6.0 \times 10^{-5}` | :math:`4.0 \times 10^{-5}` |
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | GAU_VERYTIGHT [#ff]_        |                            | :math:`2.0 \times 10^{-6}` | :math:`1.0 \times 10^{-6}` | :math:`6.0 \times 10^{-6}` | :math:`4.0 \times 10^{-6}` | 
    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+

.. rubric:: Footnotes

.. [#fa] Default
.. [#fb] Baker convergence criteria are the same.
.. [#fc] Counterpart NWCHEM convergence criteria are the same.
.. [#fd] Convergence achieved when all active criteria are fulfilled.
.. [#fe] Convergence achieved when **Max Force** and one of **Max Energy** or **Max Disp** are fulfilled.
.. [#ff] Normal convergence achieved when all four criteria (**Max Force**, **RMS Force**,
         **Max Disp**, and **RMS Disp**) are fulfilled. To help with flat 
         potential surfaces, alternate convergence achieved when 100\ :math:`\times`\ *rms force* is less 
         than **RMS Force** criterion.

For ultimate control, specifying a value for any of the five monitored options activates that
criterium and overwrites/appends it to the criteria set by |optking__g_convergence|.
Note that this revokes the special convergence arrangements detailed in notes [#fe]_ and [#ff]_ 
and instead requires all active criteria to be fulfilled to 
achieve convergence. To avoid this revokation, turn on keyword |optking__flexible_g_convergence|.

.. index::
   pair: geometry optimization; output

Output
^^^^^^

The progress of a geometry optimization can be monitored by grepping the output file for the
tilde character (``~``). This produces a table like the one below that shows
for each iteration the value for each of the five quantities and whether the criterion
is active and fulfilled (``*``), active and unfulfilled ( ),  or inactive (``o``). ::

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

The full list of keywords for optking is provided in Appendix :ref:`apdx:optking`.

Information on the Psithon function that drives geometry optimizations is provided
at :py:func:`~psi4.optimize`.

