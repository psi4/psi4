
.. index::
   single: geometry optimization
   see: optimization; geometry optimization

.. _`sec:optking`:

Geometry Optimization
=====================

.. codeauthor:: Rollin A. King
.. sectionauthor:: Rollin A. King, Lori A. Burns

|PSIfour| carries out molecular optimizations using a module called
optking.  The optking program takes as input nuclear gradients and,
optionally, nuclear second derivatives |w---w| both in Cartesian coordinates.
The default minimization algorithm employs an empirical model Hessian,
redundant internal coordinates, a RFO step, and the BFGS Hessian update.

The principal literature references include the introduction of redundant
internal coordinates by Peng et al. [Peng:1996:49]_.
The general approach employed in this code
is similar to the "model Hessian plus RF method" described and tested by Bakken and
Helgaker [Bakken:2002:9160]_. (However, for separated
fragments, we have chosen not to employ by default their "extra-redundant"
coordinates defined by their "auxiliary interfragment" bonds.  These can be
included via the option :term:`ADD_AUXILIARY_BONDS`).

The internal coordinates are generated automatically based on an assumed bond
connectivity.  The connectivity is determined by testing if the interatomic
distance is less than the sum of atomic radii times the value of
:term:`COVALENT_CONNECT`.  If the user finds that some
connectivity is lacking by default, then this value may be increased.
Otherwise, the internal coordinate definitions may be modified.  If one
desires to see or modify the internal coordinates being used, then one can set
:term:`INTCOS_GENERATE_EXIT` to true.  The internal coordinate
definitions are provided in the file named "intco.dat".  See the :ref:`sec:optkingExamples`
section for more detail.

The ongoing development of optking is providing for unique treatment of
coordinates which connect distinct molecular fragments.  Thus, several keywords
relate to "interfragment modes", though many of these capabilities are
still under development.  Presently by default, separate fragments are bonded by
nearest atoms, and the whole system is treated as if it were part of one
molecule.  However, with the option :term:`FRAG_MODE`, fragments
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
:term:`INTRAFRAG_HESS` keyword.

All the common Hessian update schemes are available.  For formulas, see
Schlegel [Schlegel:1987:AIMQC]_ and Bofill [Bofill:1994:1]_.

The Hessian may be computed during an optimization using the following
:term:`FULL_HESS_EVERY` keyword

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

* Generate the internal coordinates and then stop::

   set intcos_generate_exit true
   optimize('scf')

  The coordinates may then be found in the file "intco.dat".  In this case, the file contains::
  
     F 1 3
     R      1     2
     R      1     3
     B      2     1     3
  
  The first line indicates a fragment containing atoms 1-3.  The following lines define
  two distance coordinates (bonds) and one bend coordinate.  This file can be modified, and if present,
  is used in subsequent optimizations.  Since the multiple-fragment coordinates are still under
  development, they are not documented here.  However, if desired, one can change the value
  of :term:`FRAG_MODE`, generate the internal coordinates, and see how multiple
  fragment systems are defined.
  
  Coordinates may be frozen or fixed by adding an asterisk after the letter of the coordinate.
  To optimize with the bond lengths fixed at their initial values, it is currently necessary to
  generate and then modify the internal coordinate definitions as follows::
  
     F 1 3
     R*     1     2
     R*     1     3
     B      2     1     3

.. index:: geometry optimization; convergence criteria

Convergence Criteria
^^^^^^^^^^^^^^^^^^^^

Optking monitors five quantities to evaluate the progress of a geometry 
optimization. These are (with their keywords) the change in energy 
:term:`MAX_ENERGY_G_CONVERGENCE`), the maximum element of 
the gradient (:term:`MAX_FORCE_G_CONVERGENCE`), the root-mean-square 
of the gradient (:term:`RMS_FORCE_G_CONVERGENCE`), the maximum element
of displacement (:term:`MAX_DISP_G_CONVERGENCE`), and the 
root-mean-square of displacement (:term:`RMS_DISP_G_CONVERGENCE`), 
all in internal coordinates and atomic units. Usually, these options will not 
be set directly. Primary control for geometry convergence lies with the keyword 
:term:`G_CONVERGENCE` which sets the aforementioned in accordance 
with Table :ref:`Geometry Convergence <table:optkingconv>`.

|
|

.. _`table:optkingconv`:

.. table:: Summary of sets of geometry optimization criteria available in |PSIfour|

    +-----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+----------------------------+
    | :term:`G_CONVERGENCE`       | Max Energy                 | Max Force                  | RMS Force                  | Max Disp                   | RMS Disp                   |
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
criterium and overwrites/appends it to the criteria set by :term:`G_CONVERGENCE`.
Note that this revokes the special convergence arrangements detailed in notes [#fe]_ and [#ff]_ 
and instead requires all active criteria to be fulfilled to 
achieve convergence. To avoid this revokation, turn on keyword :term:`FLEXIBLE_G_CONVERGENCE`.

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
at :py:func:`~driver.optimize`.

