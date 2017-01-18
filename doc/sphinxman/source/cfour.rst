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

.. index:: Cfour
.. _`sec:cfour`:

Interface to CFOUR by J. Stanton & J. Gauss
===========================================

.. codeauthor:: Lori A. Burns
.. sectionauthor:: Lori A. Burns

*Module:* :ref:`Keywords <apdx:cfour>`, :ref:`PSI Variables <apdx:cfour_psivar>`, :ref:`Samples <apdx:testSuitecfour>`

|PSIfour| contains code to interface to the Cfour quantum chemistry suite of
John F. Stanton (U. Texas, Austin) and J\ |u_dots|\ rgen Gauss (U. Mainz),
which is available after a license agreement from 
`http://www.cfour.de/ <http://www.cfour.de/>`_.

Installation
~~~~~~~~~~~~

Follow the instructions provided with the Cfour download to install the
executable or to build the source. To by used by |PSIfour|, the program
binary (:program:`xcfour`) must be found in your :envvar:`PATH` or
:envvar:`PSIPATH`.  The ``GENBAS`` file containing basis sets in Cfour
format is not necessary for this interface, but if you prefer to access
basis sets the "Cfour way" using a custom ``GENBAS`` file (the distributed
one is included with the interface), it, too, must be in :envvar:`PATH` or
:envvar:`PSIPATH`. If |PSIfour| is unable to execute the binary, an error
will be reported.

.. .. caution:: The p4c4 interface isn't in the master branch nor will it be in
..    the near future. To run this code, (1) build the ``c4`` branch of psi4,
..    (2) find a copy of cfour and put it in :envvar:`PATH` or
..    :envvar:`PSIPATH`, and (3) clone https://github.com/loriab/qcdb.git
..    python module and prepend :envvar:`PYTHONPATH` with the top qcdb
..    directory (the path added to PYTHONPATH should have one "qcdb" in it;
..    the cloned qcdb is what needs to be imported in preference to the one
..    already in psi4). Execute psi4 as usual.

.. caution:: The p4c4 interface hasn't been fully adapted for the new March 2014 version.

Cfour for |PSIfour| Users
~~~~~~~~~~~~~~~~~~~~~~~~~

* Set memory as usual

* Set molecule as usual

* Set basis set as usual (Cfour only cares about orbital basis, no fitting
  bases)

* Set the task as usual, indicating Cfour as the intended code by
  prepending "c4-" to the method argument. So ``energy('scf')`` becomes
  ``energy('c4-scf')`` and ``optimize('ccsd(t)')`` becomes
  ``optimize('c4-ccsd(t)')``. Find available methods for
  :py:func:`~psi4.energy` at :ref:`Energy (CFOUR) <table:energy_cfour>`
  and for :py:func:`~psi4.optimize` at :ref:`Gradient (CFOUR)
  <table:grad_cfour>`.

* Generally, the p4c4 interface will handle best practices for path of
  execution: ``vcc``/``ecc``, derivative type, *etc.* The user is still
  responsible for setting convergence, frozen core, guess, diis, *etc.*
  For the moment, so-called "best-practices" keywords are summarized at
  :ref:`Best Practices <table:cfour_cc_program>`.

* For the type of computation intended, find appropriate options at
  :ref:`Keywords <apdx:cfour>`. These keyword summaries contain the same
  information as the `proper CFOUR options list
  <http://slater.chemie.uni-mainz.de/cfour/index.php?n=Main.ListOfKeywordsInAlphabeticalOrder>`_
  plus notes on keyword relevance when run through |PSIfour|.  Information
  at the `CFOUR manual
  <http://slater.chemie.uni-mainz.de/cfour/index.php?n=Main.Manual>`_ may
  also be useful, as may the many samples at :source:`samples/cfour`.

* Set Cfour keywords just like |PSIfour| keywords. The names of keywords
  are unchanged beyond a prepended "cfour\_". (Though be aware that common
  abbreviations like CALC and REF must be fully spelled out as
  |cfour__cfour_calc_level| and |cfour__cfour_reference| when used in
  |PSIfour|.)

* In limited trial cases, keywords nominally directed at non-Cfour modules
  are translated into their Cfour counterparts. For example, setting
  |scf__reference| will appropriately set |cfour__cfour_reference|. For a
  list of applicable keywords, see source of
  :py:func:`qcdb.cfour.muster_psi4options`.

* Consult :ref:`sec:cfourFunctionality` for information on what Cfour
  functionality is accessible through |PSIfour|.

|PSIfour| for Cfour Users
~~~~~~~~~~~~~~~~~~~~~~~~~

In the simplest use of the Psi4/Cfour interface, a |PSIfour| input file
can simply "wrap" a ``ZMAT`` file and execute :program:`xcfour`. This is
illustrated in the following example::

    cfour {
    UHF-SCF energy calculation 
    N
    H 1 R
    H 1 R 2 A
    
    R=1.008
    A=105.0
    
    *ACES2(CALC=HF,BASIS=qz2p
    MULT=2,REF=UHF
    OCCUPATION=3-1-1-0/3-0-1-0
    SCF_CONV=12
    MEMORY=20000000)
    }
    
    energy('cfour')

Here, the contents of the ``cfour {...}`` block are written directly to a
``ZMAT`` file. This is joined by a default ``GENBAS`` file
(:source:`share/basis/GENBAS`).  To preferentially use your own ``GENBAS``,
place it in :envvar:`PATH` or :envvar:`PSIPATH`. The line calling
:py:func:`~psi4.energy` with argument ``'cfour'`` invokes
:program:`xcfour`.

After execution of the ``energy('cfour')`` line completes, Cfour results
are read back into |PSIfour| format and are thereafter accessible for
further processing in the input file. See :ref:`sec:cfourOutput` for
details. This storage of results in variables and arrays *in memory* for
the duration of the |PSIfour| instance (as opposed to solely *in files*)
is the only advantage thus far incurred by the P4C4 interface. We'll call
this mode of basic utility the "sandwich" mode.

Molecule specification in |PSIfour| allows Cartesians, Z-matrices, mixed
Cartesian/Z-matrix, negation of variables, delayed specification of
variables, specification of fragments, etc., all in a whitespace-tolerant
format. See :ref:`sec:moleculeSpecification` for details and
:srcsample:`cfour/mints5` for examples. When a |PSIfour|-style molecule is
supplied, its geometry is written to ``ZMAT`` in Cartesian form and the
|cfour__cfour_coordinates|\ =CARTESIAN, |cfour__cfour_units|\ =ANGSTROM,
|cfour__cfour_charge|, and |cfour__cfour_multiplicity| keywords are set
appropriately in the ``*CFOUR(...)`` directive.

.. warning:: There exist molecules (*e.g.*, allene) where the
   inertial frame is not unique (planes along atoms or between
   atoms). The orientation reconciling machinery currently does not
   handle these cases and will fail with "Axis unreconcilable between
   QC programs". I will get to this soon.

Whenever the molecule is supplied in |PSIfour| format, the job control
keywords must be too. All :ref:`Cfour keywords <apdx:cfour>` are the usual
ones, prepended by ``cfour_`` to avoid any possible name conflicts.  As
detailed in :ref:`sec:jobControl`, setting keywords is flexible in
format. The previous example translates to::

    # UHF-SCF energy calculation 

    molecule {
    0 2                                          # multiplicity from the MULT keyword
    N
    H 1 R
    H 1 R 2 A
    
    R=1.008
    A=105.0
    }
    
    set {
    cfour_CALC_level=HF                          # only full keyword names allowed
    cfour_BASIS=qz2p
    #MULT=2                                      # now in molecule {...} block
    cfour_REFerence=UHF
    cfour_OCCUPATION [[3, 1, 1, 0], [3,0,1,0] ]  # arrays in python notation
    cfour_SCF_CONV=12
    cfour_MEMORY=20000000
    }
    
    energy('cfour')

Here, note that none of capitalization, equals sign, or whitespace matter
for the keyword commands. Specification of strings and integers requires no
translation; :ref:`booleans <op_c_boolean>` have extended freedom of
format; arrays must be translated into Python-style (square-bracket
bounded and comma delimited) of appropriate dimension. There are many
sample inputs in :source:`tests/cfour/` starting with ``sp-`` that take
examples from the Cfour manual and first run them in sandwich mode and
then run them as translated into |PSIfour| format.

.. note:: |PSIfour| only recognizes keywords by their full name, so the common
   Cfour keyword abbreviations CALC, REF, etc. must be replaced by their
   proper names of |cfour__cfour_calc_level|, |cfour__cfour_reference|, etc.

Whenever the molecule is supplied in |PSIfour| format, it is possible to
perform geometry optimizations where Cfour supplies the gradient and the
|PSIfour| module :ref:`optking <sec:optking>` drives the structural
changes. Because of the limitations on geometry specification for
optimizations in Cfour, optking-driven optimizations are the *only*
optimizations allowed in the P4C4 interface. (The exception is sandwich
mode, which, of course, permits optimizations with the Cfour optimizer.)
Below is an example of a geometry optimization::

    memory 200 mb
    
    molecule {
    O
    H 1 R
    H 1 R 2 A
    
    R=0.958
    A=104.5
    }
    
    set {
    
    cfour_CALC_level CCSD(T)
    cfour_BASIS      DZP
    cfour_CC_CONV    12
    cfour_LINEQ_CONV 12
    cfour_SCF_CONV   12
    g_convergence    cfour
    }

    optimize('cfour')

Note that the primary change is the exchange of :py:func:`~psi4.energy`
for :py:func:`~psi4.optimize` to trigger an optimization.  Setting
|optking__g_convergence|\ =CFOUR provides a good imitation of Cfour
default convergence criteria. Although Cfour produces gradients only in
its standard orientation and atom ordering, these are transformed back to
input orientation by the P4C4 interface. Several sample inputs in
:source:`tests/cfour/` starting with ``opt-`` show basic geometry
optimizations. :srcsample:`cfour/mints5-grad` shows optimizations from a
variety of molecule input formats, and :srcsample:`cfour/psi-ghost-grad`
shows an optimization with ghosted atoms. To obtain a single gradient
*sans* optimization, call instead :py:func:`~psi4.gradient`.

Note that it can be convenient to monitor the progress of a geometry
optimization by grepping the tilde ``~`` character. ::
 
   Measures of convergence in internal coordinates in au.
   Criteria marked as inactive (o), active & met (*), and active & unmet ( ).
   --------------------------------------------------------------------------------------------- ~
    Step     Total Energy     Delta E     MAX Force     RMS Force      MAX Disp      RMS Disp    ~
   --------------------------------------------------------------------------------------------- ~
     Convergence Criteria    1.00e-06 *    3.00e-04 *    1.00e-06 *    1.20e-03 *             o  ~
   --------------------------------------------------------------------------------------------- ~
       1     -76.33224285   -7.63e+01      2.41e-03      1.60e-03      1.51e-02      8.82e-03 o  ~
       2     -76.33226097   -1.81e-05      4.84e-04      4.03e-04      7.71e-04 *    7.04e-04 o  ~
       3     -76.33226140   -4.39e-07 *    4.31e-05 *    3.58e-05      9.89e-05 *    8.93e-05 o  ~
       4     -76.33226141   -4.26e-09 *    9.76e-07 *    6.58e-07 *    6.22e-06 *    3.71e-06 o  ~
   --------------------------------------------------------------------------------------------------------------- ~
    Step         Total Energy             Delta E       MAX Force       RMS Force        MAX Disp        RMS Disp  ~
   --------------------------------------------------------------------------------------------------------------- ~
       1     -76.332242848098    -76.332242848098      0.00241281      0.00160359      0.01507630      0.00881949  ~
       2     -76.332260965382     -0.000018117284      0.00048446      0.00040256      0.00077146      0.00070447  ~
       3     -76.332261404452     -0.000000439070      0.00004307      0.00003577      0.00009889      0.00008926  ~
       4     -76.332261408714     -0.000000004262      0.00000098      0.00000066      0.00000622      0.00000371  ~
   --------------------------------------------------------------------------------------------------------------- ~

The above example also shows the total memory for the computation being
set in |PSIfour| format. See :ref:`sec:memory` for details. When
specified, the memory value is passed on to Cfour by setting keywords
|cfour__cfour_memory_size| and |cfour__cfour_mem_unit|\ =MB.

|PSIfour| has an extensive :ref:`basis set library <apdx:basisElement>` in
Gaussian94 format. See :ref:`sec:basisSets` for details.  Contrasts to
Cfour basis handling include identifying basis sets by standard name
(aug-cc-pVDZ instead of AUG-PVDZ), direct handles for
diffuse-function-pruned sets (*e.g.*, jun-cc-pVDZ), case insensitivity,
appropriate setting of spherical/Cartesian depending on basis set design,
and syntax to set different basis sets to different classes of atoms
without listing each atom. All of these features are available to Cfour by
using the |mints__basis| keyword instead of |cfour__cfour_basis|
(accompanied, of course, by specifying the molecule |PSIfour|-style).
Internally, |PSIfour| processes the basis set as usual, then translates
the basis set format and writes out a ``GENBAS`` file with an entry for
each atom. The P4C4 interface sets keyword |cfour__cfour_basis|\ =SPECIAL
and |cfour__cfour_spherical| as appropriate, then writes the basis section
necessary for SPECIAL below the ``*CFOUR(...)`` block. (I'm sorry that the
name of the basis doesn't appear in ``ZMAT``, but the combination of the
~14 character basis name limit and the absence of a comment line marker
rather precludes that helpful label.)

The input below employs a |PSIfour| library basis set and also introduces
the final stage of conversion toward |PSIfour| format. Instead of the
generic ``'cfour'``, the computational method is specified as the first
argument to the :py:func:`~psi4.optimize` call. In the computational
command below, the string argument ``'c4-ccsd(t)'`` directs that a CCSD(T)
computation be run using Cfour (as opposed to ``'ccsd(t)'`` which would
call |PSIfour| CC code). Specifying computational method in this manner
sets |cfour__cfour_calc_level| from the argument and
|cfour__cfour_deriv_level| as appropriate from the function call:
:py:func:`~psi4.energy`, :py:func:`~psi4.gradient`, or
:py:func:`~psi4.optimize`. If those keywords are also set explicitly to
contradictory values, the interface will complain. ::

   memory 2 gb

   molecule CH2F2  {
     units au
     C     0.0000000000  -0.0000000000   1.0890958457
     F     0.0000000000  -2.1223155812  -0.4598161475
     F    -0.0000000000   2.1223155812  -0.4598161475
     H     1.7084139850   0.0000000000   2.1841068002
     H    -1.7084139850  -0.0000000000   2.1841068002
   }

   set basis aug-cc-pvdz
   set rms_force_g_convergence 6
   set cfour_abcdtype aobasis
   set cfour_scf_conv 12
   set cfour_cc_conv 12
   set cfour_lineq_conv 12

   optimize('c4-ccsd(t)')

The utility of this method specification is that examination can be made
of the reference, the derivative level, the excitation level, *etc.* and
some options can be set according to best practices. Practically speaking,
|cfour__cfour_cc_program| (and eventually |cfour__cfour_abcdtype|) will
always be set to the :ref:`fastest safe value <table:cfour_cc_program>`.
For example, the input above will run with |cfour__cfour_cc_program|\ =ECC
unless explicitly set to VCC.

An advantage of |PSIfours| Python driver is that any number of common
work-up procedures can be automated and wrapped around the
conventional single-point and optimization procedures at the heart of all
quantum chemistry codes. Three core "wrappers" available in |PSIfour| are
:py:func:`~driver_nbody.nbody_gufunc`,
:py:func:`~wrapper_database.database`, and
:py:func:`~driver_cbs.complete_basis_set`; read their respective sections
for details, but an overview is provided here. :py:func:`~driver_nbody.nbody_gufunc`
computes the interaction energy of a bimolecular complex (counterpoise-corrected,
not, or both). ::

   molecule dimer {
     Ne
   --
     Ne 1 R
     symmetry c1
   }
   
   Rvals=[2.5, 3.0, 4.0]
   set basis aug-cc-pVDZ
   
   for R in Rvals:
     dimer.R = R
     ecp = cp('c4-mp2')
     print_stdout('R [A] = %.1f  IE [kcal/mol] = %.3f\n' % (R, psi_hartree2kcalmol * ecp))
   
yields ::

   R [A] = 2.5  IE [kcal/mol] = 0.804
   R [A] = 3.0  IE [kcal/mol] = 0.030
   R [A] = 4.0  IE [kcal/mol] = -0.014

Next, the :py:func:`~wrapper_database.database` wrapper allows any computational
model chemistry to be applied a predefined collection of molecules. Thus
an input ::

   set {
       basis jun-cc-pvdz
       d_convergence 9
   }
   
   database('c4-mp2','nbc10',cp='on',subset='MeMe')

yields the counterpoise-corrected interaction energy for several points
along the dissociation curve of methane dimer, which is a member of the
:srcdb:`NBC10` database::

   //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
   //       Database nbc10 Results      //
   //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
 
   For each VARIABLE requested by tabulate, a 'Reaction Value' will be formed from
   'Reagent' values according to weightings 'Wt', as for the REQUESTED ENERGY below.
   Depending on the nature of the variable, this may or may not make any physical sense.
 
   ==> Requested Energy <==
 
   ----------------------------------------------------------------------------------------------
               Reaction     Reaction Energy      Error           Reagent 1           Reagent 2
                               Ref     Calc [kcal/mol]              [H] Wt              [H] Wt
   ----------------------------------------------------------------------------------------------
          NBC1-MeMe-3.2     0.0690   1.1639     1.0949     -80.72700202  1     -40.36442840 -2
          NBC1-MeMe-3.3    -0.2390   0.6709     0.9099     -80.72764911  1     -40.36435916 -2
          NBC1-MeMe-3.4    -0.4170   0.3407     0.7577     -80.72806043  1     -40.36430165 -2
          NBC1-MeMe-3.5    -0.5080   0.1244     0.6324     -80.72831099  1     -40.36425461 -2
          NBC1-MeMe-3.6    -0.5410  -0.0129     0.5281     -80.72845373  1     -40.36421659 -2
          NBC1-MeMe-3.7    -0.5390  -0.0961     0.4429     -80.72852567  1     -40.36418623 -2
          NBC1-MeMe-3.8    -0.5150  -0.1430     0.3720     -80.72855247  1     -40.36416227 -2
          NBC1-MeMe-3.9    -0.4800  -0.1659     0.3141     -80.72855167  1     -40.36414365 -2
          NBC1-MeMe-4.0    -0.4390  -0.1733     0.2657     -80.72853498  1     -40.36412938 -2
          NBC1-MeMe-4.1    -0.3960  -0.1712     0.2248     -80.72850993  1     -40.36411859 -2
          NBC1-MeMe-4.2    -0.3540  -0.1633     0.1907     -80.72848118  1     -40.36411044 -2
          NBC1-MeMe-4.3    -0.3150  -0.1525     0.1625     -80.72845143  1     -40.36410422 -2
          NBC1-MeMe-4.4    -0.2790  -0.1403     0.1387     -80.72842215  1     -40.36409932 -2
          NBC1-MeMe-4.6    -0.2170  -0.1155     0.1015     -80.72836761  1     -40.36409177 -2
          NBC1-MeMe-4.8    -0.1680  -0.0933     0.0747     -80.72831991  1     -40.36408563 -2
          NBC1-MeMe-5.0    -0.1300  -0.0747     0.0553     -80.72827951  1     -40.36408021 -2
          NBC1-MeMe-5.4    -0.0800  -0.0479     0.0321     -80.72821875  1     -40.36407122 -2
          NBC1-MeMe-5.8    -0.0500  -0.0312     0.0188     -80.72817678  1     -40.36406353 -2
   ----------------------------------------------------------------------------------------------
            Minimal Dev                         0.0188
            Maximal Dev                         1.0949
        Mean Signed Dev                         0.3509
      Mean Absolute Dev                         0.3509
                RMS Dev                         0.4676
   ----------------------------------------------------------------------------------------------

Thirdly, the :py:func:`~driver_cbs.complete_basis_set` wrapper allows any
compound computational method that can be expressed through :ref:`CBS
<eq:cbs>` to be applied to a molecule while employing the minimum number
of calculations. For example, the job below computes a
triple-quadruple-zeta Helgaker extrapolation of the mp2 correlation energy
atop a quadruple zeta reference to which is appended a double-triple-zeta
Helgaker extrapolated ccsd(t) - mp2 delta correction. Since the mp2 has
been requested through |PSIfour| and the ccsd(t) through Cfour, the
wrapper runs only MP2/cc-pVQZ (P4), CCSD(T)/cc-pVDZ (C4), and
CCSD(T)/cc-pVTZ (C4) single-points. ::

   molecule {
   H 0.0 0.0 0.0
   H 1.0 0.0 0.0
   }
   
   set mp2_type conv

   cbs('mp2', corl_basis='cc-pV[TQ]Z', delta_wfn='c4-ccsd(t)', delta_basis='cc-pV[DT]Z')

This yields::

   ==> CBS <==

   ---------------------------------------------------------------------------------------------------------
       Stage               Method / Basis                                Energy [H]   Scheme
   ---------------------------------------------------------------------------------------------------------
         scf                  scf / cc-pvqz                             -1.10245974   highest_1
        corl                  mp2 / cc-pv[tq]z                          -0.03561890   corl_xtpl_helgaker_2
       delta     c4-ccsd(t) - mp2 / cc-pv[dt]z                           0.03507767   corl_xtpl_helgaker_2
       total                  CBS                                       -1.10300098
   ---------------------------------------------------------------------------------------------------------

Note that especially for :py:func:`~driver_cbs.complete_basis_set`, the
basis set needs to be specified through |mints__basis|, not
|cfour__cfour_basis|.  Many of the wrappers can be used in combination to,
for example, apply a compound method to every molecule in a database or to
optimize a molecule with an extrapolated basis set (findif only for the
moment- analytics coming).

Finally, any number and combination of jobs can be run from a single
|PSIfour| input file.  Depending on the nature of preceding or following
jobs, it is prudent to separate them with the following::

    clean()            # removes Psi4 scratch files
    clean_variables()  # empties the PSI variables list
    cfour {}           # empties the cfour block

.. warning:: Because p4c4 does not inspect the contents of the ``cfour {...}``
   block, once the user specifies a |PSIfour|-style molecule, the
   interface cannot judge whether a sandwich mode (drop the |PSIfour| molecule
   and use the cfour block as the entirety of the ``ZMAT``) or a standard mode
   (translate the |PSIfour| molecule and append additional input from the
   cfour block) is intended. The latter is what actually occurs. If
   there is both a |PSIfour| molecule and a molecule in the cfour block,
   ``ZMAT`` *will* end up with multiple molecules and multiple ``*CFOUR(...)``
   blocks, and it *will not* run.  Therefore, if mixing sandwich and
   standard or pure-\ |PSIfour| computations in an input file, place all
   the sandwich jobs at the beginning before declaring |PSIfour|
   molecules. If necessary, clear the cfour block with ``cfour {}`` before
   commencing standard P4C4 jobs.

.. _`sec:cfourOutput`:

Output
~~~~~~

The output of :program:`xcfour` invoked from a |PSIfour| input file is
written to the |PSIfour| output file as the computation progresses.  If a
Cfour module terminates with a non-zero error code, the value will show up
in :psivar:`CFOUR ERROR CODE <CFOURERRORCODE>`.

.. rubric:: Energies & Scalars

After execution of :program:`xcfour` has completed, the output string is
extensively parsed and appropriate results are stored in :ref:`PSI
Variables <apdx:cfour_psivar>`. All gleaned variables are printed in the
Cfour output section of the |PSIfour| output file, as shown below. ::

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
  //  Cfour c4-ccsd(t) Energy Results  //
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//


  Variable Map:
  ----------------------------------------------------------------------------
  "(T) CORRECTION ENERGY"      =>      -0.007263598030
  "CCSD CORRELATION ENERGY"    =>      -0.275705492359
  "CCSD TOTAL ENERGY"          =>     -76.338453952539
  "CCSD(T) CORRELATION ENERGY" =>      -0.007263598030
  "CCSD(T) TOTAL ENERGY"       =>     -76.345717550569
  "CFOUR ERROR CODE"           =>       0.000000000000
  "CURRENT CORRELATION ENERGY" =>      -0.007263598030
  "CURRENT ENERGY"             =>     -76.345717550569
  "CURRENT REFERENCE ENERGY"   =>     -76.062748460180
  "MP2 CORRELATION ENERGY"     =>      -0.270191667755
  "MP2 OPPOSITE-SPIN ENERGY"   =>      -0.204890356651
  "MP2 SAME-SPIN ENERGY"       =>      -0.065301311104
  "MP2 TOTAL ENERGY"           =>     -76.332940127935
  "NUCLEAR REPULSION ENERGY"   =>       9.187331653300
  "SCF TOTAL ENERGY"           =>     -76.062748460180

The PSI Variables are also available from the input file for manipulation.
For instance, to compute the MBPT 2 3/4 energy from MBPT 3 results, the
following could be used. ::

   energy('c4-mp3')
   mp2p75_corl = 0.75 * get_variable('mp3 correlation energy') + \
                 0.25 * get_variable('MP2 correlation energy')
   print mp2p75_corl + get_variable('scf total energy')

.. caution:: Some features are not yet implemented. Buy a developer a coffee.

   - No PSI Variables for properties: *e.g.*, :psivar:`SCF DIPOLE X<SCFDIPOLEX>`

   - No PSI Variables for excited state energies

   The formation of further regexes for properties, excited states, etc.
   is one of the primary areas in which this interface requires further
   work.

.. rubric:: Gradients and Arrays

In addition to parsing the output stream, results are collected from files
written to the scratch directory. Presently, the ``GRD`` file is parsed
and printed to the output file, as shown below. Also printed is the Cfour
gradient after manipulation by the P4C4 interface and used by |PSIfour|
going forward. Manipulation is necessary because Cfour determinedly uses
its own internal orientation and atom ordering while |PSIfour| naturally
expects the gradient to be aligned with the active molecule. The geometry
in ``GRD`` and the geometry of |PSIfours| active molecule are shifted,
rotated, flipped, and otherwise badgered into coincidence, then the same
manipulations are applied to the gradient in ``GRD``, the result of which
is printed below and passed on to Optking. ::

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
  //   Cfour c4-scf Gradient Results   //
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//

  ...

  Irrep: 1 Size: 4 x 3 

                 1            2            3     

    1    0.0000000   -0.0122978    0.0000000
    2   -0.0051192    0.0040993   -0.0088667
    3   -0.0051192    0.0040993    0.0088667
    4    0.0102384    0.0040993    0.0000000


  CFOUR scratch file GRD has been read
    4        0.0000000000
        7.0000000000       -0.0880964705        0.0000000000        0.0000000000
        1.0000000000        0.4080144050       -0.9178691296       -1.5897959670
        1.0000000000        0.4080144050       -0.9178691296        1.5897959670
        1.0000000000        0.4080144050        1.8357382590        0.0000000001
        7.0000000000       -0.0122978407        0.0000000000        0.0000000000
        1.0000000000        0.0040992802       -0.0051191833       -0.0088666856
        1.0000000000        0.0040992802       -0.0051191833        0.0088666856
        1.0000000000        0.0040992802        0.0102383666        0.0000000000

The gradient can also be accessed from the input file as a
:py:class:`~psi4.core.Matrix` object through
:py:func:`psi4.get_gradient`.

.. rubric:: Cfour Files

The contents of all files associated with Cfour are accessible from the
input file through the Python dictionary ``P4C4_INFO``. That is,
``P4C4_INFO['zmat']`` returns a string of the input file sent to Cfour.
Accessible arguments are ``zmat``, ``output``, and any that have been
produced of ``grd``. For example, to print to the screen if CC convergence
is reached, the following could be placed in the |PSIfour| input file. ::

   energy('c4-ccsd')
   print 'miracle?', 'miracle' in P4C4_INFO['output']

.. rubric:: Scratch Files

By default, a separate subdirectory for each Cfour call is created within
the job's scratch directory. To explicitly specify the location of the
Cfour scratch, execute with, for example, ``energy('cfour',
path='/full/path/to/cfour/scratch')``. Regardless of whether the location
is specified or default, whether to preserve the scratch directory after
the computation can be specified with ``energy('cfour', keep=True)`` or
(the default) ``energy('cfour', keep=False)``. *path* and *keep* are
keyword arguments that get interpreted by the
:py:func:`~procedures.interface_cfour.run_cfour` function documented below.

.. autofunction:: procedures.interface_cfour.run_cfour(name [, keep, path])

.. _`sec:cfourFunctionality`:

Functionality
~~~~~~~~~~~~~

Through clever use of the ``cfour {...}`` block, one could run most any
Cfour computation through the P4C4 interface.  In contrast, enumerated
below are tested functionalities where results from Cfour are collected
into |PSIfour| data objects.

.. rubric:: Implemented

* Single-point :py:func:`~psi4.energy` commands for :ref:`ground state
  methods <table:energy_cfour>`. Examples:
  :srcsample:`cfour/sp-rhf-ccsd_t_-ao-ecc`, :srcsample:`cfour/scf4`,
  :srcsample:`cfour/mints5`.

* Analytic :py:func:`~psi4.gradient` and :py:func:`~psi4.optimize`
  commands for :ref:`ground state methods <table:grad_cfour>`. Real and
  Ghost atoms permitted (though the latter will naturally collapse after
  several cycles). Examples: :srcsample:`cfour/opt-rhf-ccsd_t_`,
  :srcsample:`cfour/mp2-1`, and :srcsample:`cfour/mints5-grad`.

.. warning:: There exist molecules (*e.g.*, allene) where the
   inertial frame is not unique (planes along atoms or between
   atoms). The orientation reconciling machinery currently does not
   handle these cases and will fail with "Axis unreconcilable between
   QC programs". I will get to this soon.

* Finite difference of energy :py:func:`~psi4.gradient` and
  :py:func:`~psi4.optimize` for :ref:`methods <table:energy_cfour>`.
  Force with ``gradient('name', dertype=0)``, *etc.*.

* :py:func:`~driver_nbody.nbody_gufunc` for computation of interaction energies with or
  without counterpoise correction. Example: :srcsample:`cfour/dfmp2-1`.

* :py:func:`~wrapper_database.database` for computation of a collection of molecules in a
  single input, with summarization of results. Examples:
  :srcsample:`cfour/pywrap-db1` and :srcsample:`cfour/psi-a24-grad`.

* :py:func:`~driver_cbs.complete_basis_set` for computation of compound methods involving
  basis set extrapolations and/or delta corrections with any combination
  of |PSIfour| and Cfour computational methods and |PSIfour| basis sets.
  Example: :srcsample:`cfour/pywrap-cbs1`.

.. rubric:: Not Yet Implemented

* Ground state CI energies and optimizations

* Excited state energies and optimizations

* Properties are not yet regex-ed, transformed into input frame, and
  stowed in PSI Variables.

* Property calls that required extra computation not yet translated into
  :py:func:`~psi4.property` computation command

* Frequencies

Energy methods available through P4C4 interface

.. include:: cfour_table_energy.rst

Gradient methods available through P4C4 interface

.. include:: cfour_table_grad.rst

.. _`table:cfour_cc_program`:

.. notes on preferred modules from JFS
.. comment always abcdtype = aobasis (but sometimes 
.. comment ncc does rhf ccsdt(q)
.. comment reccommended code to do with, not only code (b/c mrcc can do much of this)
.. comment .. table:: Cfour coupled-cluster program defaults by calculation type
.. comment 
.. comment     +-----------------------------------------+---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         |                                 |                       | RHF    | UHF    | ROHF    |
.. comment     |                                         |                                 |                       +--------+--------+---------+
.. comment     | Driver Call, |cfour__cfour_deriv_level| | name, |cfour__cfour_calc_level| | |cfour__cfour_excite| | |cfour__cfour_cc_program| |
.. comment     +=========================================+=================================+=======================+========+========+=========+
.. comment     | :py:func:`~psi4.energy`, zero           | cc2                             | none                  | vcc    | vcc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | _cc    | _cc    | _cc     |
.. comment     |                                         |vcc for everything               +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | _cc    | _cc    | _cc     |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | ccsd                            | none                  | ecc    | ecc    | ecc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | vcc    | vcc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | vcc    | vcc    | vcc     |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | ccsd(t)                         | none                  | ecc    | ecc    | ecc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 |        |        |         |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           |        |        |         |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | cc3                             | none                  | vcc    | vcc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | vcc    | vcc    |         |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           |        |        |         |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | ccsdt (no aobasis sp or opt, also prob for cc3, or for eomea grad) 
.. comment ecc / vcc / vcc
.. comment | none                  | ecc    | ecc    | ecc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | ecc    | mrcc   | mr_cc   |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomeano/eomipyes      | ecc    |        |         |
.. comment     +-----------------------------------------+---------------------------------+-----------------------+--------+--------+---------+
.. comment     | :py:func:`~psi4.optimize`, first        | cc2                             | none                  | vcc    | vcc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | vcc    | vcc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip not sure ask  | vcc  | vcc  | vcc     |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | ccsd                            | none                  | ecc    | ecc    | ecc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | vcc    | vcc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | vcc    | vcc    | vcc     |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | ccsd(t)                         | none                  | ecc    | ecc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 |        |        |         |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           |        |        |         |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | cc3                             | none                  | vcc    |        |         |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | vcc    |        |         |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | vcc    |        |         |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | ccsdt                           | none                  | ecc    | mrcc   | mrc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | mrcc   | mrcc   | mrcc    |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | mrcc   | mrcc   | mrcc    |
.. comment     +-----------------------------------------+---------------------------------+-----------------------+--------+--------+---------+
.. comment properties same as grad
.. comment 2nd deriv ecc only for ccsd(t)

Specification Details
~~~~~~~~~~~~~~~~~~~~~

The above narrative introduction to the P4C4 interface should be
sufficient to get started. Issues of competition between |PSIfour| and
Cfour specification format are generally resolved behind the scenes:
not according to a *simple* rule but according to sensible, sometimes
intricate, rules governed by user intent (and integration of Cfour to
behave like a |PSIfour| module). Much can be gleaned by just running
inputs and inspecting the ``ZMAT`` passed to Cfour, but when questions
arise, here are the specifics, the governing laws.

* Specifying a piece of input in |PSIfour| format is entering into
  a contract that you mean it. In particular this applies to
  molecule (including charge/multiplicity through :samp:`molecule
  {optional_molecule_name} \\{...\\}`), memory (through :samp:`memory
  {value} {unit}`), computational method (through . If Cfour keywords
  are specified with values that contradict the |PSIfour| input,
  execution is halted.

  As an example, the input below is set up to fail in four ways:
  contradictory specification of memory, multiplicity, computational
  method, and derivative level. Note, though, that the ``cfour_units
  angstrom`` setting is permissible, since it concurs with the value
  implied in the molecule block. ::

    memory 300 mb

    molecule {
    H
    H 1 0.7
    }

    set basis 6-31g
    set cfour_multiplicity 3         # clash with implicit singlet in molecule {} above
    set cfour_units angstrom         # no problem, consistent with molecule {} above
    set cfour_memory_size 100000000  # clash with 300 mb above
    set cfour_calc_level ccsd        # clash with 'c4-scf' below
    set cfour_deriv_level first      # clash with energy() below (use gradient('c4-scf') to achieve this)

    energy('c4-scf')

* Specifying anything in |PSIfour| format (molecule, basis, options,
  method call) starts building a ``*CFOUR(...)`` directive for the
  ``ZMAT`` file. Since the contents of the ``cfour {...}`` block are
  blindly appended to any input interpreted from |PSIfour| format, mixing
  of |PSIfour| and Cfour input formats likely *will* give rise to multiple
  ``*CFOUR(...)`` directives in the prospective ``ZMAT``, execution of
  which *will* be trapped and halted.  Proper uses for the ``cfour {...}``
  block are for the sandwich mode, where the entire ``ZMAT`` is enclosed,
  or for extra directives like ``%excite*``, which presently have no other
  specification route.

* Specifying the basis is perhaps the regulated piece of input. Since
  basis set names differ between |PSIfour| and Cfour and it's not
  practical to compare exponent-to-exponent, any input file with both
  |mints__basis| and |cfour__cfour_basis| keywords present will halt. Once
  a basis set has been requested through |mints__basis|, overriding the
  default spherical/Cartesian setting must be done through
  |globals__puream| (as opposed to |cfour__cfour_spherical|).

* Specifying keywords that control geometry optimization is
  straightforward. Unless the optimization is invoked in sandwich mode,
  all Cfour optimization keywords (*e.g.*, |cfour__cfour_geo_maxcyc|) are
  ineffective, as the Cfour optimizer is never invoked. |PSIfour|
  optimization keywords (*e.g.*, |optking__geom_maxiter|) instead fill
  these roles.

* Specifying the computational method (through, for instance,
  ``energy('c4-ccsd')`` instead of ``energy('cfour')``) often
  sets additional keywords consistent with best practices (*e.g.*,
  |cfour__cfour_cc_program|). Since those settings are implicit, any
  explicit setting of those those keywords, whether contradicting or
  concurring, takes priority (halts never generated). The following are
  some concrete examples. For the moment, click the source button at
  :py:func:`qcdb.cfour.muster_modelchem` for details of what keywords
  get set.

  * runs in vcc since that's Cfour's default for cc_program ::

      set cfour_calc_level ccsd
      energy('cfour')

  * runs in ecc since Cfour's default overwritten by keyword ::

      set cfour_calc_level ccsd
      set cfour_cc_program ecc
      energy('cfour')

  * runs in ecc since that's best practice for the requested ccsd ::

      energy('c4-ccsd')

  * runs in vcc since *hidden* default overwritten by keyword ::

      set cfour_cc_program vcc
      energy('c4-ccsd')

* Specifying certain keywords that are nominally applicable for pure-\
  |PSIfour| modules directs them to fulfil analogous roles
  in the Cfour program (*e.g.*, |scf__maxiter| is used to set
  |cfour__cfour_scf_maxcyc|). This keyword translation only takes place
  if the keywords are explicitly set in the input file (part of that
  contract that you mean it), meaning that |PSIfours| defaults don't
  get imposed on Cfour. Also, in the case where a translatable pure-\
  |PSIfour| keyword and its translation Cfour keyword are both set,
  the value attached to the latter is always used. Below are a few
  clarifying examples.

  * uses :math:`10^{-7}` SCF conv crit since that's Cfour's default
    for |cfour__cfour_scf_conv| ::

      energy('c4-scf')

  * uses :math:`10^{-6}` SCF conv crit since default overwritten by
    keyword ::

      set cfour_scf_conv 6
      energy('c4-scf')

  * uses :math:`10^{-5}` SCF conv crit since default overwritten by
    :ref:`SCF module<apdx:scf>` keyword ::

      set d_convergence 5
      energy('c4-scf')

  * uses :math:`10^{-6}` SCF conv crit since default overwritten by
    :ref:`SCF module<apdx:scf>` keyword (local scope works, too) where
    the |PSIfours| more flexible float input has been rounded down to
    the integer required by Cfour ::

      set scf d_convergence 5e-6
      energy('c4-scf')

  * uses :math:`10^{-6}` SCF conv crit since default overwritten
    and Cfour module keyword trumps |PSIfour| SCF module keyword ::

      set cfour_scf_conv 6
      set d_convergence 8
      energy('c4-scf')

  The keyword translation feature is still in the proof-of-principle
  stage, so only a handful (found here) of keywords participate.

.. note:: Longtime Cfour users who may consider this keyword
   translation a flaw rather than a feature can avoid it entirely by
   confining keywords to the :ref:`Cfour module<apdx:cfour>` along with
   |mints__basis| and |globals__puream| (opt, too?)

Misc. Running
~~~~~~~~~~~~~

Naturally, in |PSIfour| multiple jobs can be run in succession from the input file.

Control optimizations with optking keywords HERE. Cfour ``GRD`` file is
written to |PSIfour| output file. Gradient transformed back into the frame
in which it was shipped off to Cfour is also written to the |PSIfour|
output file and is available from input as :py:func:`~psi4.get_gradient`.

sandwich mode := molecule and cfour list within
Naturally, additional jobs can follow in the input file.
Depending on the nature of preceding or following jobs, it is prudent to
separate them with the following::

    clean()            # removes Psi4 scratch files
    clean_variables()  # empties the PSI variables list
    cfour {}           # empties

In this scheme, the contents of the ``cfour {...}`` block are tacked onto
the end of the ``ZMAT`` file that is otherwise written from psi style
format. It is by this route that, for example ``%excite*`` sections can at
present be specified.

The execution of :program:`xcfour` can be modified by a few parameters.  Setting
the option |cfour__cfour_omp_num_threads| sets the environment variable
:envvar:`OMP_NUM_THREADS` for only the duration of the Cfour computation.
That is, portions of an input file that run |PSIfour| modules are
unaffected.  Additionally, there are a few arguments to the function
:py:func:`~procedures.interface_cfour.run_cfour` that control the Cfour scratch
directory. 

.. comment Notes to Self
.. comment ~~~~~~~~~~~~~
.. comment 
.. comment Test checked-in GENBAS on installed copy
.. comment 
.. comment Reference still not factored into cc_program!
.. comment 
.. comment optimize on a sandwich calc? errors out
.. comment 
.. comment 
.. comment .. _`table:cfour_cc_program`:
.. comment 
.. comment .. table:: Cfour coupled-cluster program defaults by calculation type
.. comment 
.. comment     +-----------------------------------------+---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         |                                 |                       | RHF    | UHF    | ROHF    |
.. comment     |                                         |                                 |                       +--------+--------+---------+
.. comment     | Driver Call, |cfour__cfour_deriv_level| | name, |cfour__cfour_calc_level| | |cfour__cfour_excite| | |cfour__cfour_cc_program| |
.. comment     +=========================================+=================================+=======================+========+========+=========+
.. comment     | :py:func:`~psi4.energy`, zero           | cc2                             | none                  | vcc    | vcc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | _cc    | _cc    | _cc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | _cc    | _cc    | _cc     |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | ccsd                            | none                  | ecc    | ecc    | ecc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | ecc    | _cc    | _cc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | _cc    | _cc    | _cc     |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | ccsd(t)                         | none                  | ecc    | ecc    | ecc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | _cc    | _cc    | _cc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | _cc    | _cc    | _cc     |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | cc3                             | none                  | vcc    | vcc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | _cc    | _cc    | _cc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | _cc    | _cc    | _cc     |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | ccsdt                           | none                  | ecc    | ecc    | ecc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | _cc    | _cc    | _cc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | _cc    | _cc    | _cc     |
.. comment     +-----------------------------------------+---------------------------------+-----------------------+--------+--------+---------+
.. comment     | :py:func:`~psi4.optimize`, first        | cc2                             | none                  | _cc    | _cc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | _cc    | _cc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | _cc    | _cc    | vcc     |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | ccsd                            | none                  | _cc    | _cc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | _cc    | _cc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | _cc    | _cc    | vcc     |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | ccsd(t)                         | none                  | ecc    | _cc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | _cc    | _cc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | _cc    | _cc    | vcc     |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | cc3                             | none                  | _cc    | _cc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | _cc    | _cc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | _cc    | _cc    | vcc     |
.. comment     |                                         +---------------------------------+-----------------------+--------+--------+---------+
.. comment     |                                         | ccsdt                           | none                  | ecc    | _cc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomee                 | _cc    | _cc    | vcc     |
.. comment     |                                         |                                 +-----------------------+--------+--------+---------+
.. comment     |                                         |                                 | eomea/eomip           | _cc    | _cc    | vcc     |
.. comment     +-----------------------------------------+---------------------------------+-----------------------+--------+--------+---------+
.. comment 
.. comment 
