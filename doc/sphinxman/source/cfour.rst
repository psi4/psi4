
.. include:: autodoc_abbr_options_c.rst

.. index:: Cfour
.. _`sec:cfour`:

Interface to CFOUR by J. Stanton and J. Gauss
=============================================

.. codeauthor:: Lori A. Burns
.. sectionauthor:: Lori A. Burns

*Module:* :ref:`Keywords <apdx:cfour>`, :ref:`PSI Variables <apdx:cfour_psivar>`, :source:`CFOUR <src/bin/cfour>`, :ref:`Samples <apdx:testSuitecfour>`

|PSIfour| contains code to interface to the Cfour quantum chemistry suite of
John F. Stanton (U. Texas, Austin) and J\ |u_dots|\ rgen Gauss (U. Mainz),
which is available after a license agreement from 
`http://www.cfour.de/ <http://www.cfour.de/>`_.

Installation
~~~~~~~~~~~~

Follow the instructions provided with the Cfour download to install the
executable or to build the source. To by used by |PSIfour|, the program
binary (``xcfour``) must be found in your :envvar:`PATH` or
:envvar:`PSIPATH`.  The ``GENBAS`` file containing basis sets in Cfour
format is not necessary for this interface, but if you prefer to access
basis sets the "Cfour way", ``GENBAS``, too, must be in :envvar:`PATH` or
:envvar:`PSIPATH`. If |PSIfour| is unable to execute the binary, an error
will be reported.

.. warning:: The p4c4 interface isn't in the master branch nor will it be in
   the near future. To run this code, (1) build the ``c4`` branch of psi4, (2) 
   find a copy of cfour and put it in :envvar:`PATH` or :envvar:`PSIPATH`, and
   (3) clone https://github.com/loriab/qcdb.git python module and prepend
   :envvar:`PYTHONPATH` with the top qcdb directory (the path added to
   PYTHONPATH should have one qcdb in it; the cloned qcdb is what needs to be
   imported in preference to the one already in psi4). Execute psi4 as usual.

NOTE: Test checked-in GENBAS on installed copy

Cfour for |PSIfour| Users
~~~~~~~~~~~~~~~~~~~~~~~~~

* Set memory as usual

* Set molecule as usual

* Set basis set as usual (Cfour only cares about orbital basis, no fitting
  bases)

* For the type of computation intended, find appropriate options at
  :ref:`Keywords <apdx:cfour>`. These keyword summaries contain the same information as the
  `proper CFOUR options list
  <http://slater.chemie.uni-mainz.de/cfour/index.php?n=Main.ListOfKeywordsInAlphabeticalOrder>`_
  plus notes on keyword relevance when run through |PSIfour|.  Information
  at the `CFOUR manual
  <http://slater.chemie.uni-mainz.de/cfour/index.php?n=Main.Manual>`_ may
  also be useful, as may the many samples at :source:`samples/cfour`.

* Generally, the p4c4 interface will handle best practices for path of
  execution: ``vcc``/``ecc``, derivative type, *etc.* The user is
  still responsible for setting convergence, frozen core, guess, diis,
  *etc.*

* Set Cfour keywords just like |PSIfour| keywords. The names of keywords
  are unchanged beyond a prepended "cfour\_". (Though be aware that common
  abbreviations like CALC and REF must be fully spelled out when used in
  |PSIfour|.)

* Set the task as usual, indicating Cfour as the intended code by
  prepending "c4-" to the method argument. So ``energy('scf')`` becomes
  ``energy('c4-scf')`` and ``optimize('ccsd(t)')`` becomes
  ``optimize('c4-ccsd(t)')``.


|PSIfour| for Cfour Users
~~~~~~~~~~~~~~~~~~~~~~~~~

In the simplest version of the Psi4/Cfour interface, a |PSIfour| input
file can simply "wrap" a ``ZMAT`` file and execute ``xcfour``. This is
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

Here, the contents of the ``cfour {...}`` block are written directly
to a ``ZMAT`` file. This is joined by a default ``GENBAS`` file
(:source:`lib/basis/GENBAS`).  To preferentially use your own ``GENBAS``,
place it in :envvar:`PATH` or :envvar:`PSIPATH`. The line calling
:py:func:`~driver.energy` with argument ``'cfour'`` invokes :program:`xcfour`.

After execution of the ``energy('cfour')`` line completes, Cfour results
are read back into |PSIfour| format and are thereafter accessible for
further processing in the input file. See :ref:`sec:cfouroutput` for
details. This in-memory variable and array storage of results (as opposed to only
file storage) is the only advantage thus far incurred by the p4c4
interface. We'll call this mode of basic utility the "sandwich" mode.

Molecule specification in |PSIfour| allows Cartesians, Z-matrices, mixed
Cartesian/Z-matrix, negation of variables, delayed specification of
variables, etc., all in a whitespace-tolerant format. See
:ref:`sec:moleculeSpecification` for details and
:srcsample:`cfour/psi-mints5` for examples. When a |PSIfour|-style
molecule is supplied, its geometry is written to ``ZMAT`` in Cartesian
form and the |cfour__cfour_coordinates| =CARTESIAN, |cfour__cfour_units|
=ANGSTROM, |cfour__cfour_charge|, and |cfour__cfour_multiplicity| keywords
are set appropriately in the ``*CFOUR(...)`` directive.

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
for the keyword commands. Specifcation of strings and integers requires no
translation; :ref:`booleans <op_c_boolean>` have extended freedom of
format; arrays must be translated into Python-style (square-bracket
bounded and comma delimited) of appropriate dimension. There are many
sample inputs in :source:`tests/cfour/` starting with ``sp-`` that take
examples from the Cfour manual and first run them in sandwich mode and
then run them as translated into |PSIfour| format.

.. note:: |PSIfour| only recognizes keywords by their full name, so the common
   Cfour keyword abbreviations CALC, REF, etc. must be replaced by their
   proper names of CALC_LEVEL, REFERENCE, etc.

Whenever the molecule is supplied in |PSIfour| format, it is possible to
perform geometry optimizations where Cfour supplies the gradient and the
|PSIfour| module :ref:`optking <sec:optking>` drives the structural
changes. Because of the limitations on geometry specification for
optimizations in Cfour, optking-driven optimizations are the *only*
optimizations allowed in the p4c4 interface. (The exception is sandwich
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

Note that the primary change is the exchange of :py:func:`~driver.energy`
for :py:func:`~driver.optimize` to trigger an optimization.
Setting |optking__g_convergence| =CFOUR provides a good imitation
of Cfour default convergence criteria.  Several sample inputs
in :source:`test/cfour/` starting with ``opt-`` show basic
geometry optimizations. :srcsample:`cfour/psi-mints5-grad`
shows optimizations from a variety of molecule input formats, and
:srcsample:`cfour/psi-ghost-grad` shows an optimization with ghosted
atoms.

The above example also shows the total memory for the computation being
set in |PSIfour| format. See :ref:`sec:memory` for details. When
specified, the memory value is passed on to Cfour by setting keywords
|cfour__cfour_memory_size| and |cfour__cfour_mem_unit| =MB.

|PSIfour| has an extensive :ref:`basis set library <apdx:basisElement>` in
Gaussian94 format. See :ref:`sec:basisBuiltIn` for details.  Contrasts to
Cfour basis handling include identifying basis sets by standard name
(aug-cc-pVDZ instead of AUG-PVDZ), direct handles for diffuse-fn-pruned
sets (*e.g.*, jun-cc-pVDZ), case insensitivity, appropriate setting of
spherical/Cartesian depending on basis set design, and syntax to set
different basis sets to different classes of atoms without listing each
atom. All of these features are available to Cfour by using the
|mints__basis| keyword instead of |cfour__cfour_basis| (accompanied, of
course, by specifying the molecule |PSIfour|-style). Internally, |PSIfour|
processes the basis set as usual, then translates the basis set format and
writes out a ``GENBAS`` file with an entry for each atom. The p4c4
interface sets keyword |cfour__cfour_basis| =SPECIAL and
|cfour__cfour_spherical| as appropriate, then writes the basis section
necessary for SPECIAL below the ``*CFOUR(...)`` block. (I'm sorry that the
name of the basis doesn't appear in ``ZMAT``, but the combination of the
~14 character basis name limit and the absense of a comment line marker
rather preclude that helpful label.)

* other names to energy, optimize
* multiple jobs in file
* ref translation, more to come
* wrappers (must use psi4 basis for cbs)

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
   commencing standard p4c4 jobs.


sandwich mode := molecule and cfour list within
warning to a lesser extent with options
Naturally, additional jobs can follow in the input file.
Depending on the nature of preceeding or following jobs, it is prudent to
separate them with the following::

    clean()            # removes Psi4 scratch files
    clean_variables()  # empties the PSI variables list
    cfour {}           # empties

In this scheme, the contents of the ``cfour {...}`` block are tacked onto
the end of the ``ZMAT`` file that is otherwise written from psi style
format. It is by this route that, for example ``%excite*`` sections can at
present be spcified.

One key difficulty to consider in constructing an interface, and
particularly one such as this where input form is liberal, between cfour
and psi4 formats, is that specificity on how info should be conveyed is
lost. It is possible to specify something multiple ways, each of which
are legal, then what happens if multiples ones are used. Does the last
one get used? the psi4 one? does execution notify or stop?

The above narrative introduction to the P4C4 interface should be
sufficient to get started. Issues of competition between |PSIfour| and
Cfour specification format are generally resolved behind the scenes:
not according to a *simple* rule but according to sensible, sometimes
intricate, rules governed by user intent (and integration of Cfour to
behave like a |PSIfour| module). Much can be gleaned by just running
inputs and inspecting the ``ZMAT`` passed to Cfour, but when questions
arise, here are the specifics.

Some Governing Laws

* opt!!
* puream

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

* Specifying the basis is perhaps the regulated piece of input. Since
  basis set names differ between |PSIfour| and Cfour and it's not
  practical to compare exponent-to-exponent, any input file with both
  |mints__basis| and |cfour__cfour_basis| keywords present will halt.

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


* Specifying anything in |PSIfour| format (molecule, basis, options,
  method call) starts building a ``*CFOUR(...)`` directive for the
  ``ZMAT`` file. Since the contents of the ``cfour {...}`` block are
  blindly appended to any input interpreted from |PSIfour| format,
  mixing of |PSIfour| and Cfour input formats likely *will* give rise
  to multiple ``*CFOUR(...)`` directives in the prospective ``ZMAT``,
  execution of which *will* be trapped and halted.

  Proper uses for the ``cfour {...}`` block are for the sandwich mode,
  where the entire ``ZMAT`` is enclosed, or for extra directives like
  ``%excite*``, which presently have no other specification route.

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


Since these overridings of Cfour defaults 


.. molecule {
.. 0 2
.. N
.. H 1 R
.. H 1 R 2 A
.. 
.. R=1.008
.. A=105.0
.. }
.. 
.. set {
.. cfour_CALC_level=HF
.. cfour_BASIS=qz2p
.. cfour_REFerence=UHF
.. cfour_occupation [[3,1,1,0],[3,0,1,0]]
.. cfour_SCF_CONV=12
.. }





Naturally, in |PSIfour| multiple jobs can be run in succession from the input file.


The execution of ``xcfour`` can be modified by a few parameters.  Setting
the option |cfour__cfour_omp_num_threads| sets the environment variable
:envvar:`OMP_NUM_THREADS` for only the duration of the Cfour computation.
That is, portions of an input file that run |PSIfour| modules are
unaffected.  Additionally, there are a few arguments to the function
:py:func:`~interface_cfour.run_cfour` that control the Cfour scratch
directory. By default, a separate subdirectory for each Cfour call is
created within the job's scratch directory. To explicitly specify the
location of the the Cfour scratch, execute with, for example,
``energy('cfour', path='/full/path/to/cfour/scratch')``. Regardless of
whether the location is specified or default, whether to preserve the
scratch directory after the computation can be specified with
``energy('cfour', keep=True)`` or (the default) ``energy('cfour',
keep=False)``.






optimize on a sandwich calc?

.. _`sec:cfouroutput`:

Output
~~~~~~

The output of ``xcfour`` invoked from a |PSIfour| input file is written to
the |PSIfour| output file as the computation progresses. Additionally, the
output string is extensively parsed and appropriate results are stored in
:ref:`PSI Variables <apdx:cfour_psivar>`. The formation of further regexes
for properties, excited states, etc. is one of the primary areas in which
this interface requires further work. If a Cfour module terminates with a
non-zero error code, the value will show up in :psivar:`CFOUR ERROR CODE
<CFOURERRORCODE>`.

In addition to parsing the output stream, results are collected from files
written to the scratch directory. Presently, the ``GRD`` file is
parsed 

tilde for geom opt

.. :py:func:`~interface_cfour.run_cfour`

.. autofunction:: interface_cfour.run_cfour(name [, keep, path])

Functionality
~~~~~~~~~~~~~

.. include:: cfour_table_energy.rst

.. include:: cfour_table_grad.rst

.. _`table:cfour_cc_program`:

.. table:: Cfour coupled-cluster program defaults by calculation type

    +-----------------------------------------+---------------------------------+-----------------------+--------+--------+---------+
    |                                         |                                 |                       | RHF    | UHF    | ROHF    |
    |                                         |                                 |                       +--------+--------+---------+
    | Driver Call, |cfour__cfour_deriv_level| | name, |cfour__cfour_calc_level| | |cfour__cfour_excite| | |cfour__cfour_cc_program| |
    +=========================================+=================================+=======================+========+========+=========+
    | :py:func:`~driver.energy`, zero         | cc2                             | none                  | vcc    | vcc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomee                 | _cc    | _cc    | _cc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomea/eomip           | _cc    | _cc    | _cc     |
    |                                         +---------------------------------+-----------------------+--------+--------+---------+
    |                                         | ccsd                            | none                  | vcc    | vcc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomee                 | ecc    | _cc    | _cc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomea/eomip           | _cc    | _cc    | _cc     |
    |                                         +---------------------------------+-----------------------+--------+--------+---------+
    |                                         | ccsd(t)                         | none                  | ecc    | ecc    | ecc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomee                 | _cc    | _cc    | _cc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomea/eomip           | _cc    | _cc    | _cc     |
    |                                         +---------------------------------+-----------------------+--------+--------+---------+
    |                                         | cc3                             | none                  | vcc    | vcc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomee                 | _cc    | _cc    | _cc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomea/eomip           | _cc    | _cc    | _cc     |
    |                                         +---------------------------------+-----------------------+--------+--------+---------+
    |                                         | ccsdt                           | none                  | ecc    | ecc    | ecc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomee                 | _cc    | _cc    | _cc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomea/eomip           | _cc    | _cc    | _cc     |
    +-----------------------------------------+---------------------------------+-----------------------+--------+--------+---------+
    | :py:func:`~driver.optimize`, first      | cc2                             | none                  | _cc    | _cc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomee                 | _cc    | _cc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomea/eomip           | _cc    | _cc    | vcc     |
    |                                         +---------------------------------+-----------------------+--------+--------+---------+
    |                                         | ccsd                            | none                  | _cc    | _cc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomee                 | _cc    | _cc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomea/eomip           | _cc    | _cc    | vcc     |
    |                                         +---------------------------------+-----------------------+--------+--------+---------+
    |                                         | ccsd(t)                         | none                  | ecc    | _cc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomee                 | _cc    | _cc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomea/eomip           | _cc    | _cc    | vcc     |
    |                                         +---------------------------------+-----------------------+--------+--------+---------+
    |                                         | cc3                             | none                  | _cc    | _cc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomee                 | _cc    | _cc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomea/eomip           | _cc    | _cc    | vcc     |
    |                                         +---------------------------------+-----------------------+--------+--------+---------+
    |                                         | ccsdt                           | none                  | ecc    | _cc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomee                 | _cc    | _cc    | vcc     |
    |                                         |                                 +-----------------------+--------+--------+---------+
    |                                         |                                 | eomea/eomip           | _cc    | _cc    | vcc     |
    +-----------------------------------------+---------------------------------+-----------------------+--------+--------+---------+

