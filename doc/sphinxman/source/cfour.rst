
.. include:: autodoc_abbr_options_c.rst

.. index:: Cfour
.. _`sec:cfour`:

Interface to CFOUR by J. Stanton and J. Gauss
=============================================

.. codeauthor:: Lori A. Burns
.. sectionauthor:: Lori A. Burns

*Module:* :ref:`Keywords <apdx:cfour>`, :ref:`PSI Variables <apdx:cfour_psivar>`, :source:`CFOUR <src/bin/cfour>`

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

NOTE: Need to check in a GENBAS so tests can run?

Cfour for |PSIfour| Users
~~~~~~~~~~~~~~~~~~~~~~~~~

* Set memory as usual

* Set molecule as usual

* Set basis set as usual (Cfour only cares about orbital basis, no fitting
  bases)

* For the type of computation intended, find appropriate options at
  :ref:`Keywords <apdx:cfour>` which contains the same information as the
  `proper CFOUR options list
  <http://slater.chemie.uni-mainz.de/cfour/index.php?n=Main.ListOfKeywordsInAlphabeticalOrder>`_
  plus notes on keyword relevance when run through |PSIfour|.  Information
  at the `CFOUR manual
  <http://slater.chemie.uni-mainz.de/cfour/index.php?n=Main.Manual>`_ may
  also be useful, as may the many samples at :source:`samples/cfour`.

* Generally, the p4c4 interface will handle best practices for path of
  execution: ``vcc``/``ecc``, derivative type, *etc.* Whereas the user is
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

Here, the contents of the ``cfour {...}`` block are written directly to a
``ZMAT`` file. This is joined by a default ``GENBAS`` file (to use your
own, place it in :envvar:`PATH` or :envvar:`PSIPATH`). The line calling
:py:func:`~driver.energy` with argument ``'cfour'`` invokes ``xcfour``.

After execution of the ``energy('cfour')`` line completes, Cfour results
are read back into |PSIfour| format and are thereafter accessible for
further processing in the input file. See :ref:`sec:cfouroutput` for
details. This variable and array storage of results (as opposed to only
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

Whenever the molecule is supplied in |PSIfour| format, the job control
keywords must be too. All :ref:`Cfour keywords <apdx:cfour>` are the usual
ones, prepended by ``cfour_`` to avoid any possible name conflicts.  As
detailed in :ref:`sec:jobControl`, setting keywords is, too, flexible in
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
for :py:func:`~driver.optimize` to trigger an optimization.  Setting
|optking__g_convergence| =CFOUR provides a good imitation of Cfour default
convergence criteria.  Several sample inputs in :source:`test/cfour/`
starting with ``opt-`` show basic geometry optimizations.

The example above also shows the total memory for the computation being
set in |PSIfour| format. See :ref:`sec:memory` for details. When
specified, the memory value is passed on to Cfour by setting keywords
|cfour__cfour_memory_size| and |cfour__cfour_mem_unit| =MB.

* Basis Sets
* other names to energy, optimize
* multiple jobs in file
* ref translation, more to come
* wrappers

.. warning:: Because p4c4 does not inspect the contents of the ``cfour {...}``
   block, once the user specifies a |PSIfour|-style molecule, the
   interface cannot judge whether a sandwich (drop the |PSIfour| molecule
   and use the cfour block as the entirety of the ``ZMAT``) or a standard
   (translate the |PSIfour| molecule and append additional input from the
   cfour block) mode is intended. The latter is what actually occurs. If
   there is both a |PSIfour| molecule and a molecule in the cfour block,
   ``ZMAT`` *will* end up with multiple molecules and ``*CFOUR(...)``
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

