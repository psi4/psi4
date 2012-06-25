
.. include:: autodoc_abbr_options_c.rst

.. _`sec:handlingOptions_py`:

LibOptions: globals, locals, has_changed and all that
=====================================================

To simplify parsing of options and handling of defaults, the Options class
was created. It functions in the following way:

- Each module (or plugin) declares which options it will look for in the
  input: their name, type (string, int, double, array, etc.), and any
  default value they take.

- The input is parsed for these options, and defaults are assigned for
  those keywords not specified by the user.

- The c-side module or plugin can then query the Options object for the
  values associated with each keyword.

- The options will also be accessible py-side to the procedures that drive
  the modules. Array-type options are not available in python.

Declaring Options
-----------------

Each module needs to make itself known to the Options object, via a
read_options function. For plugins, this routine is provided by the user
in the plugin code. For native |PSIfour| modules, the entries need to
be appended to the read_options code in :source:`src/bin/psi4/read_options.cc`.
An example of such a routine is ::

    if (name == "MYMODULE"|| options.read_globals()) {
        /*- The amount of information printed
        to the output file -*/
        options.add_int("PRINT", 1);
        /*- Do save information to |mymodule__data_file| at the end of the computation? -*/
        options.add_bool("SAVE_INFO", true);
        /*- An array containing the number of doubly occupied orbitals per irrep 
        (in :ref:`Cotton order <table:irrepOrdering>`) -*/
        options.add("DOCC", new ArrayType());
        /*- The factor by which the harmonic vibrational frequencies are multiplied to
        obtain an approximation to the fundamental vibrational frequencies -*/
        options.add_double("FREQUENCY_SCALE_FACTOR", 1.0);
        /*- The filename to which data is dumped. !expert -*/
        options.add_str("DATA_FILE", "data.dat");
        /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
        options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    }

In the above example, the following options are declared (in order):

- An integer called ``PRINT`` with a default value of 1.
- A boolean called ``SAVE_INFO`` with a default of true.
- An array called ``DOCC``, no default is possible for this type.
- A double called ``FREQUENCY_SCALE_FACTOR`` with a default of 1.0.
- A string called ``DATA_FILE``, with a default of "data.dat" and any possible value.
- A string called ``AO_BASIS`` with a default of "NONE", and possible values of "NONE", "DISK", or "DIRECT".

The purpose of the "if" statement in the above read_options function is
the following. Suppose in an input file the user sets an option through
the construct ``set mymodule print 1`` or through a ``set mymodule {...}``
block. The first thing to happen is a call to read_options with name set
to "MYMODULE". (Note that all user input is converted to upper case.) This
call to read_options should tell the Options object only about those
options expected by the module called "mymodule"; this prevents overlap of
options between different modules.

Notice also that there's a special comment immediately before the
declaration of each keyword. You must provide these comments for any
options you add as they will be automatically inserted into the user
manual Providing a clear description will also help you to remember what
the keywords do and how they're used. The comments must live between the
special comment delimiters. For options that most users shouldn't need,
add an expert flag to the comment. This will place these options in a
separate section of the user manual. ::

   /*- comment -*/
   options.add_ ...
   /*- comment !expert -*/
   options.add_ ...

As is apparent from the examples above, comments can span multiple lines
(see ``PRINT``), can refer to other options (through hyperlinks; see
``SAVE_INFO``), can refer to sections of the manual (through hyperlinks;
see ``DOCC``), and can contain LaTeX notation (see ``AO_BASIS``). (To get
the LaTeX subscript command, use "@@" instead of "_".)

See `Best Practices <http://sirius.chem.vt.edu/trac/wiki/BestPractices#point1>`_ 
for guidelines on naming options.

What is *has_changed* ?
-----------------------

There are times when we need to know whether an option was provided by the
user or if the defaults are being used. For this reason, the Options
object stores a boolean *has_changed* value, in addition to the option
value itself.  A clarification of definition:

- [a] has_changed DOESN'T answer "Has option been changed by the user?"
- [b] has_changed DOESN'T answer "Is option now different from the default?"
- [c] has_changed DOES answer "Has option value been touched at all, by user or code?"

The above items notwithstanding, psi4 code should be written so that
*has_changed* DOES effectively mean, "Has option been changed by the
user?". The way to do this is to isolate and nullify any changes to
options made by the code, the difference between [a] and [c]. C-side,
there is no concern since options are read but essentially never
written-to within the modules.

Py-side is another matter since the driver's role is to take terse
instructions from the user and translate those into instructions to the
C++ modules, usually through manipulation of options.

.. comment * Reading options C-side
.. comment 
.. comment   This usually takes place within each module during options parsing,
.. comment   see section [] for details. No option is modified, so this process has
.. comment   no entanglement with the definition of has_changed.
.. comment 
.. comment * Setting options C-side
.. comment 
.. comment   This is very rarely done (what's going on in optking?). This
.. comment   document was written as background to the only case of note: PUREAM.
.. comment   This option has a clearly defined default, but it can't be set in
.. comment   read_options because its default depends on other options. This is a
.. comment   situation common to many options (and most all array options) and is
.. comment   generally handled within the module code and so is never available to
.. comment   the user. Analogously, PUREAM is handled in libmints but it is never
.. comment   reset.

In order to preserve effective definition [a], the strategy for each
python driver function is to query for the value of any option the
function may want to change and for the current has_changed status
(presumably reflecting whether the user has changed the value, as long as
no preceeding code has corrupted that definition). The python function
then makes its changes to the option and runs any c-side modules with
those changes. Finally, just before the function returns, the options are
reset to the user's value and has_changed status (which should now again
reflect only whether the user has changed the value).




.. comment options["AO_BASIS"].has_changed()
.. comment will return false if the default value is being used, and true if the user specified this keyword in the input.


.. warning:: |globals__puream| is an exception in that its value and
   ``has_changed()`` value only reflect what the user has explicitly set.
   This keyword should not be queried to find out the current
   |globals__puream| state for the active basis; use instead,
   ``PsiMod.MintsHelper().basisset().has_puream()``.




Reading Options in Module
-------------------------




Handling Options in Driver
--------------------------

This section is about the scopes of options and how best to handle them in
the python driver. There are four groups of commands available.
Options from the c-side Options object are accessible in the Python driver through four sets of commands.

- get 

  - :py:func:`~PsiMod.get_global_option()`
  - :py:func:`~PsiMod.get_local_option()`
  - :py:func:`~PsiMod.get_option()`

- set 

  - :py:func:`~PsiMod.set_global_option()` 
  - :py:func:`~PsiMod.set_local_option()`

- has_changed 

  - :py:func:`~PsiMod.has_global_option_changed()`
  - :py:func:`~PsiMod.has_local_option_changed()`
  - :py:func:`~PsiMod.has_option_changed()`

- revoke_changed 

  - :py:func:`~PsiMod.revoke_global_option_changed()`
  - :py:func:`~PsiMod.revoke_local_option_changed()`

There's a pattern here. Setting something, either a value (set) or a
negative changed status (revoke_changed), can only be done for a specific
scope, either global or local to the specified module. Querying, either a
value (get) or a changed status (has_changed), can be done in the global
scope, in a specified local scope, or in the context of "What will the
specified module use?".

.. note:: "Global" in the sense of the discussion has *nothing*
   to do with the globals section at the top of read_options.cc .  That
   section is just a convenient place for options and associated values
   that are used by most, if not all, modules.

.. comment Those options could be distributed out to
   all the modules below and the globals section dissolved with no change
   to psi's operation.
   :source:`src/bin/psi4/read_options.cc`. That section is just a

There are two primary purposes for interacting with options in the python driver.

- **Preserving User Options** (Enforcing definition [a] of has_changed)

  The first, less-interesting, use of retrieving user option values has
  been to preserve them so that they may be restored at the end after the
  procedure itself has clobbered them. By decoupling global_option and
  local_option commands, this can now be performed neatly by saving at the
  beginning the global and local values and the global and local
  has_changed values, then restoring them at the end.  Below is an example
  of this procedure; don't actually do this. ::

    g_user_scftype = PsiMod.get_global_option('SCF_TYPE')
    l_user_scftype_scf = PsiMod.get_local_option('SCF', 'SCF_TYPE')
    bg_user_scftype = PsiMod.has_global_option_changed('SCF_TYPE')
    bl_user_scftype_scf = PsiMod.has_local_option_changed('SCF', 'SCF_TYPE')

    g_user_wfn = PsiMod.get_global_option('WFN')
    l_user_wfn = PsiMod.get_local_option('MP2', 'WFN')
    bg_user_wfn = PsiMod.has_global_option_changed('WFN')
    bl_user_wfn = PsiMod.has_local_option_changed('MP2', 'WFN')

    # body of function
    # scf_type and wfn are freely changed, LOCALLY
    # PsiMod.scf() and PsiMod.mp2() are run

    PsiMod.set_global_option('SCF_TYPE', g_user_scftype)
    if not bg_user_scftype:
        PsiMod.revoke_global_option_changed('SCF_TYPE')
    PsiMod.set_local_option('SCF', 'SCF_TYPE', l_user_scftype_scf)
    if not bl_user_scftype_scf:
        PsiMod.revoke_local_option_changed('SCF', 'SCF_TYPE')

    PsiMod.set_global_option('WFN', g_user_wfn)
    if not bg_user_wfn:
        PsiMod.revoke_global_option_changed('WFN')
    PsiMod.set_local_option('MP2', 'WFN', l_user_wfn_scf)
    if not bl_user_wfn_scf:
        PsiMod.revoke_local_option_changed('MP2', 'WFN')

  Instead of cluttering the driver with the above boilerplate, use an
  :py:class:`~optproc.OptionsState` object that stores values and
  has_changed values for each keyword and module pair given to it as
  arguments. At the end of the python function, these stored user settings
  can be restored. ::

    optstash = OptionsState(
        ['SCF', 'SCF_TYPE'],
        ['MP2', 'WFN'],
        ['DF_BASIS_SCF'])

    # body of function
    # scf_type and wfn are freely changed, LOCALLY
    # puream and df_basis_scf are freely changed, GLOBALLY
    # PsiMod.scf() and PsiMod.mp2() are run

    optstash.restore()


  .. note:: Some options (BASIS, BASIS-like, and PUREAM) should always
     be used globally (no module argument) with the OptionsState objects.
     Similarly, within the body of the function, they should always be
     queried and set globally.

- **Setting-Up Calculations**

  The other types of options calls in python driver functions are (a)
  those to query what option value an upcoming c++ module is going to use
  (determined by user and defaults) and (b) those to set options to govern
  the course of a procedure. Finding out the intended option value for a
  molecule should employ the :py:func:`~PsiMod.get_option` command 
  (and :py:func:`~PsiMod.has_option_changed` for has_changed), which
  (newly) requires a module for scope. *(Previously, this command used the
  "active module", which isn't well-defined in the context of the python
  driver, and consequently, the command gave variable results, depending
  on whether a get_local/set_local command had been previously executed to
  define the active module.)* ::

    if (PsiMod.get_option('SCF', 'REFERENCE') == 'RHF'):
        PsiMod.set_local_option('SCF', 'REFERENCE', 'RKS')

  Setting of options in python should use the
  :py:func:`~PsiMod.set_local_option` command. Using the local, rather
  than global, scope will ensure that the newly set option will be used by
  the module. Otherwise, if the python procedure set in the global scope
  and the user had happened to set that option in local scope, the local
  user option will take precedence against the programmer's intent.
  *(Anyone who has heard advice to "query local, set global" should forget
  that and follow the new scheme outlined here.)*


PsiMod Options Commands
-----------------------

.. function:: PsiMod.get_global_option(keyword)

   Given a string of *keyword* name, returns the value associated with the
   keyword from the global options. Returns error if keyword is not
   recognized.

.. function:: PsiMod.get_local_option(module, keyword)

   Given a string of *keyword* name and a particular *module*, returns the
   value associated with the keyword in the module options scope. Returns
   error if keyword is not recognized for the module.

.. function:: PsiMod.get_option(module, keyword)

   Given a string of *keyword* name and a particular *module*, returns the
   local value associated with the keyword if it's been set, else the global
   value if it's been set, else the local default value. Returns error if
   keyword is not recognized globally or if keyword is not recognized for the
   module.

.. function:: PsiMod.set_global_option(keyword, value)

   Sets *value* to option *keyword* for all modules.
    
.. function:: PsiMod.set_local_option(module, keyword, value)

   Sets *value* to option *keyword* scoped only to specific *module*.

.. function:: PsiMod.has_global_option_changed(keyword)

   Returns boolean for whether the *keyword* has been touched in the global
   scope, by either user or code. Notwithstanding, code is written such that
   in practice, this returns whether the option has been touched in the
   global scope by the user.

.. function:: PsiMod.has_local_option_changed(module, keyword)

   Returns boolean for whether *keyword* has been touched in the scope of
   the specified *module*, by either user or code. Notwithstanding, code is
   written such that in practice, this returns whether the option has been
   touched in the module scope by the user.

.. function:: PsiMod.has_option_changed(module, keyword)

   Returns boolean for whether *keyword* has been touched either locally
   to specified *module* or globally, by either user or code.
   Notwithstanding, code is written such that in practice, this returns
   whether the option has been touched by the user.

.. function:: PsiMod.revoke_global_option_changed(keyword)

   Given a string of *keyword* name, sets the has_changed attribute in the
   global options scope to false. Used in python driver when a function
   sets the value of an option. Before the function exits, this command is
   called on the option so that has_changed reflects whether the user (not
   the program) has touched the option.

.. function:: PsiMod.revoke_local_option_changed(module, keyword)

   Given a string of *keyword* name and a particular *module*, sets the
   has_changed attribute in the module options scope to false. Used in
   python driver when a function sets the value of an option. Before the
   function exits, this command is called on the option so that
   has_changed reflects whether the user (not the program) has touched the
   option.

