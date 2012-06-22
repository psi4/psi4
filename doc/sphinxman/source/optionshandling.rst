
.. include:: autodoc_abbr_options_c.rst

.. _`sec:handlingOptions_py`:

LibOptions: globals, locals, has_changed and all that
=====================================================



has_changed
-----------

There are times when we need to know if an option was provided by the user
or if the defaults are being used. For this reason, the Options object
stores a boolean *has_changed* value, in addition to the option value itself. 
A clarification of definition:

a. has_changed DOESN'T mean "Has option been changed by the user?"
b. has_changed DOESN'T mean "Is option now different from the default?"
c. has_changed DOES mean "Has option value been touched at all, by user or code?"

The above items notwithstanding, psi4 code should be written so that
*has_changed* DOES effectively mean, "Has option been changed by the
user?". The way to do this is to isolate and nullify any changes to
options made by the code, the difference between *a.* and *c.*. C-side,
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

In order to preserve effective definition *a.*, the strategy for each
python driver function is to query for the value of any option the
function may want to change and for the current has_changed status
(presumably reflecting whether the user has changed the value, as long as
no preceeding code has corrupted that definition). The python function
then makes its changes to the option and runs any c-side modules with
those changes. Finally, just before the function returns, the options are
reset to the user's value and has_changed status (which should now again
reflect only whether the user has changed the value).




options["AO_BASIS"].has_changed()
will return false if the default value is being used, and true if the user specified this keyword in the input.


.. warning:: |globals__puream| is an exception in that its value and
   ``has_changed()`` value only reflect what the user has explicitly set.
   This keyword should not be queried to find out the current
   |globals__puream| state for the active basis; use instead,
   ``PsiMod.MintsHelper().basisset().has_puream()``.






Options in the Python Driver
----------------------------

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

- **Preserving User Options** (Enforcing definition (a) of has_changed)

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
        ['MP2', 'WFN'])

    # body of function
    # scf_type and wfn are freely changed, LOCALLY
    # PsiMod.scf() and PsiMod.mp2() are run

    optstash.restore()


  .. warning:: The OptionsState procedure is not currently compatible with BASIS-type options.

- **Setting-Up Calculations**

  The other types of options calls in python driver functions are (a)
  those to query what option value an upcoming c++ module is going to use
  (determined by user and defaults) and (b) those to set options to govern
  the course of a procedure. Finding out the intended option value for a
  molecule should employ the :py:func:`~PsiMod.get_option` command, which
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

  .. warning:: Stick with setting BASIS-type options in global scope for now.


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

