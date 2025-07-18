.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2025 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This file is part of Psi4.
.. #
.. # Psi4 is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU Lesser General Public License as published by
.. # the Free Software Foundation, version 3.
.. #
.. # Psi4 is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU Lesser General Public License for more details.
.. #
.. # You should have received a copy of the GNU Lesser General Public License along
.. # with Psi4; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #

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
read_options section. For plugins, this routine is provided by the user
in the plugin code. For native |PSIfour| modules, the entries need to
be appended to the read_options code in :source:`psi4/src/read_options.cc`.
An example of such a routine is

.. code-block:: cpp

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
        options.add_str_i("DATA_FILE", "data.dat");
        /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
        options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    }

In the above example, the following options are declared (in order):

- An integer called ``PRINT`` with a default value of 1.
- A boolean called ``SAVE_INFO`` with a default of true.
- An array called ``DOCC``, no default is possible for this type.
- A double called ``FREQUENCY_SCALE_FACTOR`` with a default of 1.0.
- A case-sensitive string called ``DATA_FILE``, with a default of "data.dat" and any possible value.
- A string called ``AO_BASIS`` with a default of "NONE", and possible values of "NONE", "DISK", or "DIRECT".

The purpose of the "if" statement in the above read_options function is
the following. Suppose in an input file the user sets an option through
the construct ``set mymodule print 1`` or through a ``set mymodule {...}``
block. The first thing to happen is a call to read_options with name set
to "MYMODULE". (Note that all user input is converted to upper case unless a
``add_str_i`` which should be used sparingly for files.) This
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
separate section of the user manual.

.. code-block:: cpp

   /*- comment -*/
   options.add_ ...
   /*- comment !expert -*/
   options.add_ ...

As is apparent from the examples above, comments can span multiple lines
(see ``PRINT``), can refer to other options (through hyperlinks; see
``SAVE_INFO``), can refer to sections of the manual (through hyperlinks;
see ``DOCC``), and can contain LaTeX notation (see ``AO_BASIS``). (To get
the LaTeX subscript command, use "@@" instead of "_".)

See :ref:`faq:readoptions`
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
there is no concern since options are essentially read-only
within the modules.

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
no preceding code has corrupted that definition). The python function
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
   ``psi4.MintsHelper().basisset().has_puream()``.




Reading Options in Module
-------------------------




Handling Options in Driver
--------------------------

This section is about the scopes of options and how best to handle them in
the python driver. There are four groups of commands available.
Options from the c-side Options object are accessible in the Python driver through four sets of commands.

- get 

  - :py:func:`psi4.core.get_global_option()`
  - :py:func:`psi4.core.get_local_option()`
  - :py:func:`psi4.core.get_option()`

- set 

  - :py:func:`psi4.core.set_global_option()`
  - :py:func:`psi4.core.set_local_option()`

- has_changed 

  - :py:func:`psi4.core.has_global_option_changed()`
  - :py:func:`psi4.core.has_local_option_changed()`
  - :py:func:`psi4.core.has_option_changed()`

- revoke_changed 

  - :py:func:`psi4.core.revoke_global_option_changed()`
  - :py:func:`psi4.core.revoke_local_option_changed()`

There's a pattern here. Setting something, either a value (set) or a
negative changed status (revoke_changed), can only be done for a specific
scope, either global or local to the specified module. Querying, either a
value (get) or a changed status (has_changed), can be done in the global
scope, in a specified local scope, or in the context of "What will the
specified module use?".

.. note:: "Global" in the sense of the discussion has *nothing*
   to do with the globals section at the top of :source:`psi4/src/read_options.cc`. That
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

    from psi4 import core

    g_user_scftype = core.get_global_option('SCF_TYPE')
    l_user_scftype_scf = core.get_local_option('SCF', 'SCF_TYPE')
    bg_user_scftype = core.has_global_option_changed('SCF_TYPE')
    bl_user_scftype_scf = core.has_local_option_changed('SCF', 'SCF_TYPE')

    g_user_wfn = core.get_global_option('WFN')
    l_user_wfn = core.get_local_option('MP2', 'WFN')
    bg_user_wfn = core.has_global_option_changed('WFN')
    bl_user_wfn = core.has_local_option_changed('MP2', 'WFN')

    # body of function
    # scf_type and wfn are freely changed, LOCALLY
    # core.scf() and core.mp2() are run

    core.set_global_option('SCF_TYPE', g_user_scftype)
    if not bg_user_scftype:
        core.revoke_global_option_changed('SCF_TYPE')
    core.set_local_option('SCF', 'SCF_TYPE', l_user_scftype_scf)
    if not bl_user_scftype_scf:
        core.revoke_local_option_changed('SCF', 'SCF_TYPE')

    core.set_global_option('WFN', g_user_wfn)
    if not bg_user_wfn:
        core.revoke_global_option_changed('WFN')
    core.set_local_option('MP2', 'WFN', l_user_wfn_scf)
    if not bl_user_wfn_scf:
        core.revoke_local_option_changed('MP2', 'WFN')

  Instead of cluttering the driver with the above boilerplate, use an
  :py:class:`~psi4.driver.p4util.OptionsState` object that stores values and
  has_changed values for each keyword and module pair given to it as
  arguments. At the end of the python function, these stored user settings
  can be restored. ::

    from psi4.driver import p4util

    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'],
        ['MP2', 'WFN'],
        ['DF_BASIS_SCF'])

    # body of function
    # scf_type and wfn are freely changed, LOCALLY
    # puream and df_basis_scf are freely changed, GLOBALLY
    # core.scf() and core.mp2() are run

    optstash.restore()


  .. note:: Some options (BASIS, BASIS-like, and PUREAM) should always
     be used globally (no module argument) with the OptionsState objects.
     Similarly, within the body of the function, they should always be
     queried and set globally. Same for FREEZE_CORE.

- **Setting-Up Calculations**

  The other types of options calls in python driver functions are (a)
  those to query what option value an upcoming c++ module is going to use
  (determined by user and defaults) and (b) those to set options to govern
  the course of a procedure. Finding out the intended option value for a
  molecule should employ the :py:func:`~psi4.core.get_option` command
  (and :py:func:`~psi4.core.has_option_changed` for has_changed), which
  requires a module for scope. (Programmer-supplied scope is needed Py-side,
  whereas C-side, commands use the "active module".) ::

    if (psi4.core.get_option("SCF", "REFERENCE") == "RHF"):
        psi4.core.set_local_option("SCF", "REFERENCE", "RKS")

  Setting of options in python should use the
  :py:func:`~psi4.core.set_local_option` command. Using the local, rather
  than global, scope will ensure that the newly set option will be used by
  the module. Otherwise, if the python procedure set in the global scope
  and the user had happened to set that option in local scope, the local
  user option will take precedence against the programmer's intent.

