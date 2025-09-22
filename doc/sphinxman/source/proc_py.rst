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

.. _`sec:proc_py`:

Adding Methods to Driver
========================

``proc.py``
-----------

Methods that are computable by only one module should be added to the ``procedures`` dictionary in
:source:`psi4/driver/procrouting/proc_table.py`
that associates method names with functions
to run them located in :source:`psi4/driver/procrouting/proc.py`.

The function should start with a declaration, as below. ``methodname`` is
never seen by users, so it's good to be specific to method or module.
The function must always take as arguments ``(name, **kwargs)``. ::

    # energy method
    def run_methodname(name, **kwargs):

    # gradient method
    def run_methodname_gradient(name, **kwargs):

If the function needs to test the identity of ``name`` several times, it
can be convenient to predefine the lowercase version of the variable. The
case of all other py-side options (in kwargs) has already been handled by
:py:func:`~psi4.driver.energy()`, etc. in driver.py and need not be repeated here. ::

    # include if convenient
    lowername = name.lower()

    # never include
    kwargs = kwargs_lower(kwargs)

The function often needs to set options for the
c-side modules it calls. In order that the state of the options set by the
user remains when control is returned to the user, an
:py:class:`~psi4.driver.p4util.OptionsState` object is set up. See
:ref:`sec:handlingOptions_py` for details. *All* options set by the
function need to be included here, and *only* options set by the function
should be included. Most options should be associated with a particular
module, but a few (see below) are given without module. ::

    # include if any options set
    optstash = OptionsState(
        # these and other basis options should have no associated module
        ['BASIS'],
        ['DF_BASIS_SCF'],
        ['DF_BASIS_MP2'],
        ['PUREAM'],
        ['FREEZE_CORE'],
        # all others should have an associated module
        ['SCF', 'SCF_TYPE'],
        ['SCF', 'GUESS'],
        ['DFMP2', 'MP2_OS_SCALE'],
        )

If options need to be set, set them anywhere here. Options should be set
locally to a module, except for those without a module in
:py:class:`~psi4.driver.p4util.OptionsState`. ::

    # include if necessary as globals
    psi4.set_global_option('BASIS', guessbasis)
    psi4.set_global_option('DF_BASIS_SCF', guessbasisdf)

    # include if necessary as locals
    psi4.set_local_option('TRANSQT2', 'WFN', 'MP2')
    psi4.set_local_option('CCSORT', 'WFN', 'MP2')
    psi4.set_local_option('MP2', 'WFN', 'MP2')

If the regular scf module is to be run, run it through
``psi4.driver.procrouting.proc.scf_helper`` so that cast-up can be used. Also, add
the option to pass the reference wavefunction by pre-running scf,
then running the module with the ``ref_wfn`` kwarg.  Also, if the full
two-electron integrals are necessary for the post-scf, compute them if
only the df integrals were run previously. ::

    # Bypass the scf call if a reference wavefunction is given
    
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = scf_helper(name, **kwargs)
    
        # If the scf type is DF/CD, then the AO integrals were never written to disk
        if psi4.get_option('SCF', 'SCF_TYPE') in ['DF', 'CD']:
            psi4.MintsHelper(ref_wfn.basisset()).integrals()

Direct any post-scf modules to be run. ::

    # include if further post-scf modules are needed
    psi4.transqt2()
    psi4.ccsort()
    psi4.mp2()

If an :py:class:`~psi4.driver.p4util.OptionsState` object was set up, those options
need to be returned to the original user state with the following. ::

    # include if optstash = OptionsState(...) was set up previously
    optstash.restore()

Current best practice is to store as much as possible on the wavefunction, not in globals. The
driver should handle interactions with globals. When QCVariables are stored on the wavefunction in
the module, copy to globals with the below::

    # Shove variables into global space
    for k, v in dfmp2_wfn.variables().items():
        core.set_variable(k, v)

The function should return the wavefunction, except for rare cases like EFP where no wavefunction available.
For now, ``CURRENT ENERGY`` will be set by
:py:func:`~psi4.driver.energy`, etc. In future, this will be extracted from the wavefunction. ::

    # return highest or most prominent wavefunction (like dimer for SAPT)
    return fnocc_wfn


Managed Methods
---------------

There are several conditions when a method and derivative combination should be *managed*:

* when functionality overlaps between modules, a pattern is needed to
  access each route through the code;

* when functionality doesn't overlap completely, a pattern is needed to apportion defaulting among
  the modules, taking into account reference (RHF/UHF/ROHF), calc type (CONV/DF/CD), and possibly
  |globals__freeze_core| state (AE/FC).

* for higher-level derivatives, when, say, gradient functionality for mtd+ref+type+fcae doesn't
  exactly match energy functionality, a pattern is needed to decide analytic vs. finite difference.

* when default type is not available for a method (e.g., CCD governed by |globals__cc_type| that
  defaults to ``CONV`` but only ``DF`` and ``CD`` CCD is available), an informative error message is needed.

Managed methods handle these cases through the addition of a new
keyword |globals__qc_module| and a set of type keywords analogous to
|globals__mp2_type|: |globals__mp_type|,
|globals__ci_type|, |globals__cc_type|, which can have values ``CONV``,
``DF``, and ``CD``. These are all *global* keywords, as their values are
shared among modules rather than (or in addition to) being used internally
by the module). We're sticking with |globals__scf_type| and
|globals__mp2_type| defaulting to ``DF``, while most everything higher defaults
to ``CONV``. (Exceptions are MP2.5 and MP3 that default to ``DF``.)
In :source:`psi4/driver/procrouting/proc_table.py`, a managed method calls a
"select" function rather than a "run" function. ::

    procedures = {
        'energy': {
            'scf'           : run_scf,
            'mp3'           : select_mp3,
            'dct'           : run_dct,

Then in :source:`psi4/driver/procrouting/proc.py`, the select function runs through
reference, type, and possibly freeze_core to specify the proc
function to call for any able, non-default module (*e.g.*, ``mtd_type ==
'DETCI'`` ) or able, default module (*e.g.*, ``mtd_typd == ['', 'FNOCC']`` ).
Don't worry about 'else' statements as anything that falls through will be
caught and a readable error generated. ::

    def select_mp3(name, **kwargs):
        """Function selecting the algorithm for a MP3 energy call
        and directing to specified or best-performance default modules.

        """
        reference = psi4.get_option('SCF', 'REFERENCE')
        mtd_type = psi4.get_global_option('MP_TYPE')
        module = psi4.get_global_option('QC_MODULE')
        # Considering only [df]occ/fnocc/detci

        func = None
        if reference == 'RHF':
            if mtd_type == 'CONV':
                if module == 'DETCI':
                    func = run_detci
                elif module == 'FNOCC':
                    func = run_fnocc
                elif module in ['', 'OCC']:
                    func = run_occ
            elif mtd_type == 'DF':
                if module in ['', 'OCC']:
                    func = run_dfocc
            elif mtd_type == 'CD':
                if module in ['', 'OCC']:
                    func = run_dfocc
        elif reference == 'UHF':
            if mtd_type == 'CONV':
                if module in ['', 'OCC']:
                    func = run_occ
            elif mtd_type == 'DF':
                if module in ['', 'OCC']:
                    func = run_dfocc
            elif mtd_type == 'CD':
                if module in ['', 'OCC']:
                    func = run_dfocc
        elif reference == 'ROHF':
            if mtd_type == 'CONV':
                if module in ['DETCI']:
                    func = run_detci

        if func is None:
            raise ManagedMethodError(['select_mp3', name, 'MP_TYPE', mtd_type, reference, module])

        return func(name, **kwargs)

Naturally, in the run function, you must either use the type keyword for
type switching or translate it into whatever ``DO_CD``-like keywords your
module uses. At run time with a closed-shell molecule, ::

    energy('mp3')

will run OCC, while ::

    set qc_module fnocc
    energy('mp3')

will run FNOCC mp3.

A special case is DETCI that *can* run mp3, but oughtn't to be used for such. So above, ROHF CONV mp3 has no default, but can still access the detci code with ::

    set reference rohf
    set qc_module detci
    energy('mp3')

While the below gives an error ::

    set reference rohf
    energy('mp3')


Again, whenever a single method name needs to call multiple proc.py run
functions, it should be managed. In :ref:`table:managedmethods` "Y" means method available in
module, "D" means module is default for that method, "" mean method not
available.

