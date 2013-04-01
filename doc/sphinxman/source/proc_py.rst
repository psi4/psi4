
.. include:: autodoc_abbr_options_c.rst

.. _`sec:proc_py`:

Adding Methods to Driver
========================

This is concerned at present with normal methods added first to the
procedures table in driver.py that associates method names with functions
to run them located in proc.py .

The function should start with a declaration, as below. ``methodname`` is
never seen by users, so it's good to be specific; if there's lots of
modules that can run mp2, call methodname modulenamemethodname, perhaps.
The function must always take as arguments ``(name, **kwargs)``. ::

    # energy method
    def run_methodname(name, **kwargs):

    # gradient method
    def run_methodname_gradient(name, **kwargs):

If the function needs to test the identity of ``name`` several times, it
can be convenient to predefine the lowercase version of the variable. The
case of all other py-side options (in kwargs) has already been handled by
:py:func:`~driver.energy()`, etc. in driver.py and need not be repeated here. ::

    # include if convenient
    lowername = name.lower()

    # never include
    kwargs = kwargs_lower(kwargs)

It's often necessary to The function often needs to set options for the
c-side modules it calls. In order that the state of the options set by the
user remains when control is returned to the user, an
:py:class:`~optproc.OptionsState` object is set up. See
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
:py:class:`~optproc.OptionsState`. ::

    # include if necessary as globals
    PsiMod.set_global_option('BASIS', guessbasis)
    PsiMod.set_global_option('DF_BASIS_SCF', guessbasisdf)

    # include if necessary as locals
    PsiMod.set_local_option('TRANSQT2', 'WFN', 'MP2')
    PsiMod.set_local_option('CCSORT', 'WFN', 'MP2')
    PsiMod.set_local_option('MP2', 'WFN', 'MP2')

If the regular scf module is to be run, run it through
:py:func:`~proc.scf_helper` so that cast-up can be used. Also, add the
option to bypass it by pre-running scf, then running the module with this
``bypass_scf`` kwarg.  Also, if the full two-electron integrals are
necessary for the post-scf, compute them if only the df integrals were run
previously. ::

    # include if scf module is to be run

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
         scf_helper(name, **kwargs)
 
         # include if TEI are needed beyond scf

         # If the scf type is DF, then the AO integrals were never generated
         if PsiMod.get_option('SCF', 'SCF_TYPE') == 'DF':
             mints = PsiMod.MintsHelper()
             mints.integrals()
 
Direct any post-scf modules to be run. ::

    # include if further post-scf modules are needed
    PsiMod.transqt2()
    PsiMod.ccsort()
    PsiMod.mp2()

If an :py:class:`~optproc.OptionsState` object was set up, those options
need to be returned to the original user state with the following. ::

    # include if optstash = OptionsState( was set up previously
    optstash.restore()

No function should return anything. ``CURRENT ENERGY`` will be set by
:py:func:`~driver.energy`, etc. ::

    # never include
    return returnvalue

