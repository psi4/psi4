
.. include:: autodoc_abbr_options_c.rst

.. _`sec:bestPractices_py`:

Best Practices for Python Functions
===================================

- Thy python functions shall always have final argument \*\*kwargs, that they may take in and pass on keywords meant for other functions. Yea, even the run_mcscf(), and run_ccsd() -type functions that have no use for kwargs. The exceptions are python functions that are only helpers called by a driver function.

- Python functions should read the kwargs dictionary and (possibly) add to it. Functions should not pop or remove keywords from kwargs, even those keywords meaningful only to itself. This will ensure that the complete kwargs is available for pickling and sow/reap procedures. The exception is the molecule argument, which is read by the first function that gets ahold of it. This first function activates the molecule and pops it out of kwargs, effectively setting molecule for all subsequent functions. The code below should suffice. ::

    # Make sure the molecule the user provided is the active one
    if 'molecule' in kwargs:
        activate(kwargs['molecule'])
        del kwargs['molecule']
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()

- Preferrably, the python function signature (for functions intended to be called in input files) is ``function(name, **kwargs)``. For functions that have other positional keywords, please bundle them into kwargs at earliest convenience (see :ref:`sec:db()` argument db_name for example).

- After the docstring, the first two lines of your function should be the ones below. The first provides a case insensitive handle to the name argument value. The second converts all the kwargs dictionary keys to lowercase versions of themselves, so that input files can be case insensitive. ::

    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

- Case sensitivity for kwargs dictionary values still needs to be handled. The first line below shows how to convert argument values to lowercase for matching. When not matching a whole value such that regular expressions are needed, the second line below performs a case insensitive match. ::

    if (kwargs['db_mode'].lower() == 'continuous'):
    if re.match(r'^sapt', name, flags=re.IGNORECASE):

- Match boolean keywords (db_cp in the example below) with expressions like the following, which allow case insensitive yes/true/on/1/no/false/off/0 user input. If your argument's value is a derivative level, similarly, use input.der0th, input.der1st, and input.der2nd. ::

    if input.yes.match(str(db_cp)):
    elif input.no.match(str(db_cp)):

-   For keywords that might be used in other functions as well as your own, prepend the argument name with a short representation of your function name. For example, there are keywords cp_func, db_func, and opt_func to request what python function, if not energy(), is called by cp(), database(), and optimize().

- Upon checking in a new python file, edit the file ``psi4/doc/userman/source/index.rst`` and follow the instructions therein that your file may be autodocumented here.

- Write docstrings! For a major function intended for use in input files, start with the skeleton docstring in ``psi4/lib/python/example_docstring`` and replace anything that looks like ``<this>``. For a behind-the-scenes function or if you don't want the bother of dealing with `reStructuredText <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_, just write an ordinary docstring. It will get slurped into the documentation in plain text.

- Your python function should follow `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_ conventions (without the line-length restriction). I'm aiming for files to pass the line below, unless for good reason. The second line is for database Python files.

    >>> pep8.py -r --ignore=E501 pythonfile.py
    >>> pep8.py -r --ignore=E501,E221,E222,E241,E201,E202 databasefile.py

- Your python function should not prevent any test case (``make tests``, NOT ``make longtests``) from passing. A test case(s) should be written and checked in for any major python function, so that others do not break your code. If most of your work was on the python (as opposed to c++) side, the test case prefix pywrap\_ is suggested.

- Be sure to set any new PSI variables through lines like those below. Especially if the function returns an energy, set the 'current energy' variable. This last is needed to communicate with the optimizer. ::

    PsiMod.set_variable('MP2.5 CORRELATION ENERGY', ce_mp25)
    PsiMod.set_variable('MP2.5 TOTAL ENERGY', e_mp25)
    PsiMod.set_variable('CURRENT ENERGY', e_mp25)

- Once your python function is fairly stable on its own, it's potential for interoperability with energy()/opt()/cp()/db()/cbs()/etc. should be evaluated. If it makes physical sense that it should work, you should strive to make that interoperability a reality. Some steps:

    - If any interoperability is possible, define an argument xx_func, where xx is a short name for your function. Add near the top of your function code like the below (less the final two lines). The net result of this code is that if the user specifies no \*_func arguments, then energy() gets called. If the user defines xx_func, then its value gets called. If the user defines func, then its value gets reassigned to xx_func, func itself is deleted, and xx_func() gets called. Whatever is getting called is stored in func within the function. ::

        # Establish function to call
        if not('xx_func' in kwargs):
            if ('func' in kwargs):
                kwargs['xx_func'] = kwargs['func']
                del kwargs['func']
            else:
                kwargs['xx_func'] = energy
        func = kwargs['xx_func']
        if not func:
            raise ValidationError('Function \'%s\' does not exist to be called by wrapper counterpoise_correct.' % (func.__name__))
        if (func is db):
            raise ValidationError('Wrapper xx is unhappy to be calling function \'%s\'.' % (func.__name__))

    - If specific interoperabilities are known, code them in. For example, if xx shouldn't call db, add the last two lines above to the xx function. If db shouldn't call xx, add the following two lines below to the db function. ::

        if (func is xx):
            raise ValidationError('Wrapper database is unhappy to be calling function \'%s\'.' % (func.__name__))

    - Create a multipart test case that runs some intercalls between your function and others (akin to :srcsample:`pywrap_all)`. In trials, permute the order of calls a few times to expose any calls that don't clean up after themselves and need further attention.

    - When all is validated, add your findings to the great :ref:`table:intercalls` table in the documentation.



