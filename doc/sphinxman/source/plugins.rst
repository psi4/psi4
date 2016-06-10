
.. include:: autodoc_abbr_options_c.rst

.. _`sec:plugins`:

Plugins: Adding New Functionality to |PSIfour|
==============================================

Modular Approach to Development
-------------------------------

The redesign of |PSIfour| into a single-executable changed the way that
code development is done. The standalone nature of modules in previous
versions of Psi made their development very easy, as new functionality
could be implemented almost as a standalone executable, which could easily
be ported into the Psi code when completed. The new design specifies that
these modules are now libraries, not separate executables, which are
linked into the main Psi executable. The single-executable design is
conducive to code reuse, as it allows common tasks to be implemented as a
class instead of a module; the functionality can then be easily obtained
throughout the code by creating objects as required. Examples of this are
the LibMints class, which provides similar functionality to the old cints
module, and LibTrans, which replaces the old transqt code. When codes are
wrapped in a library, they should be placed into :source:`src/lib`, and
codes that resemble modules belong in :source:`src/bin`.

The single-executable design leads to a somewhat cumbersome development
cycle, since every time a change is made, one must compile the code,
archive it into a library, and then re-link the code into the main
executable. It's also daunting to new developers if they're required to
learn the structure of the source tree, executable initialization code,
and makefile systems in the existing code in order to add new features,
which was never a problem with previous versions due to the independent
nature of the modules. To overcome these problems, |PSIfour| now has a
useful plugin feature. This allows codes to be developed as standalone
entities, which are compiled independently of the Psi source, but can
still link against Psi's vast library. The plugins can be loaded at
run-time from any location.

Creating a New Plugin
---------------------

|PSIfour| can create basic plugins for you and automatically tailor them
to your compilation environment. To create a basic plugin, run the
following while replacing ``myplugin`` with the name of your great code.
If the name you provide is not valid, |PSIfour| will complain.

   >>> psi4 --new-plugin myplugin

|PSIfour| will create a new directory with the name you specify for the
plugin. In this example, a directory named myplugin will be created.

All you need to do is cd into the directory and type make. Then execute
psi4 in the directory on the default input file.

|PSIfour| comes with a few templates that provide an excellent starting
point. These include code that demonstrates AO, MO, and SO integrals. Use
one of the following commands that meets your needs::

   >>> psi4 --new-plugin myplugin +aointegrals
   >>> psi4 --new-plugin myplugin +mointegrals
   >>> psi4 --new-plugin myplugin +sointegrals
   >>> psi4 --new-plugin myplugin +wavefunction
   >>> psi4 --new-plugin myplugin +scf
   >>> psi4 --new-plugin myplugin +ambit


Several stable sample plugin directories are available to consult in the
:source:`plugins` directory. Other plugin directories can be used as models
but are in active development. For documentation on plugin modules, see
:ref:`Available Plugins <sec:availablePlugins>`.

* :source:`plugins/aointegrals/aointegrals.cc` 
  An example that uses the LibMints library to generate and print AO basis (no symmetry) integrals.

* :source:`plugins/backtrans/backtrans.cc` 
  A test of the one- and two-particle density matrix backtransformation code.

* :source:`plugins/mointegrals/mointegrals.cc` 
  An example that uses the LibTrans library to generate and print MO basis integrals.

* :source:`plugins/mollerplesset2/mp2.cc` 
  A plugin that uses LibTrans to generate open- and closed-shell MP2 energies.

* :source:`plugins/sointegrals/sointegrals.cc` 
  An example that uses the LibMints library to generate and print SO basis (with symmetry) integrals.

.. _`sec:condaplugins`:

Creating a New Plugin Using a Conda Pre-compiled Binary
-------------------------------------------------------

|PSIfour| plugins can also be created using Conda for both |PSIfour| binary and 
development environment.
To compile a plugin it is necessary to have a compiler (*e.g.*, ``gcc``) and blas libraries
(*e.g.*, ``openblas``) installed in the Conda environment used to run |PSIfour|.
It is recommended to create a new Conda environment with packages `gcc` and `openblas` installed. ::

   >>> conda create -n p4plugenv psi4     # makes environment named p4plugenv with psi4 binary installed
   >>> source activate p4plugenv          # activate the env so its contents are first in your PATH
   >>> cd "$(dirname $(which psi4))"/..   # move into env directory
   >>> #cd $CONDA/envs/p4plugenv          # same effect as line above where $CONDA is path to Miniconda/Anaconda
   >>> conda install gcc                  # place compilers into expected place
   # Linux
   >>> conda install openblas
   # Mac
   >>> conda install boost=1.57

Once these packages are installed, plugins can be created and compiled. ::

   >>> source activate p4plugenv          # important: activate conda environment
   >>> psi4 --new-plugin testplugin       # generate new plugin
   >>> cd testplugin                      # move into plugin directory
   >>> make                               # compile the plugin to produce testplugin.so
   >>> psi4                               # run sample input.dat 

Please note that the conda enviroment must be activated before compilation and execution of
plugins created using this procedure.

Files in a Plugin Directory
---------------------------

In addition to the main ``myplugin.cc`` file, a fresh plugin directory contains the following files

* **Makefile** |w---w| Makefile for the directory. As long as you are the
  only user of the plugin, this should not need editing. After any change to
  the plugin C++ code, ``make`` must be run in the plugin directory to
  recompile the ``myplugin.so`` executable, but recompiling the main
  |PSIfour| code is not necessary.

* **input.dat** |w---w| Sample input file for the plugin.
  Since the ``__init__.py`` file makes the plugin directory look like a
  Python module, the plugin can be treated as such in an input file. The
  location of the plugin directory must be included in :envvar:`PYTHONPATH`,
  either externally in the calling shell or defined in the input file. 
  This is usually done by manipulating :envvar:`PSIPATH`. Then,
  the plugin can be loaded as ``import myplugin`` and executed as
  ``energy('myplugin')``. Any other Python functions are also available from
  the input file, *e.g.* ``myplugin.testfunction()``, note the namespace
  protection.

* **pymodule.py** |w---w| Python component of the plugin. The procedure
  for calling plugin code shown in ``input.dat`` sounds very simple, but it
  can be made simpler still. By encoding the sequence of |PSIfour| module
  calls needed to run the plugin in the ``run_myplugin()`` function in this
  file, the plugin is hooked into the main |PSIfour| driver function
  :py:func:`~driver.energy` and so can be accessed through
  ``energy('myplugin')`` in an input file. Any other Python functions can
  also be placed in this file.

* **__init__.py** |w---w| Init script for the plugin (in the sense that
  the whole plugin directory is a Python module). This file generally won't
  need editing unless additional Python files are added to the plugin
  directory (add additional lines to the ``# Load Python modules`` section)
  or the plugin depends on .so codes in other plugin directories (add
  additional plugin_load lines relative to the current plugin directory to
  the ``# Load C++ plugin`` section).

.. comment  as modeled in :source:`tests/plugin_libcim/__init__.py`).

  .. literalinclude:: @SFNX_INCLUDE@share/plugin/__init__.py.template

* **doc.rst** |w---w| Documentation file. Place in this file any notes,
  equations, warnings to users, todo lists, *etc.*. Plain text is fine,
  though reStructuredText is the ultimate goal. Remove the ``.. comment``
  text and build Sphinx documentation for samples of linking keywords,
  sections, and math. This file is absorbed into the |PSIfour|
  documentation, along with any docstrings to Python functions, and the C++
  keywords block in the ``myplugin.cc`` file. See :ref:`sec:documentation`
  for building documentation and :ref:`Available Plugins <sec:availablePlugins>` 
  for this file's final destination.

Please note that pure virtual functions in a plugin may cause undefined symbols errors when
the plugin is loaded.

To create a purely Python plugin, create a new plugin directory, then
remove the ``Makefile`` and ``myplugin.cc`` files and
erase the shared object loading portion of ``__init__.py``. Create as many .py
files as necessary (registering each one in ``__init__.py``), use
``input.dat`` as a model for loading the plugin, no recompile ever
necessary.

