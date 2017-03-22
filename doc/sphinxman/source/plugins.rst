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

.. _`sec:plugins`:

Plugins: Adding New Functionality to |PSIfour|
==============================================

Modular Approach to Development
-------------------------------

It is slightly cumbersome to the development cycle to recompile |PSIfour|
every time a change is made to the C++ code.
It's also daunting to new developers if they're required to
learn the structure of the source tree, executable initialization code,
and makefile systems in the existing code in order to add new features,
which was never a problem with Psi3 due to the independent
nature of the modules. To overcome these problems, |PSIfour| now has a
useful plugin feature. This allows codes to be developed as standalone
entities, which are compiled independently of the Psi source, but can
still link against Psi's vast library. The plugins can be loaded at
run-time from any location.

.. _`sec:newplugins`:

Creating a New Plugin
---------------------

|PSIfour| can create basic plugins for you and automatically tailor them
to your compilation environment. To create a basic plugin, run the
following while replacing ``myplugin`` with the name of your great code.
If the name you provide is not valid, |PSIfour| will complain.

   >>> psi4 --plugin-name myplugin

|PSIfour| will create a new directory with the name you specify for the
plugin. In this example, a directory named myplugin will be created.

All you need to do is cd into the directory, use |PSIfour| to generate
a Makefile, and type make. Then execute psi4 in the directory on the
default input file. ::

   >>> cd myplugin
   >>> `psi4 --plugin-compile`
   >>> make
   >>> psi4

|PSIfour| comes with a few templates that provide an excellent starting
point. These include code that demonstrates AO, MO, and SO integrals. Use
one of the following commands that meets your needs::

   >>> psi4 --plugin-name myplugin --plugin-template aointegrals
   >>> psi4 --plugin-name myplugin --plugin-template mointegrals
   >>> psi4 --plugin-name myplugin --plugin-template sointegrals
   >>> psi4 --plugin-name myplugin --plugin-template wavefunction
   >>> psi4 --plugin-name myplugin --plugin-template scf
   >>> psi4 --plugin-name myplugin --plugin-template dfmp2

..   >>> psi4 --plugin-name myplugin --plugin-template ambit

.. Several stable sample plugin directories are available to consult in the
.. :source:`plugins` directory. Other plugin directories can be used as models
.. but are in active development. For documentation on plugin modules, see
.. :ref:`Available Plugins <sec:availablePlugins>`.
.. 
.. * :source:`plugins/aointegrals/aointegrals.cc` 
..   An example that uses the LibMints library to generate and print AO basis (no symmetry) integrals.
.. 
.. * :source:`plugins/backtrans/backtrans.cc` 
..   A test of the one- and two-particle density matrix backtransformation code.
.. 
.. * :source:`plugins/mointegrals/mointegrals.cc` 
..   An example that uses the LibTrans library to generate and print MO basis integrals.
.. 
.. * :source:`plugins/mollerplesset2/mp2.cc` 
..   A plugin that uses LibTrans to generate open- and closed-shell MP2 energies.
.. 
.. * :source:`plugins/sointegrals/sointegrals.cc` 
..   An example that uses the LibMints library to generate and print SO basis (with symmetry) integrals.

.. _`sec:condaplugins`:

Creating a New Plugin Using a Conda Pre-compiled Binary
-------------------------------------------------------

|PSIfour| plugins can also be created using Conda for both |PSIfour|
binary and development environment. On Linux, one can use the ``gcc``
compiler installed alongside ``psi4`` itself in the Conda distribution
or environment (below, ``$PSI4CONDA``). ::

..     # prepare
..     >>> bash
..     >>> export PATH=$PSI4CONDA/bin:$PATH  # usually already done from psi4 installation
..     >>> cd "$(dirname $(which psi4))"/..  # move into distribution/environment directory, $PSI4CONDA
..     >>> conda install gcc                 # install compilers into expected place

    # check (yes, next line gives empty result. yes, LD_LIBRARY_PATH irrelevant)
    >>> echo $PYTHONHOME $PYTHONPATH $DYLD_LIBRARY_PATH $PSIDATADIR

    >>> which python psi4 gcc
    $PSI4CONDA/bin/python
    $PSI4CONDA/bin/psi4
    $PSI4CONDA/bin/gcc

    # create and compile plugin
    >>> psi4 --plugin-name testplugin     # generate new plugin
    -- Creating "testplugin" with "basic" template. -----------------
    ==> Created plugin files (in testplugin as basic):
      __init__.py, doc.rst, pymodule.py, plugin.cc, input.dat, CMakeLists.txt
    >>> cd testplugin                     # move into plugin directory
    >>> `psi4 --plugin-compile`           # use build info from parent psi4
    loading initial cache file $PSI4CONDA/share/cmake/psi4/psi4PluginCache.cmake
    -- The CXX compiler identification is GNU 5.2.0
    -- Check for working CXX compiler: $PSI4CONDA/bin/g++
    -- Check for working CXX compiler: $PSI4CONDA/bin/g++ -- works
    ...
    -- Generating done
    -- Build files have been written to: testplugin
    >>> make                              # compile the plugin to produce testplugin.so
    Scanning dependencies of target testplugin
    [ 50%] Building CXX object CMakeFiles/testplugin.dir/plugin.cc.o
    [100%] Linking CXX shared module testplugin.so
    [100%] Built target testplugin
    >>> psi4                              # run sample input.dat
    Attention! This SCF may be density-fitted.

Please note that the conda distribution must be in ``$PATH`` or the
conda enviroment must be activated before compilation and execution of
plugins created using this procedure.

Files in a Plugin Directory
---------------------------

In addition to the main ``myplugin.cc`` file, a fresh plugin directory contains the following files

* **CMakeLists.txt** |w---w| CMake file governing project *plugin*.
  The plugin source and CMakeLists.txt is independent of platform
  and |PSIfour| installation. You use CMake (version 3.1 or later)
  to generate a Makefile for the plugin by pointing it to a specific
  |PSIfour| installation. Run ``psi4 --plugin-compile`` to get a command
  to execute to generate the Makefile. What that command is doing is
  loading the compilers and options used to build the parent |PSIfour|
  (the ``-C psi4PluginCache`` part) which in turn can be overridden
  by passing ``-Doption=value`` commands to ``cmake`` *and* pointing
  toward a particular |PSIfour| (and probably pybind11) library to
  link against (the ``CMAKE_PREFIX_PATH`` part) *and* telling it to
  do an in-source build (the ``.`` part). Then just run ``make`` in
  your plugin directory. After any change to the plugin C++ code,
  ``make`` must be run in the plugin directory to recompile the
  ``myplugin.so`` executable, but recompiling the main |PSIfour| code
  is not necessary. Should you add additional (non-header) files to
  the plugin or need to link to additional external libraries, add that
  information here.

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

* **pymodule.py** |w---w| Python component of the plugin.
  By encoding the sequence of |PSIfour| module
  calls needed to run the plugin in the ``run_myplugin()`` function in this
  file, the plugin is hooked into the main |PSIfour| driver function
  :py:func:`~psi4.energy` and so can be accessed through
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

  .. literalinclude:: @SFNX_INCLUDE@psi4/share/psi4/plugin/__init__.py.template

* **doc.rst** |w---w| Documentation file. Place in this file any notes,
  equations, warnings to users, todo lists, *etc.*. Plain text is fine,
  though reStructuredText is the ultimate goal. Remove the ``.. comment``
  text and build Sphinx documentation for samples of linking keywords,
  sections, and math. This file is absorbed into the |PSIfour|
  documentation, along with any docstrings to Python functions, and the C++
  keywords block in the ``myplugin.cc`` file. See :ref:`sec:documentation`
  for building documentation.

.. and :ref:`Available Plugins <sec:availablePlugins>`
..  for this file's final destination.

Please note that pure virtual functions in a plugin may cause undefined symbols errors when
the plugin is loaded.

To create a purely Python plugin, create a new plugin directory, then
remove the ``Makefile`` and ``myplugin.cc`` files and
erase the shared object loading portion of ``__init__.py``. Create as many .py
files as necessary (registering each one in ``__init__.py``), use
``input.dat`` as a model for loading the plugin, no recompile ever
necessary.

