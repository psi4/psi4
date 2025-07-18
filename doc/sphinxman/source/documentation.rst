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

.. _`sec:documentation`:

Updating the |PSIfour| Users' and Programmers' Manual
=====================================================

|PSIfours| documentation is generated by `Sphinx <https://www.sphinx-doc.org/>`_
and lives in :source:`doc/sphinxman`. It is available online at
`<http://psicode.org/psi4manual/master/index.html>`_
for the latest development branch.

Installing Sphinx
^^^^^^^^^^^^^^^^^

Installing Sphinx is only necessary to build the documentation 
yourself, locally. The docs are served from
from psicode, so most users and developers won't need Sphinx
installed. Nevertheless, installation is easy.
Math is rendered through MathJax, so LaTeX and dvipng are no longer needed.
The sphinx executable should be in your path at CMake configure time for
documentation targets to be available.

* Binary: ``conda install sphinx``
* Binary: ``conda env create -f devtools/conda-envs/docs-cf.yaml``
* Binary: ``pip install -U Sphinx``
* Source: https://pypi.org/project/Sphinx/

* Check:

   >>> which sphinx-build
   //anaconda/bin/sphinx-build
   >>> sphinx-build --version  # needs >= 3.5
   Sphinx (sphinx-build) 3.5.3
   >>> cmake
   ...
    -- Documentation targets available: sphinxman (html), sphinxmini (quick html), sphinxpdf (LaTeX --> PDF)
   ...

Documentation Structure
^^^^^^^^^^^^^^^^^^^^^^^

Sphinx has nice capabilities for extracting docstrings from python files,
presenting both auto-generated and narrative documentation in the same
format, hyperlinking within and to external websites, and generating
documentation in different formats from the same source. |PSIfours|
documentation is a unified document covering information for both users
and programmers in separate sections. From the top-level object directory,
build the following target (note that a working version of the |PSIfour|
executable is a requirement for building the
documentation). Only GNU Makefiles, not Ninja, works for the docs:

.. code-block:: console

    >>> make sphinxman
    # -OR-
    >>> cmake --build . --target sphinxman

This will build a full set of documentation in the ``html`` directory that can be viewed offline through any browser. ::

    doc/sphinxman/html/index.html
    
Much of the documentation is auto-generated from the source. At present,
this covers:

* Physical Constants: :source:`psi4/include/psi4/physconst.h`
* Python Driver: docstrings from \*.py files in :source:`psi4/driver`
* Databases: docstrings from \*.py files in :source:`psi4/share/psi4/databases`
* Basis Sets: \*.gbs files in :source:`psi4/share/psi4/basis`
* C++ Keywords: :source:`psi4/src/read_options.cc`
* Sample Inputs: input.dat files in :source:`samples`
* PSI Variables: variables and associated modules extracted from code and comments in the Python and C++ source
  * Modules scraped are the sections of :source:`psi4/src/read_options.cc`
  * Variables should be all-caps, except where representing substitutions, e.g., ``ROOT n -> ROOT m`` and double-quote ``"`` delimited, even in Python
  * Scraper looks for ``Process::environment.globals``, ``set_array_variable``, ``variables_``, etc. lines and comments in the C++ code
  * C-side, the module for the variable is determined by the directory where it's found.
  * Scraper looks for ``set_variable`` together with ``# P::e MODULE`` lines and comments in the Python code
  * Py-side, the module for the variable is specified by ``MODULE`` in the comment
  * When a variable is set by code in either language, e.g., ``variables_[varname.str()]`` rather than plain string, ``variables_["FCI TOTAL ENERGY"]``, add a plain string line as a single-line comment, so the scraper can find it.
  * Add new places to scrape for variables to :source:`doc/sphinxman/document_psivariables.pl`
  * For now, we're scraping both global and Wfn variables
  * All of these show up in referenceable appendices like ``apdx:detci_psivar``
* Plugins: ``doc.rst`` text, \*.py modules, and C++ keywords in ``psi4/tests/plugin_*`` plugin directories (disabled at the moment)
* PSI Files: scratch file names and numbers in :source:`psi4/include/psi4/psifiles.h`

Some documentation is even extracted from |PSIfour| objects at runtime.

* psi4: docstrings for the C++ submodule ``psi4.core`` and the Python submodule ``psi4.driver`` that comprise |PSIfour|. C++ docstrings from "core" and "export" files in :source:`psi4/src/`, and Py docstrings from :source:`psi4/driver/`.
* DFT: functional availability and characteristics as encoded in :source:`psi4/driver/procrouting/dft`
* BasisFamily: fitting basis sets for each orbital basis as encoded in :source:`psi4/driver/qcdb/basislistdunning.py` and :source:`psi4/driver/qcdb/basislistother.py`

Building all the documentation takes ~10 minutes. There is now good
dependency structure built into the :source:`doc/sphinxman/CMakeLists.txt`
, so very long builds should be infrequent (unless you're touching
:source:`psi4/src/read_options.cc` or the driver. Note that not all dependencies are
encoded (PSI variables, for instance, depend on every .cc file in the
source tree), so for a definitive doc build, remove (in the object
directory) ``doc/sphinxman`` and start from scratch.

Even ~10 minutes of build time can be annoying when developing
documentation and testing ``rst`` files. In that situation, use the target
below which builds only the written docs (not autodocs) in
``psi4/doc/sphinxman/source`` quickly, though with a lot of warnings for
unresolved links::

    >>> make sphinxmini

reStructuredText
^^^^^^^^^^^^^^^^

Sphinx files are written in reStructuredText (\*.rst). In the html
documentation, source code is available from the sidebar. Here are a
few resources on Sphinx formatting.

* `reStructuredText <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_
* `links example <https://github.com/sphinx-doc/sphinx/issues/5208#issuecomment-736259355>`_
* `rendered test document <http://docutils.sourceforge.net/test/functional/expected/standalone_rst_html4css1.html>`_
  *vs.* `source test document <http://svn.python.org/projects/external/docutils-0.5/docs/user/rst/demo.txt>`_
* `Sphinx Docs <https://www.sphinx-doc.org/en/master/contents.html>`_

.. * `Another reStructuredText <http://people.ee.ethz.ch/~creller/web/tricks/reST.html>`_
.. * `LaTeX that Sphinx can handle <ftp://ftp.ams.org/ams/doc/amsmath/short-math-guide.pdf>`_

Math in the Codebase
^^^^^^^^^^^^^^^^^^^^

It is often useful to have mathematical expressions in docstrings or
comments in the code that are auto-documented into the manual. Such
locations include the ``#! comment`` comments at the top of test case
input files, the ``/*- comment -*/`` comments in
:source:`psi4/src/read_options.cc`, and the ``r""" comment """``
docstrings in python modules. (That ``r"""`` makes the string read
literally, so your LaTeX symbols aren't confused with escape characters.)
For the two former, math has traditionally
been written in LaTeX (with the special substitution ``@@`` for
subscripting underscore). The autodoc script has been trained to convert
inline LaTeX math to reST math, provided the expression within dollar
signs is offset from other text. That is, expressions of the form
:regexp:`^ $latex math$[., ]$` (pseudo-regex) are good, while ``H$_2$O`` and LaTeX tables
are not translated correctly. Python docstrings are absorbed as-is, so
please use reST math formatting (essentially ``$latex math$`` :math:`\Rightarrow`
``:math:`latex math```).
Starting around |PSIfour| 1.1, MathJax is used for in-browser LaTeX
rendering in place of offline PNG generation of math images. Check the
online rendering, as occasionally there will be errors even when the LaTeX
looked sound.

The Map of the Sphinx
^^^^^^^^^^^^^^^^^^^^^

* Adding a new Appendix or First-TOC-Level page

  Create your reST file and fill it with information. Add the name of your
  file to :source:`doc/sphinxman/source/appendices.rst` for an appendix or
  to :source:`doc/sphinxman/source/index.rst` for a first-TOC-level.
  Finally, add your file to the ``STATICDOC`` variable in
  :source:`doc/sphinxman/CMakeLists.txt`. Sphinx will now build with your
  new page.

* Adding a new module to "Theoretical Methods"

  Copy the file of a well-established module, like
  :source:`doc/sphinxman/source/sapt.rst`. Change the title, author, sec
  label, ref, and source labels at the top of the file to point instead to
  your code. Edit :source:`doc/sphinxman/source/methods.rst` to add the
  name of your file so that it will appear in the TOC tree. Add your file
  to the ``STATICDOC`` variable in
  :source:`doc/sphinxman/CMakeLists.txt`. Sphinx will now build with your new
  file.  Follow the models in existing methods pages to write your
  documentation. If you don't get all the keyword links, bibliography
  links, sample inputs, math, tables, etc. working in Sphinx, don't worry
  about it. A genie will probably come through and tidy up all your
  source.

