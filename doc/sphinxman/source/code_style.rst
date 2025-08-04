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

.. _`sec:code_style`:

Code style conventions
======================

It is important to keep a consistent formatting of the C++ and Python code
to avoid hard-to-read diffs and merge conflicts.
`clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_ and `yapf <https://github.com/google/yapf>`_ can be used to format C++ and Python code,
respectively, according to a predefined style file.
|PSIfour| provides the :source:`.clang-format` and :source:`.style.yapf` files in the root
folder of the project.
It is **recommended** that modifications and/or new files checked into the
repository are formatted according to these style files using ``clang-format``
and ``yapf``. It is then helpful if these tools be part of your development toolchain.
Once ``clang-format`` and ``yapf`` are installed, there are three ways in which
formatting of the code can be accomplished, in decreasing order of automation:

1. By integrating the formatters into your editor.
2. By installing Git hooks to run the formatters when committing.
3. By running the formatters manually on the modified files.

.. _`faq:editorcodestyle`:

How to impose code style through your editor
--------------------------------------------

Both ``clang-format`` and ``yapf`` can be integrated into widely used editors.
The `Neoformat <https://github.com/sbdchd/neoformat>`_ plugin can be configured
to format files when saving them to disk.

.. _`faq:githookscodestyle`:

How to impose code style through Git hooks
------------------------------------------

Git hooks are scripts that are run before or after certain Git events.
In this particular case, we want to make sure that all files that have been
added to the staging area with ``git add`` are formatted according to the style
*before* they committing them with ``git commit``.
The hook to be modified is then the *pre-commit* hook.
|PSIfour| uses the `pre-commit <https://pre-commit.com/>`_ framework, with configuration file :source:`.pre-commit-config.yaml`.
To take advantage of pre-commit hooks, you will need to install the ``pre-commit`` utility:

::
  pip install pre-commit

or using Conda:

::
  conda install pre_commit -c conda-forge

Finally, you need to install the actual hooks:

::
  pre-commit install

Pre-commit hooks will be run on every ``git commit``, but the ``--no-verify``
option can be used to skip their execution.

Hooks are powerful, but integrating the formatter into your editor will prove
to be better. Hooks need to be installed anew for every fresh clone of the
repository you are working on.

.. _`faq:manualcodestyle`:

How to run code-style tools `clang-format` and `yapf` manually
--------------------------------------------------------------

The least recommended approach to formatting your code is to run manually the
formatters. The following commands will format only the files that have been
modified:

::
  clang-format -style=file -i `git diff --relative --name-only HEAD -- *.cc *.h`
  yapf -i `git diff --relative --name-only HEAD -- *.py`

How and when to *not* apply code styling to your contributions
--------------------------------------------------------------

TODO

