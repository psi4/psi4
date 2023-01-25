.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2023 The Psi4 Developers.
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

.. .. _`sec:addAddOns`:

Git, Versioning
===============

.. _`faq:versionbump`:

How to bump a version
---------------------

0. **ACT** to check everything in

1. **OBSERVE** current versioning state

  * Be on master of (i) a direct clone or (ii) clone-of-fork with master
    up-to-date with upstream (including tags!!!) and with upstream as
    remote.

  * https://github.com/psi4/psi4/releases says ``v1.1a1`` & ``007a9b6``

  ::

    >>> git tag
    v1.0
    v1.1a1

    >>> cat psi4/metadata.py
    __version__ = '1.1a1'
    __version_long = '1.1a1+007a9b6'
    __version_upcoming_annotated_v_tag = '1.1a2'

    >>> git describe --abbrev=7 --long --always HEAD
    v1.1a1-417-gcbee32b

    >>> git describe --abbrev=7 --long --dirty
    v1.1a1-417-gcbee32b

    >>> ./psi4/versioner.py
    Defining development snapshot version: 1.1a2.dev417+cbee32b (computed)
    1.1a2.dev417 {master} cbee32b 1.0.0.999   1.0 <-- 1.1a2.dev417+cbee32b

    >>> git diff

  * Observe that current latest tag matches metadata scipt and git
    describe, that GitHub releases matches metadata script, that upcoming in
    metadata script matches current versioner version.

  * Note that current tag is ``v1.1a1``. Decide on imminent tag, say ``v1.1rc1``.

2. **ACT** to bump tag in code

  * Edit current & prospective tag in :source:`psi4/metadata.py`. Use your
    decided-upon tag ``v1.1rc1`` and a speculative next tag, say ``v1.1rc2``,
    and use 7 "z"s for the part you can't predict.

  ::

    >>> vi psi4/metadata.py

    >>> git diff
    diff --git a/psi4/metadata.py b/psi4/metadata.py
    index 5d87b55..6cbc05e 100644
    --- a/psi4/metadata.py
    +++ b/psi4/metadata.py
    @@ -1,6 +1,6 @@
    -__version__ = '1.1a1'
    -__version_long = '1.1a1+007a9b6'
    -__version_upcoming_annotated_v_tag = '1.1a2'
    +__version__ = '1.1rc1'
    +__version_long = '1.1rc1+zzzzzzz'
    +__version_upcoming_annotated_v_tag = '1.1rc2'

    >>> git add psi4/metadata.py

    >>> git commit -m "v1.1rc1"

3. **OBSERVE** undefined version state

  ::

    >>> git describe --abbrev=7 --long --always HEAD
    v1.1a1-418-g6100822

    >>>  git describe --abbrev=7 --long --dirty
    v1.1a1-418-g6100822

    >>>  psi4/versioner.py
    Undefining version for irreconcilable tags: 1.1a1 (computed) vs 1.1rc1 (recorded)
    undefined {master} 6100822 1.0.0.999   1.0 <-- undefined+6100822

  * Note 7-char git hash for the new commit, here "6100822".

4. **ACT** to bump tag in git, then bump git tag in code.

  * Use the decided-upon tag ``v1.1rc1`` and the observed hash "6100822" to
    mint a new *annotated* tag, minding that "v"s are present here.

  * Use the observed hash to edit :source:`psi4/metadata.py` and commit immediately.

  ::

    >>> git tag -a v1.1rc1 6100822 -m "v1.1rc1"

    >>> vi psi4/metadata.py
    >>> git diff
    diff --git a/psi4/metadata.py b/psi4/metadata.py
    index 6cbc05e..fdc202e 100644
    --- a/psi4/metadata.py
    +++ b/psi4/metadata.py
    @@ -1,5 +1,5 @@
     __version__ = '1.1rc1'
    -__version_long = '1.1rc1+zzzzzzz'
    +__version_long = '1.1rc1+6100822'
     __version_upcoming_annotated_v_tag = '1.1rc2'

    >>> psi4/versioner.py
    Amazing, this can't actually happen that git hash stored at git commit.
    >>> git add psi4/metadata.py
    >>> git commit -m "Records tag for v1.1rc1"

5. **OBSERVE** current versioning state

  * Nothing to make note of, this is just a snapshot.

  ::

    >>> psi4/versioner.py
    Defining development snapshot version: 1.1rc2.dev1+4e0596e (computed)
    1.1rc2.dev1 {master} 4e0596e 1.0.0.999   1.0 <-- 1.1rc2.dev1+4e0596e

    >>> git describe --abbrev=7 --long --always HEAD
    v1.1rc1-1-g4e0596e

    >>> git describe --abbrev=7 --long --dirty
    v1.1rc1-1-g4e0596e

    >>> git tag
    v1.0
    v1.1a1
    v1.1rc1

    >>> cat psi4/metadata.py
    __version__ = '1.1rc1'
    __version_long = '1.1rc1+6100822'
    __version_upcoming_annotated_v_tag = '1.1rc2'

    >>> cat metadata.out.py | head -8
    __version__ = '1.1rc2.dev1'
    __version_branch_name = 'master'
    __version_cmake = '1.0.0.999'
    __version_is_clean = 'True'
    __version_last_release = '1.0'
    __version_long = '1.1rc2.dev1+4e0596e'
    __version_prerelease = 'False'
    __version_release = 'False'

    >>> git log --oneline
    4e0596e Records tag for v1.1rc1
    6100822 v1.1rc1
    cbee32b Fixes pcmsolver/scf for py3. Moves source for libefp upstream.

6. **ACT** to inform remote of bump

  * Temporarily disengage "Include administrators" on protected master branch.

  ::

    >>> git push origin master
    >>> git push origin v1.1rc1

  * Now https://github.com/psi4/psi4/releases says ``v1.1rc1`` & ``6100822``


.. _`faq:remotetag`:

How to create and remove an annotated Git tag on a remote
---------------------------------------------------------

|PSIfour| versioning only works with *annotated* tags, not *lightweight*
tags as are created with the `GitHub interface
<https://github.com/psi4/psi4/releases/new>`_

* Create *annotated* tag::

    >>> git tag -a v1.1a1 <git hash if not current> -m "v1.1a1"
    >>> git push origin v1.1a1

* Delete tag::

    >>> git tag -d v1.1a1
    >>> git push origin :refs/tags/v1.1a1

* Pull tags::

    >>> git fetch <remote> 'refs/tags/*:refs/tags/*'


.. _`faq:psi4version`:

What Psi4 version is running
----------------------------

* Psithon / from the executable::

    >>> psi4 --version
    1.1rc2.dev17

* PsiAPI / from the library::

    >>> python -c "import psi4; print(psi4.__version__)"
    1.1rc2.dev17

* Output file header gives info like the ``print_header()`` below.

* Function ``print_header()`` returns a summary of citation, version, and
  git information about |PSIfour|. Function ``version_formatter()`` can
  return version and git information in any desired format string. ::

    >>> import psi4
    >>> psi4.print_header()

        -----------------------------------------------------------------------
              Psi4: An Open-Source Ab Initio Electronic Structure Package
                                   Psi4 1.1rc2.dev17

                             Git: Rev {condadoc} c852257 dirty


        R. M. Parrish, L. A. Burns, D. G. A. Smith, A. C. Simmonett,
        A. E. DePrince III, E. G. Hohenstein, U. Bozkaya, A. Yu. Sokolov,
        R. Di Remigio, R. M. Richard, J. F. Gonthier, A. M. James,
        H. R. McAlexander, A. Kumar, M. Saitow, X. Wang, B. P. Pritchard,
        P. Verma, H. F. Schaefer III, K. Patkowski, R. A. King, E. F. Valeev,
        F. A. Evangelista, J. M. Turney, T. D. Crawford, and C. D. Sherrill,
        submitted.

        -----------------------------------------------------------------------


        Psi4 started on: Friday, 28 April 2017 07:31PM

        Process ID:  95107
        PSIDATADIR: /Users/johndoe/psi4/objdir8/stage/usr/local/psi4/share/psi4
        Memory:     500.0 MiB
        Threads:    1

    >>> psi4.version_formatter()
    '1.1rc2.dev17'
    >>> psi4.version_formatter('all')
    '1.1rc2.dev17 {condadoc} c852257 1.0.0.999 dirty  1.0 <-- 1.1rc2.dev17+c852257'
    >>> psi4.version_formatter("""{{{branch}}} {versionlong}""")
    '{condadoc} 1.1rc2.dev17+c852257'


.. _`faq:grepascii`:

How to locate non-ascii characters in the codebase
--------------------------------------------------

Neither the Python interpreter nor Sphinx like non-ASCII characters one
bit, though the errors may be intermittant. Output files are usually ok,
so Jerome can live, for now. To aid in tracking down offenders, here's
the ``vi`` and ``grep`` search strings. In the docs, you want to use
the substitutions in :source:`doc/sphinxman/source/abbr_accents.rst`
instead of the actual characters. ::

    # vim
    :/[^\x00-\x7F]

    # bash
    grep -r --color='auto' -P -n "[^\x00-\x7F]" psi4/
    
.. _`faq:undefversion`:

How to fix "Psi4 undefined" version
-----------------------------------

When in a git repo, the versioner uses ``git describe`` and psi4/metadata.py
to compute the version. If you don't have all the latest tags, this mechanism
can't work. To solve, pull tags and remake. ::

    # upstream in `git remote -v` points to github.com/psi4/psi4.git
    >>> git fetch upstream 'refs/tags/*:refs/tags/*'
    >>> make
    # version healed
    
.. _`faq:cannotimportcoretlpd`:

How to fix "cannot import name 'core' from {top-level-psi4-dir}
---------------------------------------------------------------

First, what's happening? ``sys.path`` (where modules can be imported from in python) starts with ``''``.  If you `export PYTHONPATH={objdir}/stage/{prefix}/lib/{pymod_lib_dir}:$PYTHONPATH` to make PsiAPI easy, that inserts starting in pos'n 1 (0-indexed), so ``''`` still at the head of ``sys.path``. Now, if you try to run a psiapi/python file from ``{top-level-psi4-dir}`` that contains ``import psi4``, it will find the source tree ``psi4/__init__.py`` and fail because there's no ``core.so`` around. That is, it's finding what looks to be the psi4 module dir structure ``.`` when the one it wants is what you inserted into PYTHONPATH at pos'n 1.

The way around this is to move the python file you're running to any other directory. Or, within the file, do ``sys.path.insert(0, {objdir}/stage/{prefix}/lib/{pymod_lib_dir}``.


.. _`faq:findmissingoutputref`:

How to find tests without output.ref
------------------------------------

Ideally, each new test or much-altered test should add its own
``output.ref``. When that doesn't happen, this command helps. ::

    find tests/ -mindepth 1 -maxdepth 1 -type d '!' -exec test -e "{}/output.ref" ";" -print

.. _`faq:githubcodereview`:

How to do GitHub issue management and code review
-------------------------------------------------

a) Anyone, core-dev or not, is encouraged to review PRs. It's actually good practice for interacting with other open-source projects, where you don't have the advantage of knowing or working with the contributors. Before venturing into projects on GitHub where you don't know the maintainers, it doesn't hurt to read https://snarky.ca/setting-expectations-for-open-source-participation/ .

b) Psi4 is a learning tool for all involved, so partial reviews in areas of confidence and questions and comments on PRs in general are encouraged.

c) Approving before CI completes is fine, though it can be mildly personally embarrassing when CI catches something you didn't.

d) All main branches (master and `1.N.x` maintenance) are protected by GitHub, including administrators, so even with write access, no one can accidentally push (master) or rewrite the history (master and maintenance).

e) PR owners who also have maintainer status can merge their PRs as GitHub enforces three external reviews.

f) Unless there's been a lot of discussion on core-dev about merge order, generally the 3rd positive reviewer merges the PR. Also fine to add review and leave merge for later.

g) Presently only Travis-CI is set up as a required-to-merge service. Incomplete Azure won't block merging, but we do usually let it complete before merging unless it's a trivial PR.

h) We don't enforce branches to be up to date before merging since that'd be a lot of extra CI time and coordination when merging several PRs in a day. So, if a PR hasn't been updated in a while, and a reviewer is nervous about PR interference, fine to ask submitter to rebase. For this reason, we try to merge newer contributors first so the rebase falls on more experienced contributors.

i) Ideally a PR consists of atomic, compilable commits. When the PR instead is many successive small changes toward a single goal, consider squashing the PR. For core-dev's PRs, there's implicit permission to squash (unless otherwise noted in PR intro), whereas for new contributors, we often let the commits be messy.

j) When discussion on issue has overcome the original problem and settled on needing long-term work, fine to move the long-term item to Wish List and close issue.

