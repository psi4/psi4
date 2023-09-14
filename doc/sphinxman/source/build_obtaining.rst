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


===================
Obtaining |PSIfour|
===================

.. warning:: As of v1.8, primary binary distribution has moved from
   the psi4 channel to the conda-forge channel. Neither install docs
   for users nor compile docs for developers have been updated yet to
   reflect new patterns. Please consult :psicode:`psicode downloads
   page <installs/latest/>` for the latest guides.

.. _`faq:obtainpsi4`:

How to obtain Psi4: start with find-the-code quiz, end in ``{top-level-psi4-dir}``
----------------------------------------------------------------------------------

A better decision tree is available at :psicode:`installs/latest`,
though the below remains valid.

Take a :ref:`quiz <faq:quiz>` to find the best version of the codebase for
your needs, be it binary, tarball, or version-controlled repository. Or,
select outright among:

#. :ref:`faq:binary`
#. :ref:`faq:binarypackage`
#. :ref:`faq:clonepsi4public`
#. :ref:`faq:forkpsi4public` (only path to develop |PSIfour|)
#. :ref:`faq:tarballpsi4public`


.. _`faq:quiz`:

Find-the-code Quiz
------------------

A better decision tree is available at :psicode:`installs/latest`,
though the below remains valid.

* I just want to run the code. I may tweak the Python, but I'm not
  developing anything to contribute back to the code base.

  * Provided I still get good, threaded BLAS/LAPACK, OpenMP parallelism,
    and optimization for a variety of processor architectures, I'm willing to forgo
    architecture tuning wizardry to avoid compiling it myself.

    * I'm on Linux or Mac (Intel or Silicon chips) or Windows (native or WSL/Ubuntu Bash Shell).

      * I'm familiar with conda and want to manage |PSIfour| as an
        ordinary conda package. |w---w| :ref:`Goto Binary-Package
        <faq:binarypackage>`

      * I just want a |PSIfour| installer. |w---w| :ref:`Goto
        Binary-Installer <faq:binary>`

  * I want to compile it myself to eke out best performance on my
    computer. I accept responsibility for navigating compiler, threading,
    and BLAS/LAPACK compatibility

    * I'm willing to have minimal dealings with git (e.g., commands ``git
      clone`` and ``git pull``) in return for easy access in future to new
      features and bug fixes. |w---w| :ref:`Goto Clone-from-GitHub
      <faq:clonepsi4public>`

    * I don't want to deal with this newfangled git, just give me a
      tarball of the source code |w---w| :ref:`Goto Tarball-from-GitHub
      <faq:tarballpsi4public>`

* I want to run *and* develop in |PSIfour|.

  * In keeping with the open-source philosophy, I don't mind my code being
    as public as Psi4 itself during the development process. |w---w|
    :ref:`Goto Fork-from-GitHub <faq:forkpsi4public>`

  * I want to develop *using* |PSIfour| infrastructure and libraries, not
    *on* them; I think a plugin might do.

    * I've got a |PSIfour| compilation. Use it, then consult :ref:`plugins
      <sec:newplugins>`

    * I'd rather not compile |PSIfour| or I don't have compilers |w---w|
      :ref:`Goto Binary-Package <faq:binarypackage>` then consult
      :ref:`plugins through conda <sec:condaplugins>`

* I really like parentheses and/or DBOC, so I want Psi3. |w---w|
  Psi3 is available from `sourceforge <https://sourceforge.net/projects/psicode/files/psi/3.4.0/>`_, but you're on your own.

.. comment * I am a core |PSIfour| developer, yet I'm still taking this quiz.
.. comment 
.. comment   * I have minions whose Psi4 development work I want to supervise through this repository instance. Preferably, [Goto Fork-from-GitHub](#forkpsi4public); otherwise [Goto Fork-from-GitHub-Private](#forkpsi4private)
.. comment 
.. comment   * Just give me a repository to commit to directly. Preferably, [Goto Clone-from-GitHub](#clonepsi4public); otherwise [Goto Clone-from-GitHub-Private](#clonepsi4private)


.. _`faq:binary`:

Binary Installer
----------------

* **Get Initially**

  Just go to http://www.psicode.org/downloads.html, select "Installer",
  "Stable Release", and your choice of architecture and Python version,
  and follow the instructions there.

* **Build**

  Not applicable as binary is pre-built.

* **Get Updates** :ref:`directions <faq:updatepsi4>`

  .. code-block:: bash

     >>> conda update psi4

* **Contribute Back**

  Not applicable as not under git control.


.. _`faq:binarypackage`:

Conda Binary Package
--------------------

* **Get Initially**

  The pre-compiled conda packages at https://anaconda.org/conda-forge/psi4
  can be installed into an existing Anaconda or Miniconda distribution
  according to :ref:`directions <faq:psi4pkg>`. Locally, install into
  a conda environment as below.

  .. code-block:: bash

     >>> conda create -n p4env psi4 -c conda-forge/label/libint_dev -c conda-forge
     >>> conda activate p4env

  .. code-block:: bash

     >>> # nightly build (Linux and Windows only)
     >>> conda create -n p4env psi4/label/dev::psi4 -c conda-forge/label/libint_dev -c conda-forge
     >>> # release
     >>> conda create -n p4env                 psi4 -c conda-forge/label/libint_dev -c conda-forge


* **Build**

  Not applicable as binary is pre-built.

* **Get Updates** :ref:`directions <faq:updatepsi4>`

  .. code-block:: bash

     >>> conda update psi4 -c conda-forge

* **Contribute Back**

  Not applicable as not under git control.


.. _`faq:clonepsi4public`:

Clone from GitHub Repository
----------------------------

* **Get Initially**

  The |PSIfour| repository at https://github.com/psi4/psi4 works like
  `every other GitHub repo
  <https://help.github.com/articles/which-remote-url-should-i-use/>`_.
  Locally, clone as below.

  .. code-block:: bash

     # use https or ssh
     >>> git clone https://github.com/psi4/psi4.git
     >>> cd psi4
     # this is your {top-level-psi4-dir}

* **Build** :ref:`directions <faq:buildquick>`

* **Get Updates**

  .. code-block:: bash

     # on branch master
     >>> git pull origin master

* **Contribute Back**

  Contributions cannot be made directly to the main repository. :ref:`Fork
  instead <faq:forkpsi4public>`.

  To convert clone to fork, go to https://github.com/psi4/psi4, and
  hit the `Fork <https://help.github.com/articles/fork-a-repo/>`_
  button to store a |PSIfour| repository in your GitHub account.

  .. code-block:: bash

     >>> git remote rename origin upstream
     >>> git remote add origin https://github.com/johndoe/psi4.git


.. _`faq:forkpsi4public`:

Fork from GitHub Repository
---------------------------

* **Get Initially**

  Go to https://github.com/psi4/psi4, and hit the `Fork
  <https://help.github.com/articles/fork-a-repo/>`_ button to store a
  |PSIfour| repository in your GitHub account. Locally, proceed to clone:

  .. code-block:: bash

     # replace johndoe
     # use https or ssh
     >>> git clone https://github.com/johndoe/psi4.git
     >>> cd psi4
     # this is your {top-level-psi4-dir}

  `Set up a connection
  <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`_
  between your forked repository and the parent repository.

  .. code-block:: bash

     >>> git remote add upstream https://github.com/psi4/psi4.git

* **Build** :ref:`directions <faq:buildquick>`

* **Get Updates**

  Locally, `update your fork
  <https://help.github.com/articles/syncing-a-fork/>`_ from the parent
  repository and store on GitHub at your fork.

  .. code-block:: bash

     # on branch working_branch
     >>> git pull --rebase upstream master
     >>> git push origin working_branch

  Remember: Working in the master branch of a fork is considered bad practice.

* **Contribute Back**

  |PSIfour| contributions process :ref:`here <faq:githubworkflow>` and
  :source:`here <.github/CONTRIBUTING.md>`.
  Consider `preparing your contribution in a branch
  <http://blog.jasonmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/>`_
  then issue a `GitHub pull request
  <https://help.github.com/articles/creating-a-pull-request/>`_.


.. _`faq:tarballpsi4public`:

Tarball from GitHub Repository
------------------------------

* **Get Initially**

  Discouraged! From the |PSIfour| repository at https://github.com/psi4/psi4, hit the
  "Clone or download" then "Download ZIP" button. Locally, unpack as
  below.

  .. code-block:: bash

     >>> unzip psi4-master.zip
     >>> cd psi4-master
     # this is your {top-level-psi4-dir}

* **Build** :ref:`directions <faq:buildquick>`

* **Get Updates**

  Download new tarball and rebuild.

* **Contribute Back**

  Not applicable as source not under git control.

.. _`faq:githubworkflow`:

What is the suggested GitHub workflow
-------------------------------------

.. image:: /prflow.001.jpeg
.. image:: /prflow.002.jpeg
.. image:: /prflow.003.jpeg
