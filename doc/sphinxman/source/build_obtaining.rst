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


===================
Obtaining |PSIfour|
===================

.. _`faq:obtainpsi4`:

How to obtain Psi4: start with find-the-code quiz, end in ``{top-level-psi4-dir}``
----------------------------------------------------------------------------------

Take a :ref:`quiz <faq:quiz>` to find the best version of the codebase for
your needs, be it binary, tarball, or version-controlled repository. Or,
select outright among:

#. :ref:`faq:binary`
#. :ref:`faq:clonepsi4public` (*read-only* unless core developer)
#. :ref:`faq:forkpsi4public`
#. :ref:`faq:tarballpsi4public`
#. :ref:`faq:psi3sourceforge`


.. _`faq:quiz`:

Find-the-code Quiz
------------------

* I just want to run the code. I may tweak the Python, but I'm not
  developing anything to contribute back to the code base.

  * Provided I still get good, threaded BLAS/LAPACK, I'm willing to
    sacrifice processor architecture tuning to avoid compiling it myself.

    * I'm on Linux or Mac or Windows with Ubuntu Bash Shell.

      * I'm familiar with conda and want to manage |PSIfour| as an
        ordinary conda package. |w---w| :ref:`Goto Binary-Package
        <faq:binarypackage>`

      * I just want a |PSIfour| installer. |w---w| :ref:`Goto
        Binary-Installer <faq:binary>`

  * I want to compile it myself for best performance on my computer.

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

.. comment  * I have scientific competitors, and I don't want to get scooped. [Goto Fork-from-GitHub-Private](#forkpsi4private)

  * I want to develop *using* |PSIfour| infrastructure and libraries, not
    *on* them; I think a plugin might do.

    * I've got a |PSIfour| compilation. Use it, then consult :ref:`plugins
      <sec:newplugins>`

    * I'd rather not compile |PSIfour| or I don't have compilers |w---w|
      :ref:`Goto Binary-Package <faq:binarypackage>` then consult
      :ref:`plugins through conda <sec:condaplugins>`

* I really like parentheses and/or DBOC, so I want Psi3. |w---w|
  :ref:`Goto Psi3-from-SourceForge <faq:psi3sourceforge>`

.. comment * I am a core |PSIfour| developer, yet I'm still taking this quiz.
.. comment 
.. comment   * I have minions whose Psi4 development work I want to supervise through this repository instance. Preferably, [Goto Fork-from-GitHub](#forkpsi4public); otherwise [Goto Fork-from-GitHub-Private](#forkpsi4private)
.. comment 
.. comment   * Just give me a repository to commit to directly. Preferably, [Goto Clone-from-GitHub](#clonepsi4public); otherwise [Goto Clone-from-GitHub-Private](#clonepsi4private)


.. _`faq:binary`:

Binary Installer
----------------

* **Get Initially**

  Just go to http://www.psicode.org/downloads.html and follow the
  instructions there.

* **Build**

  Not applicable as binary is pre-built.

* **Get Updates**

  .. code-block:: bash

     >>> conda update psi4

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

  From the |PSIfour| repository at https://github.com/psi4/psi4, hit the
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


.. _`faq:psi3sourceforge`:

Psi3 from SourceForge
---------------------

* **Get Initially**

  A tarball of the most recent version of Psi3 (3.4.0 circa 2009) is
  available from `SourceForge
  <http://sourceforge.net/projects/psicode/files/psi/3.4.0/>`_

* **Build**

  Follow the ``INSTALL`` file that comes with the distribution. An old
  computer is probably handy for generating a working executable.

* **Get Updates**

  Updates are not forthcoming.

* **Contribute Back**

  This code is not under any development.


.. _`faq:githubworkflow`:

What is the suggested GitHub workflow
-------------------------------------

.. image:: /prflow.001.jpeg
.. image:: /prflow.002.jpeg
