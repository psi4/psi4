# <img src="https://github.com/psi4/psi4media/blob/master/logos-psi4/psi4square.png" height=150>

| **Status** | [![Azure DevOps builds](https://img.shields.io/azure-devops/build/psi4/e80489d7-9619-4512-8e7b-255e355b3ab8/1?logo=azure%20devops)](https://dev.azure.com/psi4/psi4/_build?definitionId=1) [![Codecov coverage](https://img.shields.io/codecov/c/github/psi4/psi4.svg?logo=Codecov&logoColor=white)](https://codecov.io/gh/psi4/psi4) |
| :------ | :------- |
| **Latest Release** | [![Last release tag](https://img.shields.io/github/release/psi4/psi4.svg)](https://github.com/psi4/psi4/releases)  [![Commits since release](https://img.shields.io/github/commits-since/psi4/psi4/v1.8.svg)](https://github.com/psi4/psi4/releases/tag/v1.8) [![python](https://img.shields.io/badge/python-3.8%2C%203.9%2C%203.10%2C%203.11-blue.svg)](https://psicode.org/psi4manual/master/introduction.html#supported-systems) |
| **Communication** | [![User site](https://img.shields.io/badge/home-Psi4-5077AB.svg)](https://psicode.org/) [![docs latest](https://img.shields.io/badge/docs-latest-5077AB.svg?logo=read%20the%20docs)](https://psicode.org/psi4manual/master/index.html) [![chat on forum](https://img.shields.io/badge/chat-on_forum-808493.svg?logo=Discourse&logoColor=white)](http://forum.psicode.org/) [![dev chat on slack](https://img.shields.io/badge/dev_chat-on_slack-808493.svg?logo=slack)](https://join.slack.com/t/psi4/shared_invite/zt-5s36s4rb-SQH6_AWyfWOqlKYN3cFs4Q) |
| **Foundation** | [![license](https://img.shields.io/github/license/psi4/psi4.svg)](https://opensource.org/licenses/LGPL-3.0) [![platforms](https://img.shields.io/badge/Platforms-Linux%2C%20MacOS%2C%20MacOS%20Silicon%2C%20Windows%2C%20Windows%20WSL-orange.svg)](https://psicode.org/psi4manual/master/introduction.html#supported-systems) [![python](https://img.shields.io/badge/python-3.8%2C%203.9%2C%203.10%2C%203.11-blue.svg)](https://psicode.org/psi4manual/master/introduction.html#supported-systems) |
| **Installation** | [![obtain latest](https://img.shields.io/badge/obtain-latest-green.svg)](https://psicode.netlify.com/installs/latest) [![Conda](https://img.shields.io/conda/v/conda-forge/psi4.svg)](https://anaconda.org/conda-forge/psi4) [![Anaconda-Server Badge](https://anaconda.org/conda-forge/psi4/badges/latest_release_relative_date.svg)](https://anaconda.org/conda-forge/psi4) |
| **Demo** | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/psi4/psi4/56fbc7787af67dabdf1897d0dfe4263d8d97e241?urlpath=lab%2Ftree%2Fdoc%2Fsphinxman%2Fsource%2Fpsiapi.ipynb) |

<!--  -->
<!-- [![Last release date](https://img.shields.io/github/release-date/psi4/psi4.svg)](https://github.com/psi4/psi4/releases) -->
<!-- [![Anaconda-Server Badge](https://anaconda.org/psi4/psi4/badges/version.svg)](https://anaconda.org/psi4/psi4) -->

<!--<a href="https://psi4.slack.com/messages"> <img src="https://img.shields.io/badge/dev_chat-on_slack-808493.svg" /></a>
<a href="mailto:psi4aiqc+slackinvite@gmail.com?subject=request slack invite (incl. who, where, email)"> <img src="https://img.shields.io/badge/dev_chat-invite-808493.svg" /></a> -->

<!--[![Anaconda-Server Badge](https://anaconda.org/psi4/psi4/badges/installer/conda.svg)](https://anaconda.org/psi4/psi4) 
[![Anaconda-Server Badge](https://anaconda.org/psi4/psi4/badges/platforms.svg)](https://anaconda.org/psi4/psi4) -->

<!--
| **PR Activity** | 
[![commit activity](https://img.shields.io/github/commit-activity/y/psi4/psi4.svg)](https://github.com/psi4/psi4/graphs/contributors) 
[![issues-pr-closed](https://img.shields.io/github/issues-pr-closed-raw/psi4/psi4.svg)](https://github.com/psi4/psi4/pulls)
-->

Psi4 is an open-source suite of *ab initio* quantum chemistry programs
designed for efficient, high-accuracy simulations of
molecular properties. We routinely perform computations with >2500 basis functions on multi-core machines.

With computationally demanding portions written in C++, exports
of many C++ classes into Python via Pybind11, and a flexible Python driver, Psi4
strives to be friendly to both users and developers.

* **Users' Website**  www.psicode.org

* **Downloading and Installing Psi4** https://psicode.org/psi4manual/master/build_faq.html (for the CMake adept, see [CMakeLists.txt](CMakeLists.txt)

* **Manual**  [http://bit.ly/psi4manual](https://psicode.org/psi4manual/master/index.html) (built nightly from master branch) or https://psicode.org/psi4manual/1.4.0/index.html (last release)

* **Tutorial** https://psicode.org/psi4manual/master/tutorial.html for Psithon (``psi4 job.in``), https://psicode.org/psi4manual/master/psiapi.html for PsiAPI (``python job.py``)

* **Forum** http://forum.psicode.org

* **Communication & Support** https://psicode.org/psi4manual/master/introduction.html#technical-support

* **GitHub**  https://github.com/psi4/psi4 (authoritative repository)

* **Continuous Integration Status** [![Azure DevOps builds](https://img.shields.io/azure-devops/build/psi4/e80489d7-9619-4512-8e7b-255e355b3ab8/1/master.svg?logo=azure%20devops)](https://dev.azure.com/psi4/psi4/_build?definitionId=1) on Linux and Windows

* **Anaconda**  https://anaconda.org/psi4 (binary available for Linux, Mac, Mac Silicon, Windows, and WSL Windows [![Binstar Badge](https://anaconda.org/psi4/psi4/badges/downloads.svg)](https://anaconda.org/psi4/psi4) ) [![Binstar Badge](https://anaconda.org/conda-forge/psi4/badges/downloads.svg)](https://anaconda.org/conda-forge/psi4) ) [instructions](https://psicode.org/psi4manual/master/conda.html#how-to-install-a-psi4-binary-with-the-psi4conda-installer-download-site)

* **Coverage** Python and C++ source code lines hit by running most of the test suite. [![codecov](https://img.shields.io/codecov/c/github/psi4/psi4.svg?logo=Codecov&logoColor=white)](https://codecov.io/gh/psi4/psi4)

* **Interested Developers**  https://psicode.org/developers.php (replacement page needed) (welcome to fork psi4/psi4 and follow [GitHub contribution procedure](https://psicode.org/psi4manual/master/build_obtaining.html#faq-githubworkflow)) [![PRs welcome](https://img.shields.io/badge/PRs-welcome-yellow.svg)](http://makeapullrequest.com)

* **Sample Inputs**  http://www.psicode.org/psi4manual/master/testsuite.html (also in [`samples/`](samples))

* **Download Tarball** https://github.com/psi4/psi4/releases 

<!--* **Build Dashboard** https://testboard.org/cdash/index.php?project=Psi

* **YouTube Channel** https://www.youtube.com/psitutorials-->


License [![license](https://img.shields.io/github/license/psi4/psi4.svg)](https://opensource.org/licenses/LGPL-3.0)
=======

Psi4: an open-source quantum chemistry software package

Copyright (c) 2007-2023 The Psi4 Developers.

The copyrights for code used from other parties are included in
the corresponding files.

Psi4 is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, version 3.

Psi4 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with Psi4; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

The full text of the GNU Lesser General Public License (version 3) is included in the
COPYING.LESSER file of this repository, and can also be found
[here](https://www.gnu.org/licenses/lgpl.txt).


Citation [![doi](https://img.shields.io/badge/doi-10.1063/5.0006002-5077AB.svg)](https://doi.org/10.1063/5.0006002)
========

The journal article reference describing Psi4 is:

D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
"Psi4 1.4: Open-Source Software for High-Throughput Quantum Chemistry",
J. Chem. Phys. 152(18) 184108 (2020).

* [![doi](https://img.shields.io/badge/doi-10.1021/acs.jctc.7b00174-5077AB.svg)](https://doi.org/10.1021/acs.jctc.7b00174) for Psi4 v1.1
* [![doi](https://img.shields.io/badge/doi-10.1021/acs.jctc.8b00286-5077AB.svg)](https://doi.org/10.1021/acs.jctc.8b00286) for Psi4NumPy
* [![doi](https://img.shields.io/badge/doi-10.1002/wcms.93-5077AB.svg)](https://doi.org/10.1002/wcms.93) for Psi4 alpha releases
* [![doi](https://img.shields.io/badge/doi-10.1002/jcc.20573-5077AB.svg)](https://doi.org/10.1002/jcc.20573) for Psi3
