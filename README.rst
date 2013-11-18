PSI4: Ab Initio Quantum Chemistry
---------------------------------

* Users' Website: www.psicode.org

* Download Tarball: http://sourceforge.net/projects/psicode/ (NSF likes to see download tallies which GitHub doesn't provide)

* Developers' GitHub Wiki: https://github.com/psi4/psi4/wiki (currently private, transition soon)

* Manual: http://sirius.chem.vt.edu/psi4manual/master/index.html (built nightly from master branch)

PSI4 is an open-source suite of ab initio quantum chemistry programs designed for efficient, 
high-accuracy simulations of a variety of molecular properties. We can routinely perform 
computations with more than 2500 basis functions running serially or with modest parallelism.

With computationally demanding portions written in C++, Boost exports of many C++ classes into 
Python, and a flexible Python driver, PSI4 strives to be friendly to both users and developers.

[11 Nov 2013] We are transitioning from main development on a private GitHub repository toward
main development on this public repository. Sensitive code (prone to scientific scooping) will
remain on the private branch, but otherwise we anticipate that contributions from both the primary
developers and any interested persons shall be made here. As of today, the public and private 
repos are synced, but this isn't necessarily the most stable version. A new stable, fully
tested beta6 should appear before the end of the year.
