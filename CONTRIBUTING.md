# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into Psi4.  As of February 2016, the
procedure for contributing code is exactly the same for the core development
team and for external contributors.

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free).
* Fork the psi4/psi4 repository on GitHub.
* Optionally, sign up for an account on [Travis CI](https://travis-ci.org/),
  which you can point to your fork of Psi4 to automatically test each of your
  commits to that fork.

## Making Changes

* On your local machine, clone your fork of the psi4 repo, which is located at
  github.com:yourusername/psi4.  More detailed instructions for forking Psi4
  can be found
  [here](https://github.com/psi4/psi4/wiki/1_Obtaining#fork-from-public-github-repository).
* Add some really awesome code.  It's usually a good idea to make changes on
  a new branch, but that is not mandatory.
* At any time during the development of your new feature, navigate to your fork
  of Psi4 on GitHub and open a pull request. Note that after you launch a pull
  request from one of your fork's branches, all subsequent commits to that
  branch will be added to the open pull request automatically.  Each commit
  added to the pull request will be validated for mergability, compilation and
  test suite compliance; the results of these tests will be visible on the pull
  request page.
* If you're providing a new feature, you must add test cases and documentation.
* When the code is ready to go, make sure you run the full test suite on your
  local machine to check that nothing is broken.
* When you see that all of the validation checks have passed on the pull
  request page, check the "Ready to go" box, to let the Psi4 team know that the
  changes are complete; the code will not be considered for merging until this
  box is checked.

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [GitHub pull request documentation](https://help.github.com/send-pull-requests/)
* [A guide to contributing to software packages](http://www.contribution-guide.org)
