# GHA for Psi4

### docs.yml

* Runs on: push to master
* Results
  * automated commit of built html docs to psi4/psi4docs:master, which in turn is served up by netlify to https://psi4manual.netlify.app/, which in turn is [redirected by psicode](https://github.com/psi4/psicode-hugo-website/blob/master/netlify.toml) into https://psicode.org/psi4manual/master/index.html 
* Goals
  * doxygen docs build
  * sphinx docs build ✔️, build without warnings ✔️, build without link errors ❌, build nit-picky ❌
  * nightly-build docs are available
  
----

### docs-pr.yml

* Runs on: push to PR (actually, [`pull_request_target`](https://securitylab.github.com/research/github-actions-preventing-pwn-requests/), so to edit workflow, change docs-pr.yml on LHS of "psi4/psi4:<branch> <- <user>/psi4:<prbranch>")
* Results
  * tarball of html submitter can download and inspect offline
  * automated commit of repo changes to PR (`samples/` and a few C++/Py sync files like psifiles and physconst) to PR branch. this commit is still to be reviewed, so not a security risk
  * automated PR of built docs changes to psi4/psi4docs:p4-<prnumber>, which generates a docs website preview. this PR is submitted as draft and never merged, so not a security risk
* Goals
  * PR changes don't break docs by syntax or references
  * submitter can see how PR changes appear in docs without building locally
  * psi4 repo remains up-to-date for generated files
  * psi4 builds against conda-forge environments
