# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into Psi4.

**Working on your first Pull Request?** You can learn how from
this *free* series [How to Contribute to an Open Source Project on
GitHub](https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project)

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free).
* [Fork](https://help.github.com/articles/fork-a-repo/) the
  [psi4/psi4](https://github.com/psi4/psi4) repository on GitHub.
* On your local machine,
  [clone](https://help.github.com/articles/cloning-a-repository/) your fork of
  the Psi4 repository.
* More detailed instructions for interacting with your Psi4 fork can be found
  [here](https://psicode.org/psi4manual/master/build_obtaining.html#faq-forkpsi4public).
  and [here](https://psicode.org/psi4manual/master/build_obtaining.html#faq-githubworkflow).

## Making Changes

* Add some really awesome code to your local fork.  It's usually a [good
  idea](https://jmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/)
  to make changes on a
  [branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/)
  with the branch name relating to the feature you are going to add.
* When you are ready to see CI tests and for others to examine and comment on your new feature,
  navigate to your fork of Psi4 on GitHub and open a [pull
  request](https://help.github.com/articles/using-pull-requests/) (PR). Note that
  after you launch a PR from one of your fork's branches, all
  subsequent commits to that branch will be added to the open pull request
  automatically.  Each commit added to the PR will be validated for
  mergability, compilation and test suite compliance; the results of these tests
  will be visible on the PR page.
* If you're providing a new feature, you must add test cases and documentation.
* When the code is ready to go, make sure you run the full or relevant portion of the
  [test suite](https://psicode.org/psi4manual/master/build_planning.html#faq-subsettests)
  on your local machine to check that nothing is broken.
* When you're ready to be considered for merging, issue a PR comment with only `/review-ready`
  as content to let the Psi4 team know that the changes are ready for review.
  The code will not be merged until the continuous
  integration (GiHub Actions and Azure DevOps) returns checkmarks,
  and multiple core developers give "Approved" reviews.

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [PR best practices](https://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
* [A guide to contributing to software packages](https://www.contribution-guide.org)
