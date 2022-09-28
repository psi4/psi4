import ast
import shutil
import getpass
import pathlib
import operator


# How do pytest tests in `test_standard_suite.py` communicate with docs?
# * after calling psi4 for each test, `standard_suite_runner.py` appends a one-line dict describing the result to local file "stdsuite_psi4.txt".
# * thus, local "stdsuite_psi4.txt" files with unordered entries end up in pytest's basetemp/popen-*/ (if run parallel; otherwise, single file in basetemp/).
#   * to combine these multiple files, can use bash: `cat /tmp/pytest-of-$USER/pytest-current/*/stdsuite_psi4.txt > stdsuite_psi4.txt`
#   * or this script globs them.
# * if stdsuite_psi4.txt file is to be diff-ed (say, checked into a repository), it should be sorted to minimize churn. This script performs the sort.
# * `document_capabilities.py` reads a single arbitrary-order stdsuite_psi4.txt file (specify location with --stdsuite) and formats tables.

# * Thus, this script is largely unnecessary and piecemeal.
#   * if you're working locally, run pytest, merge the record files with `cat`, and run `document_capabilities` on the result.
#   * but if you're running the full standard suite to update the repo, this script is a handy intermediary or has the pieces to run separately.
#   * rationalization for redundancy of test_standard_suite, stdsuite_psi4.txt, and autodoc_capabilities tables all in repository
#     * from a single-source-of-truth argument, one should re-run stdsuite at docs generation time and never save stdsuite_psi4.txt or generated tables.
#     * but even run with filters and `-n 20`, the standard suite takes tens of minutes (and new capabilities aren't exactly weekly occurrences), so commit stdsuite_psi4.txt to repo.
#     * still from a fewest-sources-of truth argument, one should


# edit as needed
PYTEST_SCRATCH = f"/tmp/pytest-of-{getpass.getuser()}/pytest-current/"

# concatenate pytest fragments into a single local file
with open("stdsuite_psi4.txt","wb") as wfd:
    for fl in pathlib.Path(PYTEST_SCRATCH).glob("**/stdsuite_psi4.txt"):
        with open(fl, "rb") as fd:
            shutil.copyfileobj(fd, wfd)

# read lines of dicts into list of dicts and sort it. only stability matters, not the actual order
with open("stdsuite_psi4.txt", "r") as fp:
    contents = fp.readlines()
stuffs = [ast.literal_eval(ln) for ln in contents]
stuffs.sort(key=operator.itemgetter("method", "driver", "reference", "fcae", "scf_type", "corl_type", "module", "status", "note"))

# write out single file with sorted lines of dicts
with open("stdsuite_psi4.txt", "w") as fp:
    for stuff in stuffs:
        fp.write(f"{stuff!r}\n")
