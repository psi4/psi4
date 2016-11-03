

* Making a release or prerelease

  * update all three fields of metadata.py (see below ex)
  * commit
  * on master, make annotated tag starting with v. note the hash.
  * run enough of build to make sure tag formatted properly and "Defining {} version" sane

# MUST update metadata.py on same commit at which make tag

git tag -a v1.0 d2243ef
git push origin v1.0
git push [remote] [tagname]

do NOT use the GitHub tagging interface - it creates lightweight tags

:source:`psi4/metadata.py`

__version__ = '1.0'
__version_long = '1.0+d2243ef'
__version_upcoming_annotated_v_tag = '1.1a1'

# Example current and upcoming tag pairs
# * '1.0', '1.1a1'
# * '1.0.4', '1.0.5' on maintenance branch 1.0.x
# * '1.1a1', '1.1b1'
# * '1.1rc2', '1.1'
# * '1.1', '1.2'


