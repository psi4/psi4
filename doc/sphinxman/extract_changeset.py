#!/usr/bin/python
import os
import sys
# Extracts last git changelog item into html for psicode

gitlast = 'feed/latest_trac_changeset.txt'  # path to output file with latest changeset
githist = 'feed/history_trac_changeset.txt'  # path to output file with last Nhist changesets
Nhist = 100

# So that this can be called from anywhere
writedir = sys.argv[1] if len(sys.argv) == 2 else './'
thisdir = os.path.abspath(os.path.dirname(__file__))
os.chdir(thisdir)

log_cmd = 'git log --branches=master -%s --no-merges --pretty=format:\'<title> %%cr [%%h] by %%cn: <br/> %%s<br/></title>\' > %s'\
          % (Nhist, writedir + os.sep + githist)
os.system(log_cmd)

log_cmd = 'git log --branches=master -1 --no-merges --pretty=format:\'<title> %%cr by %%cn: <br/> %%s<br/></title>\' > %s'\
          % (writedir + os.sep + gitlast)
os.system(log_cmd)
