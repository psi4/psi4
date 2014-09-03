import sys
import os

if len(sys.argv) != 3:
    print 'ERROR: Needs top_srcdir and objdir as arguments.'
    sys.exit()

src = sys.argv[1]  # @top_srcdir@/doc/sphinxman
obj = sys.argv[2]  # objdir/doc/sphinxman
# first '/' indicates "absolute" path rel to sfnx source
# final '/../' is extra buffer directory btwn sfnx/source and sfnx
print '/' + os.path.relpath(src, obj) + '/../'
