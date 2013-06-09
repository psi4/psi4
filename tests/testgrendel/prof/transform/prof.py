
import sys
import os
import cProfile
import nose

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel_tests import profile_suite, profile_subsuites

suite_to_run = profile_subsuites[sys.argv[1]]

file = sys.argv[1] if len(sys.argv) > 1 else "nosetests"

cProfile.run(
    """nose.core.run(suite=suite_to_run)
    """, 
    file+'.pyprof'
)
print("logged to {}".format(file+'.pyprof'))

os.system('pyprof2calltree -i {} -o {}'.format(file+'.pyprof', file+'.calltree')) 
