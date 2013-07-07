
import sys, os

import cProfile
import nose


sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

cProfile.runctx('nose.core.runmodule("grendel_tests")', globals(), locals(), filename="all_nosetests.profile")
