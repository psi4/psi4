
import sys
import os


# Append the directory with the grendel_tests module in it
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

# Append the directory with grendel in it
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))


import nose
import grendel_tests

this_scripts_suite = grendel_tests.standard_suite

if __name__ == "__main__":
    nose.core.run(
        suite=this_scripts_suite
    )

