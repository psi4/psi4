import sys, os

import cProfile
import nose


sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

config = nose.config.Config(
        plugins=nose.config.NoPlugins()
        )

cProfile.runctx(
    """nose.core.runmodule("grendel_tests",
                           argv=["--attr=profile"],
                           config=config)
    """, 
    globals(), 
    locals(), 
    filename=sys.argv[0] if len(sys.argv) > 0 else "nosetests.profile"
)
