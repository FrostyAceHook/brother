"""
Runs the rocket analysis.
"""

import sys

from .. import build
from . import rocket
from .. import paths

def main():
    try:
        build.build(must=False)
    except build.BuildError:
        return 2
    return rocket.run()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        raise RuntimeError(f"{paths.pyname()} does not have command line args")
    sys.exit(main())
