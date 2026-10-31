"""
Builds (if needed) and runs the main program.
"""

import sys

from . import build
from . import frontend
from . import paths

def main():
    # Compile me (if needed).
    try:
        build.build(must=False)
    except build.BuildError:
        return 2
    # Run me.
    return frontend.run()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        raise RuntimeError(f"{paths.pyname()} does not have command line args")
    sys.exit(main())
