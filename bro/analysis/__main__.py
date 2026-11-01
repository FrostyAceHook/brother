"""
Runs './rocket.py' main.
"""

from .. import geez
from .rocket import main

if __name__ == "__main__":
    with geez.instance():
        main()
