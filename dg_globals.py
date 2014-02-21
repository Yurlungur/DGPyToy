"""
DGPyToy/dg_globals.py


Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
Time-stamp: <2014-02-20 21:50:54 (jonah)>

This module contains global constants, etc. for the DG methods toy.
"""


# Global Constants
# ----------------------------------------------------------------------
NUM_FACES = 2 # We're in one dimension, so each element has two faces,
              # a left and right.
NODES_PER_FACE = 1 # Each face is made up of one vertex.
NODETOL = 1E-15 # Tolerance within which we call something zero. Set
                # to machine epsilon.

# ----------------------------------------------------------------------

