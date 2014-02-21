#!/usr/bin/env python2

"""
DTPyToy/meshgen.py

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
Time-stamp: <2014-02-20 23:03:09 (jonah)>

This small script contains a functions to initialise and run a grid for DG methods. 

The setup is needlessly complicated in one dimension, but it will
prove useful to learn how to use it.
"""

# Imports
# ----------------------------------------------------------------------
import numpy as np
from numpy import linalg
from dg_globals import *
# ----------------------------------------------------------------------


def generate_mesh(xmin,xmax,k):
    """
    Generates a grid for DG methods in one dimension. The grid
    describes a set of k non-overlapping elements that partition the
    interval [xmin,xmax]. The elements are described by a pair of vertices
    [xk_l, xk_r] which define the left and right edges of an element.

    This method returns two numpy arrays.

    The first array, called vx, contains the physical positions of all
    vertices. We call the length num_vert for the number of vertices.

    The second array, called e_to_v, maps each element to the vertices
    which form its faces. This array is k x num_vert. The row represents the
    element number and the left and right columns are the indices in
    vx of the left and right faces respectively.
    """
    # We always have K+1 vertices
    num_vert = k+1
    # Generate vx.
    vx = np.linspace(xmin,xmax,num_vert)
    # Generate e_to_v
    e_to_v = np.array([(i,i+1) for i in range(len(vx)-1)])
    # Give it all back
    return vx,e_to_v


def generate_x(vx,e_to_v,r):
    """
    Takes the vector of physical positions of the vertives, vx, the
    array mapping the elements to their vertices, e_to_v, and the
    elementwise affine vector of nodes, r, compute the physical
    position of every node in the grid.

    The output array, x, is Np x K where Np is the number of grid
    points per element and K is the number of elements. Each column
    corresponds to an element and each row corresponds to the index of
    the node in the element.

    We can think of x as a mapping from affine coordinates (indexed by
    element) into physical positions.
    """
    # The number of nodes per element
    num_points = len(r)
    # The left and right vertices for each element
    x_left = e_to_v[...,0]
    x_right = e_to_v[...,1]
    x = np.outer(np.ones(num_points),vx[x_left]) + 0.5*np.outer((r+1),(vx[x_right]-vx[x_left]))
    return x


def geometric_factors(x,Dr):
    """
    Takes the differentiation matrix Dr and positions array x and
    calculates the local metric on each node on each element. The
    local metric is rx, which is one over the Jacobian of the
    transformation. The jacobian J is also returned for good measure.
    """
    J = Dr*x
    rx = 1.0/J
    return rx,J


def make_face_map(r,x):
    """
    Generates an array that maps from an element to the left or right
    face node of an element. The array that's returned is 2xK. So the
    top row is left boundaries and the bottom node is right
    boundaries.
    """
    face_mask_left = np.argwhere(np.abs(r+1)<NODETOL)[0,0]
    face_mask_right = np.argwhere(np.abs(r-1)<NODETOL)[0,0]
    return x[[face_mask_left,face_mask_right],...]


def make_normals(num_elements):
    """
    Generates an array of outward pointing normals at each face. In
    one dimension, an outward pointing normal is just a scalar
    value. It's -1 if the face is at the left end of the element and
    +1 if the face is at the right end of the element.
    """
    # This will be the vector of normals
    nx = zeros(NODES_PER_FACE*NUM_FACES,num_elements)
    nx[0,...] = -1.0
    nx[1,...] = 1.0
    return nx



