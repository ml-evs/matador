# encoding: utf-8
""" This file implements some useful geometric functions
for the construction and manipulation of convex hulls.
"""
import numpy as np

def points2plane(points):
    """ Convert points (xi, yi, zi) for i=1,..,3 into the
    equation of the plane spanned by the vectors v12, v13.
    For unit vectors e(i):

    v12 x v13 = n = i*e(1) + j*e(2) + k*e(3)

    and so the equation of the plane is

    i*x + j*y + k*z + d = 0.

    Input:

        points = [np.array([x1, y1, z1]), ...,  np.array([x3, y3, z3])].

    Returns:

        equation = np.array([i, j, k, d]).

    """
    v12 = points[1] - points[0]
    v13 = points[2] - points[0]
    normal = np.cross(v12, v13)
    d = -np.sum(normal)
    equation = np.append(normal, d)
    return equation


def barycentric2cart(structures):
    """ Convert ternary (x, y) in A_x B_y C_{1-x-y}
    to positions projected onto 2D plane.

    Input structures array is of the form:

        [
         [l(1)_0, l(2)_0 Eform_0],
         ...,
         [l(1)_n, l(2)_n, Eform_n]
        ]


    where l3 = 1 - l2 - l1 are the barycentric coordinates of the point
    in the triangle defined by the chemical potentials.
    """
    structures = np.asarray(structures)
    cos30 = np.cos(np.pi/6)
    cos60 = np.cos(np.pi/3)
    coords = np.zeros_like(structures)
    coords[:, 0] = structures[:, 0] + structures[:, 1] * cos60
    coords[:, 1] = structures[:, 1] * cos30
    coords[:, 2] = structures[:, -1]
    return coords
