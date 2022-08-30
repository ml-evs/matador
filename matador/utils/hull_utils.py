# coding: utf-8
# Distributed under the terms of the MIT License.

""" This file implements some useful geometric functions for the
construction and manipulation of convex hulls.

"""

import warnings

import numpy as np

EPS = 1e-12


def vertices2plane(points):
    """Convert points (xi, yi, zi) for i=1,..,3 into the
    equation of the plane spanned by the vectors v12, v13.
    For unit vectors e(i):

    v12 x v13 = n = i*e(1) + j*e(2) + k*e(3)

    and so the equation of the plane is

    i*x + j*y + k*z + d = 0.

    Parameters:
        points (list of np.ndarray): list of 3 3D numpy arrays containing
            the points comprising the vertex.

    Returns:
        callable: a function which will return the vertical distance between
            the point and the plane:

    """
    v12 = points[1] - points[0]
    v13 = points[2] - points[0]
    normal = np.cross(v12, v13)
    d = -np.sum(np.dot(normal, points[0]))
    # check other points are on the plane, to some precision
    assert np.abs(np.dot(normal, points[2]) + d) < 0 + EPS
    assert np.abs(np.dot(normal, points[1]) + d) < 0 + EPS

    def get_height_above_plane(structure):
        """Find the z-coordinate on the plane matching
        the (x, y) coordinates of the structure, then calculate
        the difference between this z and the z of the point given.
        """
        x = structure[0]
        y = structure[1]
        z = structure[2]
        if np.abs(normal[2]) < EPS:
            warnings.warn(
                f"Normal of plane {normal} is ill-defined. Returning 0 for height above plane."
            )
            return 0
        z_plane = -((x * normal[0] + y * normal[1] + d) / normal[2])
        height = z - z_plane
        return height

    return get_height_above_plane


def vertices2line(points):
    """Perform a simple linear interpolation on
    two points.

    Parameters:
        points (list of np.ndarray): list of two 2D numpy arrays.
            of form [[x1, E1], [x2, E2]].

    Returns:
        (float, float): a tuple containing the gradient and
            intercept of the line intersecting the two points.

    """
    energy_pair = [points[0][1], points[1][1]]
    comp_pair = [points[0][0], points[1][0]]
    gradient = (energy_pair[1] - energy_pair[0]) / (comp_pair[1] - comp_pair[0])
    intercept = (
        (energy_pair[1] + energy_pair[0]) - gradient * (comp_pair[1] + comp_pair[0])
    ) / 2

    return gradient, intercept


def is_point_in_triangle(point, triangle, preprocessed_triangle=False):
    """Check whether a point is inside a triangle.

    Parameters:
        point (np.ndarray): 3x1 array containing the coordinates of the
            point.
        triangle (np.ndarray): 3x3 array specifying the coordinates of
            the triangle vertices.

    Keyword arguments:
        preprocessed_triangle (bool): if True, treat the input triangle
            as already processed, i.e. the array contains the inverse of
            the barycentric coordinate array.

    Returns:
        bool: whether or not the point is found to lie inside the
            triangle. If all vertices of the triangle lie on the same
            line, return False.

    """

    if not preprocessed_triangle:
        cart_planes = barycentric2cart(triangle).T
        cart_planes[-1, :] = 1
        if np.linalg.det(cart_planes) == 0:
            return False
        cart_plane_inv = np.linalg.inv(cart_planes)
    else:
        cart_plane_inv = triangle

    barycentric_structure = barycentric2cart(point.reshape(1, 3)).T
    barycentric_structure[-1, :] = 1
    plane_barycentric_structure = cart_plane_inv @ barycentric_structure

    return (plane_barycentric_structure >= 0 - EPS).all()


def barycentric2cart(structures):
    """Convert ternary (x, y) in A_x B_y C_{1-x-y}
    to positions projected onto 2D plane.

    Input structures array is of the form:

        [
            [l(1)_0, l(2)_0, Eform_0],
            [l(1)_n, l(2)_n, Eform_n]
        ]

    where l3 = 1 - l2 - l1 are the barycentric coordinates of the point
    in the triangle defined by the chemical potentials.

    Parameters:
        structures (list of np.ndarray): list of 3D numpy arrays containing
            input points.

    Returns:
        list of np.ndarray: list of numpy arrays containing converted
            coordinates.

    """
    structures = np.asarray(structures)
    cos30 = np.cos(np.pi / 6)
    cos60 = np.cos(np.pi / 3)
    coords = np.zeros_like(structures)
    coords[:, 0] = structures[:, 0] + structures[:, 1] * cos60
    coords[:, 1] = structures[:, 1] * cos30
    coords[:, 2] = structures[:, -1]

    return coords


class FakeHull:
    """Implements a thin class to mimic a ConvexHull object
    that would otherwise be undefined for two points."""

    def __init__(self):
        """Define the used hull properties."""
        self.vertices = [0, 1]
        self.simplices = []
