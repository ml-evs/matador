""" This file implements some useful functions
for cell manipulation.
"""

from math import pi, cos, sin, sqrt, acos, log10
import numpy as np


def abc2cart(lattice_abc):
    """ Convert lattice_abc=[[a,b,c],[alpha,beta,gamma]]
    (in degrees) to lattice vectors
    lattice_cart=[[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]].
    """
    a = lattice_abc[0][0]
    b = lattice_abc[0][1]
    c = lattice_abc[0][2]
    deg2rad = pi/180
    alpha = lattice_abc[1][0] * deg2rad
    beta = lattice_abc[1][1] * deg2rad
    gamma = lattice_abc[1][2] * deg2rad
    lattice_cart = []
    lattice_cart.append([a, 0.0, 0.0])
    # vec(b) = (b cos(gamma), b sin(gamma), 0)
    bx = b*cos(gamma)
    by = b*sin(gamma)
    tol = 1e-12
    if abs(bx) < tol:
        bx = 0.0
    if abs(by) < tol:
        by = 0.0
    cx = c*cos(beta)
    if abs(cx) < tol:
        cx = 0.0
    cy = c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
    if abs(cy) < tol:
        cy = 0.0
    cz = sqrt(c**2 - cx**2 - cy**2)
    lattice_cart.append([bx, by, 0.0])
    lattice_cart.append([cx, cy, cz])
    assert(cart2abc(lattice_cart) != lattice_abc)
    return lattice_cart


def cart2abc(lattice_cart):
    """ Convert lattice_cart =[[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]]
    to lattice_abc=[[a,b,c],[alpha,beta,gamma]].
    """
    vec_a = lattice_cart[0]
    vec_b = lattice_cart[1]
    vec_c = lattice_cart[2]
    lattice_abc = []
    a = 0
    b = 0
    c = 0
    for i in range(3):
        a += vec_a[i]**2
        b += vec_b[i]**2
        c += vec_c[i]**2
    a = sqrt(a)
    b = sqrt(b)
    c = sqrt(c)
    lattice_abc.append([a, b, c])
    # cos(alpha) = b.c /|b * c|
    cos_alpha = 0
    cos_beta = 0
    cos_gamma = 0
    for i in range(3):
        cos_alpha += vec_b[i] * vec_c[i]
        cos_beta += vec_c[i] * vec_a[i]
        cos_gamma += vec_a[i] * vec_b[i]
    cos_alpha /= b*c
    cos_beta /= c*a
    cos_gamma /= a*b
    alpha = 180.0*acos(cos_alpha)/pi
    beta = 180.0*acos(cos_beta)/pi
    gamma = 180.0*acos(cos_gamma)/pi
    lattice_abc.append([alpha, beta, gamma])
    return lattice_abc


def real2recip(real_lat):
    """ Convert the real lattice in Cartesian basis to
    the reciprocal space lattice.
    """
    real_lat = np.asarray(real_lat)
    recip_lat = np.zeros((3, 3))
    recip_lat[0] = (2*pi)*np.cross(real_lat[1], real_lat[2]) / \
        (np.dot(real_lat[0], np.cross(real_lat[1], real_lat[2])))
    recip_lat[1] = (2*pi)*np.cross(real_lat[2], real_lat[0]) / \
        (np.dot(real_lat[1], np.cross(real_lat[2], real_lat[0])))
    recip_lat[2] = (2*pi)*np.cross(real_lat[0], real_lat[1]) / \
        (np.dot(real_lat[2], np.cross(real_lat[0], real_lat[1])))
    return recip_lat


def calc_mp_spacing(real_lat, mp_grid):
    """ Convert real lattice in Cartesian basis and the
    kpoint_mp_grid into a grid spacing.
    """
    recip_lat = real2recip(real_lat)
    recip_len = np.zeros((3))
    recip_len = np.sqrt(np.sum(np.power(recip_lat, 2), axis=1))
    max_spacing = 0
    for j in range(3):
        spacing = recip_len[j] / (2*pi*mp_grid[j])
        max_spacing = (spacing if spacing > max_spacing else max_spacing)
    exponent = round(log10(max_spacing) - 1)
    return round(max_spacing + 0.5*10**exponent, 2)
