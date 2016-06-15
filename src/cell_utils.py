""" This file implements some useful functions
for cell manipulation.
"""

from math import pi, cos, sin, sqrt, acos


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
    if cx < tol:
        cx = 0.0
    cy = c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
    if abs(cy) < tol:
        cy = 0.0
    cz = sqrt(c**2 - cx**2 - cy**2)
    lattice_cart.append([bx, by, 0.0])
    lattice_cart.append([cx, cy, cz])
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
        cos_alpha = vec_b[i] * vec_c[i]
        cos_beta = vec_c[i] * vec_a[i]
        cos_gamma = vec_a[i] * vec_b[i]
    cos_alpha /= b*c
    cos_beta /= c*a
    cos_gamma /= a*b
    alpha = 180.0*acos(cos_alpha)/pi
    beta = 180.0*acos(cos_beta)/pi
    gamma = 180.0*acos(cos_gamma)/pi
    lattice_abc.append([alpha, beta, gamma])
    return lattice_abc
