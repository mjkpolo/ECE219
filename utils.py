import numpy as np
from functools import reduce
from sympy import (
    cos,
    sin,
    atan2,
    acos,
    pi,
    simplify,
    Matrix,
    integrate,
    Derivative,
    sqrt,
    Symbol,
)
from sympy.vector import CoordSys3D, gradient, divergence, curl, point, express, ParametricRegion, vector_integrate
from sympy.abc import r,phi,theta,x,y,z

height = Symbol('height')
radius = Symbol('radius')

"""
SPHERICAL
"""

###### system ######
sph = CoordSys3D(
    "sph",
    variable_names=("r", "theta", "phi"),
)
###### vectors ######
sph_to_cyl = Matrix(
    (
        (sin(sph.theta), cos(sph.theta), 0),
        (0, 0, 1),
        (cos(sph.theta), -sin(sph.theta), 0),
    )
)
sph_to_cart = Matrix(
    (
        (
            sin(sph.theta) * cos(sph.phi),
            cos(sph.theta) * cos(sph.phi),
            -sin(sph.phi),
        ),
        (
            sin(sph.theta) * sin(sph.phi),
            cos(sph.theta) * sin(sph.phi),
            cos(sph.phi),
        ),
        (cos(sph.theta), -sin(sph.theta), 0),
    )
)

"""
CARTESIAN
"""

###### system ######
cart = CoordSys3D(
    "cart",
    rotation_matrix=sph_to_cart,
    variable_names=("x", "y", "z"),
    parent=sph,
)
###### vectors ######
cart_to_sph = Matrix(
    (
        (
            sin(sph.theta) * cos(sph.phi),
            sin(sph.theta) * sin(sph.phi),
            cos(sph.theta),
        ),
        (
            cos(sph.theta) * cos(sph.phi),
            cos(sph.theta) * sin(sph.phi),
            -sin(sph.theta),
        ),
        (-sin(sph.phi), cos(sph.phi), 0),
    )
)
cart_to_cyl = Matrix(
    (
        (cos(sph.phi), sin(sph.phi), 0),
        (-sin(sph.phi), cos(sph.phi), 0),
        (0, 0, 1),
    )
)

"""
CYLINDRICAL
"""

###### system ######
cyl = CoordSys3D(
    "cyl",
    rotation_matrix=sph_to_cyl,
    variable_names=("r", "phi", "z"),
    parent=sph,
)
###### vectors ######
cyl_to_sph = Matrix(
    (
        (sin(sph.theta), 0, cos(sph.theta)),
        (cos(sph.theta), 0, -sin(sph.theta)),
        (0, 1, 0),
    )
)
cyl_to_cart = Matrix(
    (
        (cos(sph.phi), -sin(sph.phi), 0),
        (sin(sph.phi), cos(sph.phi), 0),
        (0, 0, 1),
    )
)

"""
Variable Transforms
"""
var_subs = {
    cart: {
        sph: lambda x, y, z: (
            sqrt(x**2 + y**2 + z**2),
            acos(z / sqrt(x**2 + y**2 + z**2)),
            atan2(y, x),
        ),
        cyl: lambda x, y, z: (
            sqrt(x**2 + y**2),
            atan2(y, x),
            z,
        ),
        cart: lambda x, y, z: (x, y, z),
    },
    sph: {
        cart: lambda r, theta, phi: (
            r * cos(theta) * sin(phi),
            r * sin(theta) * sin(phi),
            r * cos(phi),
        ),
        cyl: lambda r, theta, phi: (
            r * sin(theta),
            phi,
            r * cos(theta),
        ),
        sph: lambda r, theta, phi: (r, theta, phi),
    },
    cyl: {
        cart: lambda r, phi, z: (
            r * cos(phi),
            r * sin(phi),
            z,
        ),
        sph: lambda r, phi, z: (
            sqrt(r**2 + z**2),
            acos(z / sqrt(r**2 + z**2)),
            phi,
        ),
        cyl: lambda r, phi, z: (r, phi, z),
    },
}


def var(vector, sys):
    return vector.subs(
        tuple(zip(cart.base_scalars(), var_subs[sys][cart](*sys.base_scalars())))).subs(tuple(zip(sph.base_scalars(), var_subs[sys][sph](*sys.base_scalars())))).subs(tuple(zip(cyl.base_scalars(), var_subs[sys][cyl](*sys.base_scalars()))))


def sub(vector, sys, **kwargs):
    vec_sys = reduce(
        lambda a, b: a if a == b else None,
        (key.system for key in vector.components.keys()),
    )
    if vec_sys is None:
        raise TypeError("pass a valid vector")

    if "r" in kwargs:
        if "theta" in kwargs and "phi" in kwargs:
            var_sys = sph
        elif "phi" in kwargs and "z" in kwargs:
            var_sys = cyl
    elif "x" in kwargs and "y" in kwargs and "z" in kwargs:
        var_sys = cart
    else:
        raise TypeError("Pass a valid coordinate set")

    coordinate = lambda new_sys: var_subs[var_sys][new_sys](
        *tuple(kwargs[var.name.split(".")[1]] for var in var_sys.base_scalars())
    )

    return var(
        express(
            vector,
            sys,
        ),
        vec_sys,
    ).subs(tuple(map(tuple, zip(vec_sys.base_scalars(), coordinate(vec_sys)))))


"""
MISC
"""
cylinder = ParametricRegion((r * cos(phi), r * sin(phi), z), (r, 0, radius), (phi, 0, 2 * pi), (z, 0, height))
shell = ParametricRegion(
    (
        radius * sin(theta) * cos(phi),
        radius * sin(theta) * sin(phi),
        radius * cos(theta),
    ),
    (theta, 0, pi),
    (phi, 0, 2 * pi))
sphere = ParametricRegion(
    (
        r * sin(theta) * cos(phi),
        r * sin(theta) * sin(phi),
        r * cos(theta),
    ),
    (theta, 0, pi),
    (phi, 0, 2 * pi),
    (r,0,radius))

