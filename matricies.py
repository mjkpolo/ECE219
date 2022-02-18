import numpy as np
import sympy

theta, phi = sympy.symbols("theta phi")


def lineIntegral(charge_per_meter, dl, lower, upper, rprime, r):
    r = r.cart()
    rprime = rprime.cart(r)
    diff = (r[0] - rprime[0], r[1] - rprime[1], r[2] - rprime[2])
    magnitude = sympy.sqrt(diff[0] ** 2 + diff[1] ** 2 + diff[2] ** 2)
    coeff = charge_per_meter / (4 * sympy.pi * 8.85e-12 * magnitude**3)
    return (
        float(sympy.integrate(diff[0] * coeff, (dl, lower, upper))),
        float(sympy.integrate(diff[1] * coeff, (dl, lower, upper))),
        float(sympy.integrate(diff[2] * coeff, (dl, lower, upper))),
    )


def getPhi(x, y):
    if x == 0:
        return sympy.pi / 2 if y > 0 else -sympy.pi / 2
    elif x < 0:
        return sympy.atan(y / x) + sympy.pi
    else:
        return sympy.atan(y / x)


def evaluate_vector(vec, coord):
    v = list(vec)
    for i in range(len(v)):
        x, y, z = coord.cart()
        r, theta, phi = coord.sph()
        v[i] = v[i].subs([("x", x), ("y", y), ("z", z), ("phi", phi), ("theta", theta)])
        if coord.sys == "cyl":
            r = coord.cyl()[0]
        v[i] = v[i].subs("r", r)
    return tuple(float(f) for f in v)


class Coord:
    def __init__(self, coord, sys=None):
        if sys != "cart" and sys != "sph" and sys != "cyl":
            raise ValueError("Pass sys as 'cart', 'sph', or 'cyl'")
        self.sys = sys
        self.coord = np.array(coord)

    def __iter__(self):
        return iter(self.coord)

    def __getitem__(self, index):
        return self.coord[index]

    def __repr__(self):
        return f"Coord({self.coord}, sys={self.sys})"

    def __str__(self):
        return str(self.coord)

    def cart(self):
        if self.sys == "cart":
            return self
        elif self.sys == "sph":
            r, theta, phi = self.coord
            x = r * sympy.sin(theta) * sympy.cos(phi)
            y = r * sympy.sin(theta) * sympy.sin(phi)
            z = r * sympy.cos(theta)
        elif self.sys == "cyl":
            r, phi, z = self.coord
            x = r * sympy.cos(phi)
            y = r * sympy.sin(phi)
        else:
            raise ValueError("Yell at Marco")
        return Coord((x, y, z), sys="cart")

    def sph(self):
        if self.sys == "sph":
            return self
        elif self.sys == "cart":
            x, y, z = self.coord
            theta = sympy.acos(z / sympy.sqrt(x**2 + y**2 + z**2))
            phi = getPhi(x, y)
            r = sympy.sqrt(x**2 + y**2 + z**2)
        elif self.sys == "cyl":
            r, phi, z = self.coord
            theta = sympy.acos(z / sympy.sqrt(r**2 + z**2))
            r = sympy.sqrt(r**2 + z**2)
        else:
            raise ValueError("Yell at Marco")
        return Coord((r, theta, phi), sys="sph")

    def cyl(self):
        if self.sys == "cyl":
            return self
        elif self.sys == "cart":
            x, y, z = self.coord
            r = sympy.sqrt(x**2 + y**2)
            phi = getPhi(x, y)
        elif self.sys == "sph":
            r, theta, phi = self.coord
            z = r * sympy.cos(theta)
            r = r * sympy.sin(theta)
        else:
            raise ValueError("Yell at Marco")
        return Coord((r, phi, z), sys="cyl")


class Vec:
    def __init__(self, vec, sys=None):
        if sys != "cart" and sys != "sph" and sys != "cyl":
            raise ValueError("Pass sys as 'cart', 'sph', or 'cyl'")
        self.sys = sys
        self.vec = np.array(tuple(sympy.sympify(v) for v in vec), dtype="O")

    def __iter__(self):
        return iter(self.vec)

    def __getitem__(self, index):
        return self.vec[index]

    def __repr__(self):
        return f"Vec({self.vec}, sys={self.sys})"

    def __str__(self):
        return str(self.vec)

    def cart(self, coord):
        if type(coord) is not Coord:
            raise ValueError("Pass a coordinate object")
        _, theta, phi = coord.sph()

        if self.sys == "cart":
            return Vec(self.vec, sys="cart")
        elif self.sys == "sph":
            return Vec(
                np.array(
                    [
                        [
                            sympy.sin(theta) * sympy.cos(phi),
                            sympy.cos(theta) * sympy.cos(phi),
                            -sympy.sin(phi),
                        ],
                        [
                            sympy.sin(theta) * sympy.sin(phi),
                            sympy.cos(theta) * sympy.sin(phi),
                            sympy.cos(phi),
                        ],
                        [sympy.cos(theta), -sympy.sin(theta), 0],
                    ]
                ).dot(self.vec),
                sys="cart",
            )
        elif self.sys == "cyl":
            return Vec(
                np.array(
                    [
                        [sympy.cos(phi), -sympy.sin(phi), 0],
                        [sympy.sin(phi), sympy.cos(phi), 0],
                        [0, 0, 1],
                    ]
                ).dot(self.vec),
                sys="cart",
            )

    def sph(self, coord):
        if type(coord) is not Coord:
            raise ValueError("Pass a coordinate object")
        _, theta, phi = coord.sph()

        if self.sys == "sph":
            return Vec(v, sys="sph")
        elif self.sys == "cart":
            return Vec(
                np.array(
                    [
                        [
                            sympy.sin(theta) * sympy.cos(phi),
                            sympy.sin(theta) * sympy.sin(phi),
                            sympy.cos(theta),
                        ],
                        [
                            sympy.cos(theta) * sympy.cos(phi),
                            sympy.cos(theta) * sympy.sin(phi),
                            -sympy.sin(theta),
                        ],
                        [-sympy.sin(phi), sympy.cos(phi), 0],
                    ]
                ).dot(self.vec),
                sys="sph",
            )
        elif self.sys == "cyl":
            return Vec(
                np.array(
                    [
                        [sympy.sin(theta), 0, sympy.cos(theta)],
                        [sympy.cos(theta), 0, -sympy.sin(theta)],
                        [0, 1, 0],
                    ]
                ).dot(self.vec),
                sys="sph",
            )

    def cyl(self, coord):
        if type(coord) is not Coord:
            raise ValueError("Pass a coordinate object")
        _, theta, phi = coord.sph()

        if self.sys == "cyl":
            return Vec(v, sys="cyl")
        elif self.sys == "cart":
            return Vec(
                np.array(
                    [
                        [sympy.cos(phi), sympy.sin(phi), 0],
                        [-sympy.sin(phi), sympy.cos(phi), 0],
                        [0, 0, 1],
                    ]
                ).dot(self.vec),
                sys="cyl",
            )
        elif self.sys == "sph":
            return Vec(
                np.array(
                    [
                        [sympy.sin(theta), sympy.cos(theta), 0],
                        [0, 0, 1],
                        [sympy.cos(theta), -sympy.sin(theta), 0],
                    ]
                ).dot(self.vec),
                sys="cyl",
            )
