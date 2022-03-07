import numpy as np
import sympy

theta, phi = sympy.symbols("theta phi")
coords = {"cart":("x","y","z"),
        "cyl":("r","phi","z"),
        "sph":("r","theta","phi")}

def parse_eq(eq):
    inf = tuple(arg for arg in eq.args if str(arg)[0] == "d")
    if not inf and str(eq)[0] == "d":
        inf = (eq,)
    parsed = eq.subs(tuple((i,1) for i in inf))
    return (parsed, tuple(str(i)[1:] for i in inf))

def flux(ds,e,bounds):
    sum = 0
    for j,d in enumerate(ds):
        for i,n in enumerate(d):
            a,inf = parse_eq(n)
            if inf:
                constants = {c:bounds[j][c] for c in coords[d.sys] if c not in inf}
                coordinate = Coord(tuple(constants[c] if c in constants else c for c in coords[d.sys]),sys=d.sys)
                sub = tuple((c,getattr(coordinate,e.sys)()[c]) for c in coords[e.sys])
                integrand = a*getattr(e,d.sys)()[i].subs(sub)
                int_bounds = tuple((f,bounds[0][f],bounds[1][f]) for f in inf)
                sum += sympy.integrate(integrand,
                        *int_bounds).subs(tuple((f[0],f[1]) for f in constants.items())).evalf()
    return sum


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
        self.coord = np.array(tuple(map(sympy.sympify, coord)),dtype='O')

    def __iter__(self):
        return iter(self.coord)

    def __getitem__(self, index):
        if index == "x":
            if self.sys == "cart":
                return self.coord[0]
            else:
                raise ValueError("x is only cartesian")

        elif index == "y":
            if self.sys == "cart":
                return self.coord[1]
            else:
                raise ValueError("y is only cartesian")

        elif index == "z":
            if self.sys == "cart" or self.sys == "cyl":
                return self.coord[2]
            else:
                raise ValueError("z is only cartesian or cylindrical")
        elif index == "r":
            if self.sys == "sph" or self.sys == "cyl":
                return self.coord[0]
            else:
                raise ValueError("r is only spherical or cylindrical")
        elif index == "phi":
            if self.sys == "sph":
                return self.coord[2]
            elif self.sys == "cyl":
                return self.coord[1]
            else:
                raise ValueError("phi is only spherical or cylindrical")
        elif index == "theta":
            if self.sys == "sph":
                return self.coord[1]
            else:
                raise ValueError("theta is only spherical")
        else:
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
            phi = sympy.atan2(y, x)
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
            phi = sympy.atan2(y, x)
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

    def cart(self):
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

    def sph(self):
        if self.sys == "sph":
            return Vec(self.vec, sys="sph")
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

    def cyl(self):
        if self.sys == "cyl":
            return Vec(self.vec, sys="cyl")
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
