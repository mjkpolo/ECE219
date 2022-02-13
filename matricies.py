import numpy as np
from sympy import sympify

def get_phi(x=None,y=None):
    if not x or not y:
        raise ValueError("pass x and y")
    return np.arctan(y/x) + np.pi if x<0 else np.arctan(y/x)

def get_theta(x=None,y=None,z=None,r=None):
    if (not x or not y or not z) and (not r or not z):
        raise ValueError("pass x and y and z, or r and z")
    if not x and not y:
        return np.arccos(z/np.sqrt(r**2+z**2))
    else:
        return np.arccos(z/np.sqrt(x**2+y**2+z**2))

def evaluate_vector(vec, coord, sys=None):
    if sys != 'cart' and sys != 'sph' and sys != 'cyl':
        raise ValueError("Pass sys as 'cart', 'sph', or 'cyl'")
    v = np.array(vec,dtype='O')

    for i in range(v.size):
        v[i] = sympify(v[i])

    if sys == 'cart':
        c = coord.cart()
        for i in range(v.size):
            v[i] = v[i].subs('x',c[0])
            v[i] = v[i].subs('y',c[1])
            v[i] = v[i].subs('z',c[2])
    elif sys == 'cyl':
        c = coord.cyl()
        for i in range(v.size):
            v[i] = v[i].subs('r',c[0])
            v[i] = v[i].subs('phi',c[1])
            v[i] = v[i].subs('z',c[2])
    elif sys == 'sph':
        c = coord.sph()
        for i in range(v.size):
            v[i] = v[i].subs('r',c[0])
            v[i] = v[i].subs('theta',c[1])
            v[i] = v[i].subs('phi',c[2])
    return v


class Coord:
    def __init__(self,coord,sys=None):
        if sys != 'cart' and sys != 'sph' and sys != 'cyl':
            raise ValueError("Pass sys as 'cart', 'sph', or 'cyl'")
        self.sys=sys
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
        if self.sys == 'cart':
            return self
        elif self.sys == 'sph':
            r,theta,phi = self.coord

            x = r*np.sin(theta)*np.cos(phi)
            y = r*np.sin(theta)*np.sin(phi)
            z = r*np.cos(theta)
        elif self.sys == 'cyl':
            r,phi,z = self.coord

            x = r*np.cos(phi)
            y = r*np.sin(phi)
        else:
            raise ValueError("This shouldn't happen ever...")
        return Coord((x,y,z), sys='cart')

    def sph(self):
        if self.sys == 'sph':
            return self
        elif self.sys == 'cart':
            x,y,z = self.coord

            theta = get_theta(x=x,y=y,z=z)
            phi = get_phi(x=x,y=y)
            r = np.sqrt(x**2 + y**2 + z**2)
        elif self.sys == 'cyl':
            r,phi,z = self.coord

            theta = get_theta(z=z,r=r)
            r = np.sqrt(r**2+z**2)
        else:
            raise ValueError("This shouldn't happen ever...")
        return Coord((r,theta,phi), sys='sph')

    def cyl(self):
        if self.sys == 'cyl':
            return self
        elif self.sys == 'cart':
            x,y,z = self.coord
            
            r = np.sqrt(x**2+y**2)
            phi = get_phi(x=x,y=y)
        elif self.sys == 'sph':
            r,theta,phi = self.coord

            z = r*np.cos(theta)
            r = r*np.sin(theta)
        else:
            raise ValueError("This shouldn't happen ever...")
        return Coord((r,phi,z), sys='cyl')


class Vec:
    def __init__(self,vec,sys=None):
        if sys != 'cart' and sys != 'sph' and sys != 'cyl':
            raise ValueError("Pass sys as 'cart', 'sph', or 'cyl'")
        self.sys=sys
        self.vec = vec

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
        v = evaluate_vector(self.vec, coord, sys=self.sys)
        if self.sys == 'cart':
            return v
        elif self.sys == 'sph':
            _,theta,phi = coord.sph()
            return Vec(
                    np.array(
                    [[np.sin(theta)*np.cos(phi), np.cos(theta)*np.cos(phi), -np.sin(phi)],
                    [np.sin(theta)*np.sin(phi), np.cos(theta)*np.sin(phi), np.cos(phi)],
                    [np.cos(theta), -np.sin(theta), 0]]
                    ).dot(v),
                    sys='cart')
        elif self.sys == 'cyl':
            _,phi,_ = coord.cyl()
            return Vec(
                    np.array(
                    [[np.cos(phi),-np.sin(phi),0],
                    [np.sin(phi), np.cos(phi), 0],
                    [0,0,1]]
                    ).dot(v),
                    sys='cart')

    def sph(self,coord):
        if type(coord) is not Coord:
            raise ValueError("Pass a coordinate object")
        v = evaluate_vector(self.vec, coord, sys=self.sys)
        if self.sys == 'sph':
            return v
        elif self.sys == 'cart':
            _,theta,phi = coord.sph()
            return Vec(
                    np.array(
                    [[np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)],
                    [np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)],
                    [-np.sin(phi), np.cos(phi), 0]]
                    ).dot(v),
                    sys='sph')
        elif self.sys == 'cyl':
            _,theta,_ = coord.sph()
            return Vec(
                    np.array(
                    [[np.sin(theta),0,np.cos(theta)],
                    [np.cos(theta), 0, -np.sin(theta)],
                    [0,1,0]]).dot(v),
                    sys='sph')

    def cyl(self,coord):
        v = evaluate_vector(self.vec, coord, sys=self.sys)
        if type(coord) is not Coord:
            raise ValueError("Pass a coordinate object")
        if self.sys == 'cyl':
            return self
        elif self.sys == 'cart':
            _,phi,_ = coord.cyl()
            return Vec(
                    np.array(
                    [[np.cos(phi), np.sin(phi), 0],
                    [-np.sin(phi), np.cos(phi), 0],
                    [0, 0, 1]]).dot(v),
                    sys='cyl')
        elif self.sys == 'sph':
            _,theta,_ = coord.sph()
            return Vec(
                    np.array(
                    [[np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1],
                    [np.cos(theta), -np.sin(theta), 0]]).dot(v),
                    sys='cyl')
