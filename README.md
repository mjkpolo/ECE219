# Basically Cheating for ECE219

### Setup:

run this to install mamba-forge:
```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```
after mamba-forge is installed, run:
```
mamba install sympy numpy ipython
```
finally, navigate to this folder and run `ipython` and follow the example usage below

`matricies.py`: Coordinate and Vector transformations
----

Coordinates and Vectors can easily be converted with `.sph(), .cyl(), .cart()`

use `evaluate_vector(vector, coordinate)` to plug in values from the coordinate. r depends on the system of the coordinate passed, but x,y,z,theta,phi are all independent of coordinate system


### Example usage in iPython:

```python
In [1]: from matricies import Coord, Vec, evaluate_vector

In [2]: c = Coord((10, 4 * pi / 6, pi / 6), sys="sph")
---------------------------------------------------------------------------
NameError                                 Traceback (most recent call last)
Input In [2], in <module>
----> 1 c = Coord((10, 4 * pi / 6, pi / 6), sys="sph")

NameError: name 'pi' is not defined

In [3]: from sympy import pi

In [4]: c = Coord((10, 4 * pi / 6, pi / 6), sys="sph")

In [5]: c.sph()
Out[5]: Coord([10 2*pi/3 pi/6], sys=sph)

In [6]: c.cyl()
Out[6]: Coord([5*sqrt(3) pi/6 -5], sys=cyl)

In [7]: c.cart()
Out[7]: Coord([15/2 5*sqrt(3)/2 -5], sys=cart)

In [8]: v = Vec(("sin(theta)**2*cos(phi)", "cos(phi)**2", "-sin(phi)"), sys="sph")

In [9]: v.sph()
Out[9]: Vec([sin(theta)**2*cos(phi) cos(phi)**2 -sin(phi)], sys=sph)

In [10]: v.cyl()
Out[10]:
Vec([sin(theta)**3*cos(phi) + cos(phi)**2*cos(theta) -sin(phi)
 sin(theta)**2*cos(phi)*cos(theta) - sin(theta)*cos(phi)**2], sys=cyl)

In [11]: v.cart()
Out[11]:
Vec([sin(phi)**2 + sin(theta)**3*cos(phi)**2 + cos(phi)**3*cos(theta)
 sin(phi)*sin(theta)**3*cos(phi) + sin(phi)*cos(phi)**2*cos(theta) - sin(phi)*cos(phi)
 sin(theta)**2*cos(phi)*cos(theta) - sin(theta)*cos(phi)**2], sys=cart)

In [12]: evaluate_vector(v.cyl(), c)
Out[12]: (0.1875, -0.5, -0.9742785792574935)

```

### Example without script:

Calculate the z component of E⃗ at at R⃗ (x,y,z)=(0,0,2) [meters] due to a circular disk of radius r=6 [meters] and uniform density ρs=3 [C/m2] in the z=0 plane

 ```python
In [2]: from sympy.abc import *

In [3]: from sympy import *

In [4]: bounds = ((phi, 0, 2 * pi), (r, 0, 6))

In [5]: diff = (-r * cos(phi), -r * sin(phi), 2)

In [6]: mag = sqrt(diff[0] ** 2 + diff[1] ** 2 + diff[2]
   ...: ** 2) ** 3

In [7]: p = 3

In [8]: integrate(p * r / (4 * pi * 8.85e-12 * mag) * diff[2], *bounds).evalf()
Out[8]: 115893598980.197
```
