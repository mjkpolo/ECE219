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

Coordinates can be converted, and Vectors can be transformed with a coordinate parameter to rotate around.

use `evaluate_vector(vector, coordinate)` to plug in values from the coordinate. r depends on the system of the coordinate passed


### Example usage in iPython:

```python
In [20]: c = Coord((10, 4 * pi / 6, pi / 6), sys="sph")

In [21]: c.cyl()
Out[21]: Coord([5*sqrt(3) pi/6 -5], sys=cyl)

In [22]: c.cart()
Out[22]: Coord([15/2 5*sqrt(3)/2 -5], sys=cart)

In [23]: v = Vec(("sin(theta)**2*cos(phi)", "cos(phi)**2", "-sin(phi)"), sys="sph")

In [24]: v.cyl(c)
Out[24]:
Vec([sqrt(3)*sin(theta)**2*cos(phi)/2 - cos(phi)**2/2 -sin(phi)
 -sin(theta)**2*cos(phi)/2 - sqrt(3)*cos(phi)**2/2], sys=cyl)

In [25]: evaluate_vector(v.cyl(c), c)
Out[25]: (0.1875, -0.5, -0.9742785792574935)

```

