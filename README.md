## Basically Cheating for ECE219
--
### Requirements:
- [ ] sympy
- [ ] numpy


**Coordinate and Vector transformations**
- Implements all of the neccessary equations on Math Formula 1 sheet

--

### Example usage (in ipython):

```console
In [1]: from matricies import Coord, Vec

In [2]: from numpy import pi

In [3]: coordinate = Coord((10, 4 * pi / 6, pi / 6), sys="sph")

In [4]: vector = Vec(("sin(theta)**2*cos(phi)", "cos(phi)**2", "-sin(phi)"), sys
   ...: ="sph")

In [5]: vector.cyl(coordinate)
Out[5]: Vec([0.187500000000000 -0.500000000000000 -0.974278579257494], sys=cyl)

```
