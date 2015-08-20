from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.la import *
from ngsolve.solve import *
from ngsolve.utils import *

from numpy import linspace
from timeit import Timer
from time import sleep
from math import pi
from math import sqrt

from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from libbcip import *

geom = SplineGeometry("circleInCircle.in2d")
mp = MeshingParameters (maxh=0.12)
ng_mesh = geom.GenerateMesh (mp)
mesh = Mesh(ng_mesh)

coef_idx = bcip.idx()
coef_h = bcip.h()
coef_nr = bcip.nr()

vtk = VTKOutput(mesh, coefs=[coef_idx], names=["coef_idx"], subdivision=0)
vtk.Do()

Draw(coef_idx, mesh, "idx")
Draw(coef_h, mesh, "h")
Draw(coef_nr, mesh, "nr")

curveorder=2
mesh.Curve(curveorder)

