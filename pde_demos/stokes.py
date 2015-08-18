from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.la import *
from ngsolve.solve import *
from ngsolve.utils import *

from numpy import linspace
from timeit import Timer
from math import pi

from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from libbcip import *

geom = SplineGeometry("circleInCircle.in2d")
mp = MeshingParameters (maxh=0.08)
mesh = Mesh(geom.GenerateMesh (mp))

curveorder=3
mesh.Curve(curveorder)
#mesh.Refine()

Vx = FESpace("h1ho", mesh, order=2, dirichlet=[2])
Vy = FESpace("h1ho", mesh, order=2, dirichlet=[2])

Q = FESpace("h1ho", mesh, order=1)

X = FESpace([Vx,Vy,Q])

u,v,p = X.TrialFunction()
wu,wv,wp = X.TestFunction()

gradu = u.Deriv()
gradv = v.Deriv()
gradwu = wu.Deriv()
gradwv = wv.Deriv()

e1 = VariableCF("(1,0)")
e2 = VariableCF("(0,1)")
# e1 = CoefficientFunction((1,0))
# e2 = CoefficientFunction((0,1))

a = BilinearForm(X, symmetric=True)
a += SymbolicBFI( gradu*gradwu + gradv*gradwv + u*wu + v*wv 
                  - (gradu[0] +  gradv[1]) * wp 
                  - (gradwu[0] + gradwv[1]) * p)

# bval = DomainConstantCF([0,1])

# a.components[2] += BFI (name="BoundLaplace", dim=2, coef=bval) 
a.Assemble()

#source = ConstantCF(1)
f = LinearForm(X)
f.Assemble()

sol = GridFunction(X)

sol.components[0].Set(VariableCF("x"))
sol.components[1].Set(VariableCF("y"))

c = Preconditioner(a, type="direct", flags={ "inverse": "pardiso"} )
c.Update()
bvp = BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=20)
bvp.Do()

vel = e1 * sol.components[0] + e2 * sol.components[1]
Draw(mesh=mesh,cf=vel,name="velocity")


def AddStab(param=1.0):
    bval = DomainConstantCF([0,param])
    a.components[2] += BFI (name="BoundLaplace", dim=2, coef=bval) 
    X.Update()
    sol.Update()
    a.Assemble()
    f.Assemble()
    sol.components[0].Set(VariableCF("x"))
    sol.components[1].Set(VariableCF("y"))
    c.Update()
    bvp.Do()
    Redraw()


def MyRefine():
    mesh.Refine()
    mesh.Curve(curveorder)
    X.Update()
    sol.Update()
    a.Assemble()
    f.Assemble()
    sol.components[0].Set(VariableCF("x"))
    sol.components[1].Set(VariableCF("y"))
    c.Update()
    bvp.Do()
    Redraw()

#Draw(sol)

# print (a.mat)
# print (f.vec)


