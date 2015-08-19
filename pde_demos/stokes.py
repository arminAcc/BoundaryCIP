from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.la import *
from ngsolve.solve import *
from ngsolve.utils import *

from numpy import linspace
from timeit import Timer
from time import sleep
from math import pi

from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from libbcip import *

geom = SplineGeometry("circleInCircle.in2d")
mp = MeshingParameters (maxh=0.08)
mesh = Mesh(geom.GenerateMesh (mp))

curveorder=2
mesh.Curve(curveorder)
#mesh.Refine()

Vx = FESpace("h1ho", mesh, order=2, dirichlet=[2]) #, flags = {"dgjumps" : True} )
Vy = FESpace("h1ho", mesh, order=2, dirichlet=[2]) #, flags = {"dgjumps" : True} )

Q = FESpace("h1ho", mesh, order=1) #, flags = {"dgjumps" : True})

X = FESpace([Vx,Vy,Q] , flags = {"dgjumps" : True} )

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

bval = DomainConstantCF([0,1])
a.components[2] += BFI (name="BoundLaplace", dim=2, coef=bval) 
a.components[2] += BFI (name="cip", dim=2, coef=-1) 

# blp = BFI (name="BoundLaplace", dim=2, coef=bval)
# for elidx in range(mesh.GetNE(BND)):
#     el = ElementId(BND,elidx)
#     fel = Q.GetFE(el)
#     print ("el: ", el)
#     mat = blp.CalcElementMatrix(fel, mesh.GetTrafo(el))
#     print ("Element matrix of element ", fel, ":\n", mat)

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
#Draw(mesh=mesh,cf=vel,name="velocity")
Draw(sol.components[2])


# def AddStab(param=1.0):
#     bval = DomainConstantCF([0,param])
#     a.components[2] += BFI (name="BoundLaplace", dim=2, coef=bval) 
#     X.Update()
#     sol.Update()
#     a.Assemble()
#     f.Assemble()
#     sol.components[0].Set(VariableCF("x"))
#     sol.components[1].Set(VariableCF("y"))
#     c.Update()
#     bvp.Do()
#     Redraw()

l= []

ref_vel = VariableCF("(0.04/(x*x+y*y) * x, 0.04/(x*x+y*y) * y)")
Draw(mesh=mesh,cf=ref_vel,name="ref_vel")

ref_pressure = VariableCF("(0.04/(x*x+y*y))")
pressure = sol.components[2]
err = (ref_pressure - pressure)*(ref_pressure - pressure)

l.append ( Integrate (err, mesh, VOL) )
print(l)

Draw(mesh=mesh,cf=ref_vel,name="ref_vel")

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
    l.append ( Integrate (err, mesh, VOL) )
    print(l)
    Redraw()

for i in range(2):
    sleep(5)
    MyRefine()
    
# Draw(sol)

# print (a.mat)
# print (f.vec)


