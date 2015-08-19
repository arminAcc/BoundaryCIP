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

def PlotErrors(l_vel, l_pre) :
    import matplotlib.pyplot as plt
    plt.yscale('log')
    plt.plot(l_vel, "-*", label="vel error")
    plt.plot(l_pre, "-+", label="pre error")
    plt.legend()
    plt.ion()
    plt.show()
    input("after plot - please continue with ENTER")


def AddEdgeFunctions(mesh, fes, sumorder=1, bidx = 1):
    print ("before fes.ndof = ", fes.ndof)
    for sel in mesh.Elements(BND):
        if sel.index == bidx:
            facnr = bcip.GetSElEdge(mesh,sel.nr)
            print ("facnr:",facnr)
            bcip.SetEdgeOrder(fes,facnr, sumorder)
    bcip.UpdateDofTables(fes)
    print ("after fes.ndof = ", fes.ndof)

def stokesCIP(baseorderQ = 1, baseorderV = 2, bonusorderV = 0, boundaryStab = 1.0, cipStab = 1.0, refinements = 0):
    
    geom = SplineGeometry("circleInCircle.in2d")
    mp = MeshingParameters (maxh=0.12)
    ng_mesh = geom.GenerateMesh (mp)
    mesh = Mesh(ng_mesh)
    
    curveorder=2
    mesh.Curve(curveorder)
    
    V = FESpace("h1ho", mesh, order=baseorderV, dirichlet=[2])

    if bonusorderV > 0:
        AddEdgeFunctions(mesh, V, baseorderV+bonusorderV, bidx=1)
        print("V.ndof = ", V.ndof)
        
    Q = FESpace("h1ho", mesh, order=baseorderQ)
    
    X = FESpace([V,V,Q] , flags = {"dgjumps" : True} )
    
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
                      - ( gradu[0] +  gradv[1]) * wp 
                      - (gradwu[0] + gradwv[1]) *  p)
    
    if boundaryStab > 0: 
        bval = DomainConstantCF([0,boundaryStab])
        a.components[2] += BFI (name="BoundLaplace", dim=2, coef=bval) 

    a.components[2] += BFI (name="cip", dim=2, coef=-cipStab) 
    
    
    a.Assemble()
    
    
    #source = ConstantCF(1)
    f = LinearForm(X)
    f.Assemble()
    
    #set global for keeping in memory for gui
    global sol 
    sol = GridFunction(X)
    sol.components[0].Set(VariableCF("5*x"))
    sol.components[1].Set(VariableCF("5*y"))
    
    c = Preconditioner(a, type="direct", flags={ "inverse": "pardiso"} )
    c.Update()
    bvp = BVP(bf=a, lf=f, gf=sol, pre=c, maxsteps=20)
    bvp.Do()
    
    global vel
    vel = e1 * sol.components[0] + e2 * sol.components[1]
    #Draw(mesh=mesh,cf=vel,name="velocity")

    ## Visualize pressure solution and reference pressure
    global ref_pressure
    ref_pressure = VariableCF("(-0.04)") #-0.2*log(sqrt(x*x+y*y)))")
    Draw(mesh=mesh,cf=ref_pressure,name="ref_pressure")
    global pressure
    pressure = sol.components[2]
    err_pre = (ref_pressure - pressure)*(ref_pressure - pressure)
    Draw(pressure)
    # Draw(sol.components[2])
    
    ## Visualize velocity solution and reference velocity
    global ref_vel
    ref_vel= VariableCF("(0.04/(x*x+y*y) * x, 0.04/(x*x+y*y) * y)")
    # Draw(mesh=mesh,cf=ref_vel,name="ref_vel")
    
    global err_vel
    err_vel = (ref_vel - vel)*(ref_vel - vel)


    global l_vel
    global l_pre
    l_vel= []
    l_pre= []
    l_vel.append ( sqrt( Integrate (err_vel, mesh, VOL)) )
    l_pre.append ( sqrt( Integrate (err_pre, mesh, VOL)) )
    print("vel errors:",l_vel)
    print("pre errors:", l_pre)
    
    
    def MyRefine():
        mesh.Refine()
        mesh.Curve(curveorder)
        ## Update:
        X.Update()
        if bonusorderV > 0:
            AddEdgeFunctions(mesh, V, baseorderV+bonusorderV, bidx=1)
            print("V.ndof = ", V.ndof)
        
        sol.Update()
        a.Assemble()
        f.Assemble()
        ## Set Boundary conditions
        sol.components[0].Set(VariableCF("5*x"))
        sol.components[1].Set(VariableCF("5*y"))
        c.Update()
        bvp.Do()
        ## Measure Error
        l_vel.append ( sqrt( Integrate (err_vel, mesh, VOL)) )
        l_pre.append ( sqrt( Integrate (err_pre, mesh, VOL)) )
        print("vel errors:", l_vel)
        print("pre errors:", l_pre)
        Redraw()
    
    for i in range(refinements):
        #sleep(5)
        MyRefine()
        

stokesCIP()

print("you may call\n -stokesCIP(...)\n -PlotErrors(l_vel,l_pre)") 

######### Comments for debugging:
    # blp = BFI (name="BoundLaplace", dim=2, coef=bval)
    # for elidx in range(mesh.GetNE(BND)):
    #     el = ElementId(BND,elidx)
    #     fel = Q.GetFE(el)
    #     print ("el: ", el)
    #     mat = blp.CalcElementMatrix(fel, mesh.GetTrafo(el))
    #     print ("Element matrix of element ", fel, ":\n", mat)
