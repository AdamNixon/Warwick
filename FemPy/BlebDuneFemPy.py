from __future__ import print_function

import math
from ufl import *

import dune.create as create
import dune.fem as fem
from dune.fem.function import integrate
#from dune.grid import cartesianDomain
from dune.ufl import Space, NamedConstant
from deformation import identity, Ball, Bloodcell, Pear
from functions import ComputeArea, A1, A2,  norm, InvNormNxN, Nu, metric, cortexCalc, cortexDistance

########################################################
#########################################################

gridSelect = "cell.dgf"
# Select how you wish to deform your grid (dgf file)
case = 0		# number selects which surfaceOption will be used
Graphics = 1	# 0: outputs (u,w) over grid, 1: outputs  u over the grid warped by u

surfaceOptions = { 0 : lambda x: identity(x), # Full description of functions can be found in deformation.py
				   1 : lambda x: Ball(x),
				   2 : lambda x: Bloodcell(x),
				   3 : lambda x: Pear(x)
				 }
# Variables
deltaT = 0.002
finalT=  2
K_PSI = 1
K_B = 0.005
X_0 = 0.95
P_0 = 7.5
K_0 =  4
U_B = 0.056
L_0 = 0.04
Omega = 1


############################################
############################################


# polynomial order of surface approximation
order = 1

fem.parameter.append({"fem.verboserank": 0,
			   "fem.solver.verbose": 0,
                           "istl.preconditioning.method": "ilu",
                           "istl.preconditioning.iterations": 1,
                           "istl.preconditioning.relaxation": 1})


    # set up grid 

grid = create.grid("ALUSimplex", gridSelect, dimgrid=2, dimworld=3)
grid.hierarchicalGrid.globalRefine(2)
print("Number of elements:",grid.size(0))
print("Number of vertices:",grid.size(grid.dimension))

    # set up a lagrange scalar space with polynomial order 2 over that grid
spc = create.space("Lagrange", grid, dimrange=grid.dimWorld, order=order, store="istl")

    # non spherical initial surgface

uDimWorld=6
positions = spc.interpolate(lambda x: surfaceOptions[case](x) , name="position")
surface   = create.view("geometry",positions)
SolutionSpace = create.space("Lagrange", surface, dimrange=uDimWorld, order=order, storage="istl")

if Graphics==1:
	positionsG = spc.interpolate(lambda x: surfaceOptions[case](x) , name="position")
	surfaceG   = create.view("geometry",positionsG)
	GraphicalSpace =  create.space("Lagrange", surfaceG, dimrange=3, order=order)

totalArea= ComputeArea(surface)
print ( totalArea)

    # set up initial conditions
solution = SolutionSpace.interpolate(lambda x: [x[0], x[1], x[2], 0, 0, 0 ], name="solution")
if Graphics==1:
	solutionGraphical = GraphicalSpace.interpolate( as_vector([solution[0], solution[1], solution[2]]) , name="solution")
	surfaceG.writeVTK("Graphical", pointdata=[solution], number=0)
elif Graphics==0:
	surface.writeVTK("heat", pointdata=[solution], number=0)

	# get a discrete function to hold the old solution and tell the model to use that for the coefficient u_n
solution_n = solution.copy();
solution_0 = solution.copy();

    # now define the actual pde to solve:
cell = SolutionSpace.cell()
u = TrialFunction(SolutionSpace)
v = TestFunction(SolutionSpace)
x = SpatialCoordinate(cell)
u_n = Coefficient(SolutionSpace)
u_0 = Coefficient(SolutionSpace)
tau = NamedConstant(cell, "tau")
x_0 = NamedConstant(cell, "x_0")
omega = NamedConstant(cell, "omega")
m_psi = NamedConstant(cell, "k_psi")
k_b = NamedConstant(cell, "k_b")
k_0 = NamedConstant(cell, "k_0")
l_0 = NamedConstant(cell, "l_0")
u_B = NamedConstant(cell, "u_B")
p_0 = NamedConstant(cell, "p_0")
IntVolume = NamedConstant(cell, "IntVolume")


###################
# Bilinear form terms
###################

# for full description of functions called see functions.py
Bending = k_b*inner ( grad(A2(u)) ,  grad(A1(v))  )
f = inner ( grad(A1(u)) ,  grad(A2(v))  )
g = inner( A2(u), A2(v) )
Tension_im =  m_psi*inner( grad(A1(u)), grad(A1(v)) ) 
Tension_ex =  m_psi*x_0*sqrt(inner(grad(A1(u_0)), grad(A1(u_0)))) *InvNormNxN( u_n )*inner( grad(A1(u_n)) , grad(A1(v)) ) 
Pressure = (p_0/(IntVolume+1e-10))*inner( Nu(u_n, u_n) , A1(v) )
Connection = conditional( norm( cortexDistance(u_n, u_0, l_0) )<u_B,  1 ,  0 ) 
Compression = conditional( norm( cortexDistance(u_n, u_0, l_0) )<0.02,  101 ,  1 ) 
NuC =(1/norm(cortexDistance(u_n, u_0, l_0)))*(A1(cortexDistance(u_n, u_0, l_0))) 
Linkers_im = k_0 * inner(Connection*Compression* A1( u) , A1(v) )
Linkers_ex = k_0 * inner(Connection*Compression* (A1(cortexCalc(u_n, u_0, l_0))+l_0*NuC), A1(v) )

# Model equation

a_ex = ( omega*inner( A1(u_n), A1(v) ) +tau*(Tension_ex+Pressure+Linkers_ex) )*dx
a_im = ( omega*inner( A1(u), A1(v) )   +tau*(Tension_im+Bending+Linkers_im) +g-f )*dx

equation = a_im == a_ex


# now generate the model code and compile
model = create.model("elliptic", surface, equation, coefficients={u_n:solution_n,u_0:solution_0})

# Create volume
Volume= (1/3)* inner( Nu(solution ,solution_0) , A1(solution) ) 
intVolume = integrate(surface, Volume, order=1)[0]


    # create the solver using a standard fem scheme
scheme = create.scheme("h1", SolutionSpace, model,  solver="bicgstab")
scheme.model.tau = deltaT
scheme.model.x_0 = X_0
scheme.model.omega = Omega
scheme.model.k_psi = K_PSI
scheme.model.k_b = K_B
scheme.model.p_0 = P_0
scheme.model.IntVolume = intVolume
scheme.model.k_0 = K_0
scheme.model.l_0 = L_0
scheme.model.u_B = U_B

    # now loop through time and output the solution after each time step
steps = int(finalT/ deltaT)
print("Begin time step" )
for n in range(1,steps+1):  
	solution_n.assign(solution)
	scheme.solve(target=solution) 
	if Graphics==1:
		solutionGraphical = GraphicalSpace.interpolate( as_vector([solution[0], solution[1], solution[2]]),  name="solution" )
		positionsG.dofVector.assign(solutionGraphical.dofVector)
		totalAreaNew = ComputeArea( surfaceG)
		print("Area: ", totalAreaNew )
		surfaceG.writeVTK("Graphical", pointdata=[solutionGraphical], number=n)
	elif Graphics==0:
		surface.writeVTK("heat", pointdata=[solution], number=n)
	print(n*deltaT)
Volume= (1/3)* inner( Nu(solution ,solution_0) , A1(solution) ) 
intVolume = integrate(surface, Volume, order=1)[0]
#scheme.model.IntVolume = intVolume
print("Volume: ",  intVolume )
R8=sqrt(inner(A1(solution),A1(solution)))
R8I=integrate(surface, R8, order=1)[0]
print("Radius: ",  R8I/totalArea )
