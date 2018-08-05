import math
from ufl import *

import dune.create as create
import dune.fem as fem
#from dune.grid import cartesianDomain
from dune.ufl import Space


def ComputeArea(surface):
	totalArea = 0
	for element in surface.elements:
		totalArea += element.geometry.volume
	return totalArea

def A1(u):
    return as_vector([u[0], u[1], u[2]])

def A2(u):
    return as_vector([ u[3], u[4], u[5] ])

def norm(u):
    return sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])

def InvNormNxN( eta ):
    S1=grad(A1(eta)[0] )
    S2=grad(A1(eta)[1] )
    S3=grad(A1(eta)[2] )
    S=0
    for i in range(0,3):
    	S+=S1[i]*S1[i]+S2[i]*S2[i]+S3[i]*S3[i]
    return 1/sqrt(S)

def Nu(eta, psi):  #used in pressure
	w1=[grad(A1(eta)[0] )[0],  grad(A1(eta)[1] )[0], grad(A1(eta)[2] )[0] ]
	w2=[grad(A1(eta)[0] )[1],  grad(A1(eta)[1] )[1], grad(A1(eta)[2] )[1] ]
	Nux=w1[1]*w2[2]-w1[2]*w2[1]
	Nuy=-(w1[0]*w2[2]-w1[2]*w2[0])
	Nuz=w1[0]*w2[1]-w1[1]*w2[0]

	Nu= [Nux, Nuy, Nuz ]
	norm=0
	for i in range(0, 3):
		norm+=Nu[i]*Nu[i]
	
	for i in range(0, 3):
		Nu[i]/=sqrt(norm)
	direction=conditional( Nux*psi[0] + Nuy*psi[1] + Nuz*psi[2]>0  ,  1 ,  -1 ) 
	return direction*as_vector([Nu[0], Nu[1], Nu[2]]) # note cross product returns inward normal

def metric(eta):  #used in pressure
    w1=grad(A1(eta)[0] )
    w2=grad(A1(eta)[1] )
    w3=grad(A1(eta)[2] )
    a11=0
    a22=0
    a33=0
    a12=0
    a13=0
    a23=0
    for i in range(0, 3):   
        a11+=w1[i]*w1[i]
        a22+=w2[i]*w2[i]
        a33+=w3[i]*w3[i]
        a12+=w1[i]*w2[i]
        a13+=w1[i]*w3[i]
        a23+=w2[i]*w3[i]
    mu=sqrt(a11*a11+a22*a22+a33*a33+2*(a12*a12+a13*a13+a23*a23))
    return mu


def cortexCalc(u, u_0, L_0): # defines the cortex
	y=[0,0,0,0,0,0]
	w1=[grad(A1(u_0)[0] )[0],  grad(A1(u_0)[1] )[0], grad(A1(u_0)[2] )[0] ]
	w2=[grad(A1(u_0)[0] )[1],  grad(A1(u_0)[1] )[1], grad(A1(u_0)[2] )[1] ]
	Nux=w1[1]*w2[2]-w1[2]*w2[1]
	Nuy=-(w1[0]*w2[2]-w1[2]*w2[0])
	Nuz=w1[0]*w2[1]-w1[1]*w2[0]
	Nu= [Nux, Nuy, Nuz ]
	norm=0
	for i in range(0, 3):
		norm+=Nu[i]*Nu[i]
	
	for i in range(0, 3):
		Nu[i]/=sqrt(norm)
	
	for i in range(0, 3):
		y[i]=u_0[i]-L_0*Nu[i]
	return y

def cortexDistance(u, u_0, L_0):	# u-u_c
	y=cortexCalc(u, u_0, L_0)
	for i in range(0, 3):
		y[i]=u[i]-y[i]
	return y
