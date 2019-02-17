import math
from ufl import *

import dune.create as create
import dune.fem as fem
#from dune.grid import cartesianDomain
from dune.ufl import Space



def identity(x):
	return x


def Ball(x):
	y=[0, 0, 0]
	
	r1 = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
	r0 = 1.
		
	y[0] = r0*x[ 0 ]/r1
	y[1] = r0*x[ 1 ]/r1
	y[2] = r0*x[ 2 ]/r1
	return y

def Bloodcell(x):
	y=[0, 0, 0]
	r1 = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
	r0 = 4.
	y[0] = r0*x[ 0 ]/r1
	y[1] = r0*x[ 1 ]/r1	
	
	r2 = sqrt(x[0]*x[0]+x[1]*x[1])
	r3 = sqrt(y[0]*y[0]+y[1]*y[1])
	y[2] = 1.5-0.5*cos(pi*r3/2.)
	if r3 > 2.:
		y[2]=sqrt(4.-(r3-2.)*(r3-2.))
	if x[2] < 0.:
		y[2]=-y[2]	
	return y

def Pear(x):
	y=[0, 0, 0]
	y[0] = r0*x[0]/r1
	y[1] = r0*x[1]/r1
	if y[2]>0:
		if r3<2.5:
			y[2]= sqrt(r0*r0-2.5*2.5)+2-2*cos(pi*pow( (2.5-r3)/2.5,1.5));
			if r3 <1.35:
				y[2]=sqrt(r0*r0-2.5*2.5)+2-2*cos(pi*pow( (1.15)/2.5,1.5))+0.5*pow(1.35*1.35-r3*r3, 1./1.25);
	return y
