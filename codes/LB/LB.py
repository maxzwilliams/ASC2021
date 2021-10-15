"""
Two dimensional lattice boltzman code for a convecting fluid

We will use the D2Q9 lattice for this.
"""

import numpy as np

import math




def main():

	nx, ny = 10**3, 10**3
	nt = 10**3
	
	
	
	cs = np.array[[0,0], [0,-1], [0,1], [-1,0],[-1,-1], [-1,1], [1,0], [1,-1], [1,1]   ]
	
	opposites = np.array([0, 6, 3, 2, 8, 7, 1, 5, 4])
	## index i in cs has opposite that is index opposites[i] in cs
	
	## properties of the D2Q9 lattice
	D = 2
	Q = 9
	
	## lattice weights
	ws = []
	for el in cs:
		if (abs(el[0]) == 1 and abs(el[1]) == 1):
			ws.append(1/36)
		elif (abs(el[0]) == 1 or abs(el[1]) == 1):
			ws.append(1/9)
		else:
			ws.append(4/9)
	ws = np.array(ws)
	
	c1 = 1.0
	c2 = 3.0/(ls**2)
	c3 = 9.0/(2.0*Ls**4)
	c4 = -3.0/(2.0*Ls**2)
	
	visc = 0.1
	vf = 3 
	vg = 0.6
	tau_f = vf/(Ls*dt) + 0.5
	tau_g = vg/(Ls*dt) + 0.5

			
	
	dt = 1
	dx = 1
	Ls = dt/dx ## the lattice speed in all directions
	
	## fields
	
	f = np.zeros((Q, nx, ny))
	f_s = np.zeros((Q, nx, ny))
	f_eq = np.zeros((Q, nx, ny))
	
	rho = np.zeros((D, nx, ny))
	u = np.zeros((D, nx, ny))
	
	
	
	
	
	
	
	
	
	
	
	

	

