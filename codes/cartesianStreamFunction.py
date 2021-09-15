"""
Simple program for simulating a fluid in a box.

This code was written as an exercise towards getting a convecting fluid in a box. The top plate of the box is moving to drive flow.	
"""
import matplotlib.pyplot as plt
import time
import copy
import math

visc = 0.1
dt = 0.005

nx=17
ny=17

sf = [[0 for x in range(nx)] for y in range(ny)]
vt = [[0 for x in range(nx)] for y in range(ny)]
w =  [[0 for x in range(nx)] for y in range(ny)]

h=1/(nx-1) ## we assume that nx and ny are the same here

timeSteps = math.floor(1/dt * 1.2)

errorStandard = 0.001
Beta=0.5
t = 0

solverIterations = 1000
for timeStep in range(timeSteps):
	
	for solverIteration in range(solverIterations):
		w = copy.deepcopy(sf)
		
		
		oldsf = copy.deepcopy(sf)
		for i in range(1, nx-1):
			for j in range(1, ny-1):
				##sf[i][j] = 0.25*Beta*(sf[i+1][j] + sf[i-1][j]+sf[i][j+1] + sf[i][j-1] + h**2 *vt[i][j]) + (1-Beta)*sf[i][j]
				sf[i][j]=0.25*Beta*(oldsf[i+1][j]+oldsf[i-1][j]+oldsf[i][j+1]+oldsf[i][j-1]+h*h*vt[i][j])+ (1.0-Beta)*oldsf[i][j]
				
		totalError = 0
		for i in range(nx):
			for j in range(ny):
				totalError = totalError + abs(w[i][j] - sf[i][j])
		if (totalError < errorStandard):
			print("converged")
			break;
			
	## here we update the vorticity on the boundries
	
	for i in range(1, nx-1):
		vt[i][0] = -2*sf[i][1]/(h**2)
		
		vt[i][ny-1] = -2*sf[i][ny-2]/(h**2) - 2/h
		
		vt[0][i] = -2*sf[1][i]/(h**2)
		
		vt[nx-1][i]=  -2*sf[nx-2][i]/(h**2) 
	
	
	for i in range(1, nx-1):
		for j in range(1, ny-1):
			w[i][j]=-0.25*((sf[i][j+1]-sf[i][j-1])*(vt[i+1][j]-vt[i-1][j])-(sf[i+1][j]-sf[i-1][j])*(vt[i][j+1]-vt[i][j-1]))/(h*h)+visc*(vt[i+1][j]+vt[i-1][j]+vt[i][j+1]+vt[i][j-1]-4.0*vt[i][j])/(h*h)

	for i in range(1,nx-1):
		for j in range(1, ny-1):
			vt[i][j] = vt[i][j] + dt*w[i][j]
			
	t += dt

	
	plt.contour(vt)
	plt.savefig("vts//file"+str(timeStep)+".png")
	plt.clf()
	
	plt.contour(sf)
	plt.savefig("sfs//file"+str(timeStep)+".png")
	plt.clf()
	
	
	
	
	
	## then here we 
				


