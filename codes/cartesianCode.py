

import matplotlib.pyplot as plt
import time
import copy
import math
import random
import numpy as np


## im sure this works now
def updateSf(vt, sf, T, q, ss):
	for iteration in range(ss['poissonIterations']):
		ogsf =  np.copy(sf)
		
		sfij1 = np.roll(sf, -1, axis=1)
		sfijm1 = np.roll(sf, 1, axis=1)
		sfi1j = np.roll(sf, -1, axis=0)
		sfim1j = np.roll(sf, 1, axis=0)
		
		
		term1 = vt + (sfi1j-sfim1j)/(ss['dx']**2) + (sfij1-sfijm1)/(ss['dy']**2)
	
		sf = ss['SORparam'] * ((ss['dx']*ss['dy'])**2)/(2*((ss['dx'])**2 + (ss['dy'])**2)) * term1 + (1- ss['SORparam']) * sf
		
		for xIndex in range(ss['xSteps']):
			sf[0][xIndex] = 0
			sf[ ss['ySteps']+1 ][xIndex] = 0
		
		ep = 0
		for xIndex in range(ss['xSteps']):
			for yIndex in range(ss['ySteps']):
				ep += abs(ogsf[yIndex][xIndex] - sf[yIndex][xIndex])
				
		if (ep < ss['poissonError']):
			break;
		if (iteration == ss['poissonIterations'] - 1):
			raise Exception("did not converge")
			
	return sf
	
	
def updateT(vt, sf, T, q, ss):
	Tij1 = np.roll(T, -1, axis=1)
	Tijm1 = np.roll(T, 1, axis=1)
	Ti1j = np.roll(T, -1, axis=0)
	Tim1j = np.roll(T, 1, axis=0)
	
	sfij1 = np.roll(sf, -1, axis=1)
	sfijm1 = np.roll(sf, 1, axis=1)
	sfi1j = np.roll(sf, -1, axis=0)
	sfim1j = np.roll(sf, 1, axis=0)
	
	
	
	u = -(sfi1j - sfim1j)/(2*ss['dy'])
	v = (sfij1 - sfijm1)/(2*ss['dx'])
	
	umax = np.amax(np.absolute(u))
	vmax = np.amax(np.absolute(v))
	
	if ( umax*ss['dt']/ss['dx'] >= 1 or vmax * ss['dt']/ss['dy'] >= 1):
		return vt, "time error"
	
	
	advectionx = 1/(ss['dx']) *( np.multiply(np.absolute(u), 0.5*(Tij1-Tijm1)) - np.multiply(u, (0.5*Tij1 - T +0.5*Tijm1))  )
	advectiony = 1/(ss['dy']) *( np.multiply(np.absolute(v), 0.5*(Ti1j-Tim1j)) - np.multiply(v, (0.5*Ti1j - T + 0.5*Tim1j)))
	
	diffusionTerm = ss['kappa'] * ( ( Tij1-2*T+Tijm1 )/(ss['dx']**2) +  (Ti1j-2*T+Tim1j )/(ss['dy']**2))
	source = q/(ss['rho0'] * ss['cv'])
	
	T = ss['dt'] * (- advectionx - advectiony + diffusionTerm + source) + T
	
	
	for xIndex in range(ss['xSteps']):
		T[0][xIndex] = T[2][xIndex]
		T[ss['ySteps']+1][xIndex] = T[ss['ySteps']-1][xIndex]
		
	
	return T, None
	
	
def updateVt(vt, sf, T, q, ss):

	vtij1 = np.roll(vt, -1, axis=1)
	vtijm1 = np.roll(vt, 1, axis=1)
	vti1j = np.roll(vt, -1, axis=0)
	vtim1j = np.roll(vt, 1, axis=0)
	
	Tij1 = np.roll(T, -1, axis=1)
	Tijm1 = np.roll(T, 1, axis=1)
	Ti1j = np.roll(T, -1, axis=0)
	Tim1j = np.roll(T, 1, axis=0)
	
	
	
	for xIndex in range(ss['xSteps']):
		vt[1][xIndex] =  -2*sf[2][xIndex]/(ss['dy']**2)
		vt[ss['ySteps']][xIndex] =  -2*sf[ss['ySteps']-1][xIndex]/(ss['dy']**2)
	
	pressureTerm = ss['g']*(-ss['alpha'])/(ss['rho0'])* (Tij1-Tijm1)/(2*ss['dx'])
	laplacianTerm = ss['nu'] * (  (vtij1-2*vt+vtijm1)/(ss['dx']**2)  + (vti1j-2*vt+vtim1j)/(ss['dy']**2)  )
	
	vt = ss['dt']*(pressureTerm + laplacianTerm) + vt
	
	for xIndex in range(ss['xSteps']):
		vt[1][xIndex] =  -2*sf[2][xIndex]/(ss['dy']**2)
		vt[ss['ySteps']][xIndex] =  -2*sf[ss['ySteps']-1][xIndex]/(ss['dy']**2)
		
		vt[0][xIndex] = 0
		vt[ss['ySteps']+1][xIndex] = 0
	
	## need to set the boundry conditions still
	return vt
	
	
def printFields(vt, sf, T, ss, time):
	xs = np.arange(0, ss['xSize'] + ss['dx'], ss['dx'])
	ys = np.arange(-ss['dy'], ss['ySize'] + 2*ss['dy'], ss['dy'])
	
	xMesh, yMesh = np.meshgrid(xs, ys)
	
	fig, ax = plt.subplots(dpi=120)
	ax.contourf(xMesh, yMesh, T, 100)
	plt.xlabel("x (horizontal)")
	plt.ylabel("y (vertical)")
	plt.savefig("t/"+str(time)+".png")
	plt.close()
	
	fig, ax = plt.subplots(dpi=120)
	ax.contourf(xMesh, yMesh, sf, 100)
	plt.xlabel("x (horizontal)")
	plt.ylabel("y (vertical)")
	plt.savefig("sfs/"+str(time)+".png")
	plt.close()
	
	fig, ax = plt.subplots(dpi=120)
	ax.contourf(xMesh, yMesh, vt, 100)
	plt.xlabel("x (horizontal)")
	plt.ylabel("y (vertical)")
	plt.savefig("vts/"+str(time)+".png")
	plt.close()
	
	


def main():
	
	
	
	
	ss = dict()
	ss['xSteps'] = 60
	ss['ySteps'] = 60
	ss['xSize'] = 1
	ss['ySize'] = 1
	
	ss['dx'] = ss['xSize']/(ss['xSteps'] - 1)
	ss['dy'] = ss['ySize']/(ss['ySteps'] - 1)
	
	ss['poissonIterations'] = 10000
	ss['poissonError'] = 0.0001
	ss['SORparam'] = 0.5
	ss['dt'] = 0.001
	ss['timeSteps'] = 100000
	ss['kappa'] = 0
	ss['rho0']= 1000
	ss['cv'] = 4000
	ss['g'] = -100
	ss['alpha'] = 0.1
	ss['nu']=0.0006
	
	
	
	
	T = np.zeros((ss['ySteps']+2, ss['xSteps']))
	sf = np.zeros((ss['ySteps']+2, ss['xSteps']))
	q = np.zeros((ss['ySteps']+2, ss['xSteps']))
	vt = np.zeros((ss['ySteps']+2, ss['xSteps']))
	
	
	
	for xIndex in range(int(0.4*ss['xSteps']), int(0.6*ss['xSteps'])):
		for yIndex in range(int(0.4*ss['ySteps']), int(0.6*ss['ySteps'])):
			q[yIndex][xIndex] = 0
	
	
	
	
	step = 0
	while step < ss['timeSteps']:
	
		oldsf = np.copy(sf)
		oldT = np.copy(T)
		oldvt = np.copy(vt)
		
		
		print(step)
		if (step % 100 == 0):
			printFields(vt, sf, T, ss, step)
		
		sf = updateSf(vt, sf, T, q, ss)
		
		T, error = updateT(vt, sf, T, q, ss)
		
		if (error == "time error"):
			sf = oldsf
			T = oldT
			vt = oldvt
			ss['dt'] = ss['dt']/2
			print("new dt", ss['dt'])
			continue;
			
			
		
		for xIndex in range(int(0.4*ss['xSteps']), int(0.6*ss['xSteps'])):
			for yIndex in range(int(0.4*ss['ySteps']), int(0.6*ss['ySteps'])):
				T[yIndex][xIndex] = 100
		
		
		vt = updateVt(vt, sf, T, q, ss)
		
		step += 1
		
		
		
		
		
		
main()
		
		
	
	
	
	

	





