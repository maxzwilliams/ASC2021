"""
Fluids confined to a box with a thermal source at the bottom
"""

import matplotlib.pyplot as plt
import time
import copy
import math


def updateStreamFunction(vorticity, streamfunction, n, iterations, epsilon, omega,h):
	"""
	Function that finds the streamfunction using the SOR method in cartesian coordiantes
	
	Parameters:
		vorticity (n lists of size n): vorticity at each point in the domein
		streamfunction (n lists of size n): streamfunction at each point in the domain
		n (int): size of the domain
		iterations (int): number of iterations used for the SOR solver
		epsilon (float): acceptable error in streamfunction solution
	 	omega (float64): SOR parameter
	 
	 Returns:
	 	streamfunction: new streamfunction
	"""
	
	for iteration in range(iterations):
		oldstreamfunction = copy.deepcopy(streamfunction)
		for i in range(1,n-1):
			for j in range(1,n-1):
				streamfunction[i][j] = 0.25*omega*(oldstreamfunction[i+1][j]+oldstreamfunction[i-1][j]+oldstreamfunction[i][j+1]+oldstreamfunction[i][j-1]+h*h*vorticity[i][j])+ (1.0-omega)*oldstreamfunction[i][j]
				
		epsilonDash = 0
		for i in range(n):
			for j in range(n):
				epsilonDash += abs(oldstreamfunction[i][j] - streamfunction[i][j])
				
		if (epsilonDash < epsilon):
			break;
		if (iteration == iterations-1):
			raise Exception("when solving for streamfunction the method did not converge")
		
	return streamfunction
	
	
def updateTemperature(streamfunction, temperature, heat, h, n, Cp, kappa, Touter, dt):
	"""
	Function that updates the temperature field
	"""
	
	update = [[0 for i in range(n)] for j in range(n)]
	newTemperature = [[0 for i in range(n)] for j in range(n)]
	
	for i in range(n):
		temperature[0][i] = Touter
		newTemperature[0][i] = Touter
		temperature[n-1][i] = Touter
		newTemperature[n-1][i] = Touter
		temperature[i][0] = Touter
		newTemperature[n-1][i] = Touter
	
	
	
	for i in range(1, n-1):
		for j in range(1, n-1):
			term1 = (streamfunction[i][j+1] - streamfunction[i][j-1])/(2*h) * (temperature[i+1][j] - temperature[i-1][j])/(2*h)
			term2 = - (streamfunction[i+1][j] - streamfunction[i-1][j])/(2*h) * (temperature[i][j+1] - temperature[i][j-1])/(2*h)
			term3 = kappa * (  (temperature[i+1][j] - 2*temperature[i][j] + temperature[i-1][j])/(h**2) + (temperature[i][j+1] - 2*temperature[i][j] + temperature[i][j-1])/(h**2)  )
			term4 = heat[i][j]/Cp
			
			update[i][j] = term1 + term2 + term3 + term4
			
			newTemperature[i][j] = temperature[i][j] + dt*update[i][j]
			
			
	## what sort of boundry conditions should I do??
	

	
	## now I need to do the surface:
	## I will hold the surface at a constant temperature		
	return newTemperature
	
## now we need to write something that updates the vorticity equation for creeping flow

def updateVorticity(vorticity, temperature, streamfunction, nu, alpha, g, rho0,h,n,dt):
	newVorticity = [[0 for i in range(n)] for j in range(n)]
	update = [[0 for i in range(n)] for j in range(n)]	
	## Now we need to look at eh bounrdy conditions.
	for i in range(1, n-1):
		vorticity[i][0] = -2*streamfunction[i][1]/(h**2)
		
		vorticity[i][n-1] = -2*streamfunction[i][n-2]/(h**2)
		
		vorticity[0][i] = -2*streamfunction[1][i]/(h**2)
		
		vorticity[n-1][i]=  -2*streamfunction[n-2][i]/(h**2) 
		
		newVorticity[i][0] = -2*streamfunction[i][1]/(h**2)
		
		newVorticity[i][n-1] = -2*streamfunction[i][n-2]/(h**2)
		
		newVorticity[0][i] = -2*streamfunction[1][i]/(h**2)
		
		newVorticity[n-1][i]=  -2*streamfunction[n-2][i]/(h**2) 
		
		

	for i in range(1, n-1):
		for j in range(1, n-1):
			update[i][j] = g*alpha/rho0 * ( temperature[i+1][j] - temperature[i-1][j])/(2*h) + nu*( (vorticity[i+1][j] - 2*vorticity[i][j] + vorticity[i-1][j]  )/(h**2) + (vorticity[i][j+1] -2*vorticity[i][j] + vorticity[i][j-1])/(h**2) ) 
			
	
	for i in range(1, n-1):
		for j in range(1, n-1):
			newVorticity[i][j] = vorticity[i][j] + dt*update[i][j]
			

		
	return  newVorticity
	
	
	
def main():

	## function that ties everything together
	n = 30
	timesteps = 10000
	dt = 0.001
	
	iterations = 1000
	epsilon = 0.01
	omega=0.5
	
	h = 1/(n-1)
	Cp=1
	kappa=1
	Touter=0
	rho0=1000
	alpha=0.01
	nu=0.0005
	g=10
	t = 0
	
	
		
	
	temperature = [[0 for i in range(n)] for i in range(n)]
	streamfunction = [[0 for i in range(n)] for i in range(n)]
	vorticity = [[0 for i in range(n)] for i in range(n)]
	
	heat = [[0 for i in range(n)] for i in range(n)]
	
	for i in range(n):
		heat[i][3] =  1
	
	for timestep in range(timesteps):
		print("loop: " + str(timestep))
		
		streamfunction = updateStreamFunction(vorticity, streamfunction, n, iterations, epsilon, omega,h)
		
		temperature = updateTemperature(streamfunction, temperature, heat, h, n, Cp, kappa, Touter,dt)
		
		vorticity = updateVorticity(vorticity, temperature, streamfunction, nu, alpha, g, rho0,h,n,dt)
		
		t = t + dt
		plt.contour(vorticity)
		plt.savefig("vts//file"+str(timestep)+".png")
		plt.clf()
		
		plt.contour(streamfunction)
		plt.savefig("sfs//file"+str(timestep)+".png")
		plt.clf()
		
		
main()
		
		

	
	
	
			


	
			
	
	
	
	
	
	
	

	
	
			
			

				
				
	
		
	
	
	

