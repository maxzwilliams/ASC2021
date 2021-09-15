"""
Fluids confined to a box with a thermal source at the bottom
"""

import matplotlib.pyplot as plt
import time
import copy
import math


def updateStreamFunction(vorticity, streamfunction, n, iterations, epsilon, omega):
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
	
	
def updateTemperature(streamfunction, temperature, heat, h, n, Cp, kappa, Ttop):
	"""
	Function that updates the temperature field
	"""
	
	update = [[0 for i in range(n)] for j in range(n)]
	newTemperature = [[0 for i in range(n)] for j in range(n)]
	
	for i in range(1, n-1):
		for j in range(1, n-1):
			term1 = (streamfunction[i][j+1] - streamfunction[i][j-1])/(2 h) * (temperature[i+1][j] - temperature[i-1][j])/(2*h)
			term2 = - (streamfunction[i+1][j] - streamfunction[i-1][j])/(2*h) * (temperature[i][j+1] - temperature[i][j-1])/(2*h)
			term3 = kappa * (  (temperature[i+1][j] - 2*temperature[i][j] + temperature[i-1][j])/(h**2) + (temperature[i][j+1] - 2*temperature[i][j] + temperature[i][j-1])/(h**2)  )
			term4 = heat[i][j]/Cp
			
			update[i][j] = term1 + term2 + term3 + term4
			
			newTemperature[i][j] = temperature[i][j] + dt*update[i][j]
			
			
	## what sort of boundry conditions should I do??
	
	## I will hold the surface at a constant temperature
	for i in range(n):
		newTemperature[i][n-1] = Ttop
	
	## I will have the left, right and bottom sides as insulating
	for i in range(1, n-1):
		newTemperature[0][i] = newTemperature[1][i]
		newTemperature[n-1][i] = newTemperature[n-2][i]
		newTemperature[i][0] = newTemperature[i][1]
		
	newTemperature[0][0] = 0.5*( newTemperature[0][1] + newTemperature[1][0] )
	newTemperature[n-1][0] = 0.5*( newTemperature[n-2][0] + newTemperature[n-1][1] )
	
	
	newTemperature[0][n-1] = 0.5 * (newTemperature[0][n-2] + newTemperature[1][n-1])
	newTemperature[n-1][n-1] = 0.5 * (newTemperature[n-1][n-2] + newTemperature[n-2][n-1])
	
	## now I need to do the surface:
	## I will hold the surface at a constant temperature		
	return newTemperature[i][j]
	
## now we need to write something that updates the vorticity equation for creeping flow


	
			
	
	
	
	
	
	
	

	
	
			
			

				
				
	
		
	
	
	

