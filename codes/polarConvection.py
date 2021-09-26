"""
Polar Convection code

Written By Maximilian Williams
"""

import matplotlib.pyplot as plt
import time
import copy
import math


def updateStreamFunction(vorticity, streamfunction, simulationSettings):
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
	dPhi = simulationSettings['dPhi']
	dRadius = simulationSettings['dRadius']
	iterations = simulationSettings['poissonIterations']
	epsilon = simulationSettings['poissonError']
	omega = simulationSettings['SORParam']
	
	for interation in range(iterations):
		oldstreamfunction = copy.copy(streamfunction)
		for rIndex in range(1,simulationSettings['radiusSteps']-1):
			for phiIndex in range(simulationSettings['phiSteps']):
				r = rIndexToR(simulationSettings, rIndex)
				term1 = 1/r * (oldstreamfunction[rIndex+1][phiIndex] - oldstreamfunction[rIndex-1][phiIndex])/(2*dRadius)
				term2 = (oldstreamfunction[rIndex+1][phiIndex] + oldstreamfunction[rIndex-1][phiIndex])/(dRadius**2)
				term3 = vorticity[rIndex][phiIndex]
				term4 = 1/r**2 *( oldstreamfunction[rIndex][(phiIndex+1)%simulationSettings['phiSteps']]+oldstreamfunction[rIndex][(phiIndex-1)%simulationSettings['phiSteps']] )/(dPhi**2)
				streamfunction[rIndex][phiIndex] = (1-omega)*oldstreamfunction[rIndex][phiIndex] + omega * (dRadius * r * dPhi)**2/(2*(r**2 * dPhi**2 + dRadius**2)) *(term1 + term2 + term3 + term4)
				
		epsilonDash = 0
		for rIndex in range(simulationSettings['radiusSteps']):
			for phiIndex in range(simulationSettings['phiSteps']):
				epsilonDash += abs(oldstreamfunction[rIndex][phiIndex] - streamfunction[rIndex][phiIndex])
				
		if (epsilonDash < epsilon):
			break;
		if (interation == iterations - 1):
			raise Exception("When solving for streamfunction there was no convergence")
			
	return streamfunction
	

	
def updateTemperature(streamfunction, temperature, heat, simulationSettings):
	Cp = simulationSettings['Cp']
	kappa = simulationSettings['kappa']
	dt = simulationSettings['dt']
	dPhi = simulationSettings['dPhi']
	dRadius = simulationSettings['dRadius']
	iterations = simulationSettings['poissonIterations']
	epsilon = simulationSettings['poissonError']
	omega = simulationSettings['SORParam']
	
	for phiIndex in range(simulationSettings['phiSteps']):
		temperature[0][phiIndex] = temperature[1][phiIndex]
		temperature[simulationSettings['radiusSteps']-1][phiIndex] = temperature[simulationSettings['radiusSteps']-2][phiIndex]

	oldTemperature = copy.copy(temperature)
	for rIndex in range(1,simulationSettings['radiusSteps']-1):
		for phiIndex in range(simulationSettings['phiSteps']):
			r = rIndexToR(simulationSettings, rIndex)
			term1 = (oldTemperature[rIndex+1][phiIndex] - 2*oldTemperature[rIndex][phiIndex] + oldTemperature[rIndex-1][phiIndex])/(dRadius**2)
			term2 = 1/r *(oldTemperature[rIndex+1][phiIndex] - oldTemperature[rIndex-1][phiIndex])/(2*dRadius)
			term3 = 1/r**2 *(oldTemperature[rIndex][(phiIndex+1)%simulationSettings['phiSteps']] - 2*oldTemperature[rIndex][phiIndex] + oldTemperature[rIndex][ (phiIndex-1)%simulationSettings['phiSteps']]  )/(dPhi**2)
			term4 = heat[rIndex][phiIndex]/Cp
			term5 = -1/r * (streamfunction[rIndex][(phiIndex+1)%simulationSettings['phiSteps']] - streamfunction[rIndex][(phiIndex-1)%simulationSettings['phiSteps']])/(2*dPhi) *(oldTemperature[rIndex+1][phiIndex] - oldTemperature[rIndex-1][phiIndex])/(2*dRadius)
			term6 = 1/r * (streamfunction[rIndex+1][phiIndex] - streamfunction[rIndex-1][phiIndex])/(2*dRadius) *(oldTemperature[rIndex][(phiIndex+1)%simulationSettings['phiSteps']] - oldTemperature[rIndex][(phiIndex-1)%simulationSettings['phiSteps']])/(2*dPhi)
			
			update = kappa*(term1 + term2 + term3 ) + term4 + term5 + term6
			
			temperature[rIndex][phiIndex]  = oldTemperature[rIndex][phiIndex] + dt*update
	return temperature

def updateVorticity(vorticity, temperature, streamfunction, simulationSettings):
	## here we update the vorticity
	
	nu = simulationSettings['nu']
	dt = simulationSettings['dt']
	G = simulationSettings['G']
	rho0 = simulationSettings['rho0']
	dPhi = simulationSettings['dPhi']
	dRadius = simulationSettings['dRadius']
	alpha = simulationSettings['alpha']
	
	## we need to update the boundries
	
	for index in range(simulationSettings['phiSteps']):
		vorticity[0][index] = (streamfunction[0][index] - streamfunction[1][index])*2/(dRadius**2)
		vorticity[simulationSettings['radiusSteps']-1][index] = ( streamfunction[simulationSettings['radiusSteps']-1][index] -  streamfunction[simulationSettings['radiusSteps']-2][index])* 2/(dRadius**2)
		
	oldVorticity = copy.copy(vorticity)
	
	for rIndex in range(1, simulationSettings['radiusSteps']-1):
		for phiIndex in range(simulationSettings['phiSteps']):
			
			r = rIndexToR(simulationSettings, rIndex)
			
			g = 4/3 *math.pi * G* rho0*r
			term1 = -g/(rho0 * r) * (-rho0*alpha)*(temperature[rIndex][(phiIndex+1)%simulationSettings['phiSteps']] - temperature[rIndex][(phiIndex-1)%simulationSettings['phiSteps']])/(2*dPhi)
			term2 = (oldVorticity[rIndex+1][phiIndex] - 2*oldVorticity[rIndex][phiIndex] + oldVorticity[rIndex-1][phiIndex])/(dRadius**2)
			term3 = 1/r * (oldVorticity[rIndex+1][phiIndex] - oldVorticity[rIndex-1][phiIndex])/(2*dRadius)
			term4 = 1/r**2 * (oldVorticity[rIndex][(phiIndex+1)%simulationSettings['phiSteps']] - 2*oldVorticity[rIndex][phiIndex] + oldVorticity[rIndex][(phiIndex-1)%simulationSettings['phiSteps']])/(dPhi**2)
			
			update = term1 + nu * (term2 + term3 + term4)
			
			vorticity[rIndex][phiIndex] = oldVorticity[rIndex][phiIndex] + dt*update
	
	return vorticity
	
def rIndexToR(simulationSettings, rIndex):
	return simulationSettings['innerRadius'] + rIndex * simulationSettings['dRadius']
	
def phiIndexToPhi(simulationSettings, phiIndex):
	return phiIndex * simulationSettings['dPhi']
	
def generateFields(simulationSettings):
	streamfunction = [[0 for phiIndex in range(simulationSettings['phiSteps'])] for rIndex in range(simulationSettings['radiusSteps'])]
	
	vorticity = [[0 for phiIndex in range(simulationSettings['phiSteps'])] for rIndex in range(simulationSettings['radiusSteps'])]
	
	temperature = [[0 for phiIndex in range(simulationSettings['phiSteps'])] for rIndex in range(simulationSettings['radiusSteps'])]
	
	heat = [[0 for phiIndex in range(simulationSettings['phiSteps'])] for rIndex in range(simulationSettings['radiusSteps'])]
	
	return streamfunction, vorticity, temperature, heat
	
	
def main():
	
	## settings for the simulation, how large our space is, how fine our grid is
	simulationSettings = dict()
	simulationSettings['outerRadius'] = 1
	simulationSettings['innerRadius'] = 0.1
	simulationSettings['radiusSteps'] = 100
	simulationSettings['phiSteps'] = 100
	simulationSettings['dRadius'] = simulationSettings['radiusSteps']/(simulationSettings['outerRadius'] - simulationSettings['innerRadius'])
	simulationSettings['dPhi'] = 2*math.pi/simulationSettings['phiSteps']
	simulationSettings['poissonIterations'] = 100
	simulationSettings['poissonError'] = 0.001
	simulationSettings['SORParam'] = 0.5
	simulationSettings['dt'] = 0.0001
	simulationSettings['Cp'] = 1
	simulationSettings['kappa'] = 0.00001
	simulationSettings['nu'] = 0.01
	simulationSettings['G'] = 1
	simulationSettings['rho0'] = 1000
	simulationSettings['alpha'] = 10
	streamfunction, vorticity, temperature, heat = generateFields(simulationSettings)
	
	for index in range(simulationSettings['phiSteps']):
		heat[int(simulationSettings['phiSteps']/2)][index] = 100
	
	for index in range(10000):
		streamfunction = updateStreamFunction(vorticity, streamfunction, simulationSettings)
		
		temperature = updateTemperature(streamfunction, temperature, heat, simulationSettings)
		
		vorticity = updateVorticity(vorticity, temperature, streamfunction, simulationSettings)
		
		plt.imshow(temperature)
		if (index % 100 == 0):
			plt.savefig('tpolar/'+str(index)+".png")
		
		
	
	
	
	
	

		
		
		
		
main()
		
		

	
	
	
			


	
			
	
	
	
	
	
	
	

	
	
			
			

				
				
	
		
	
	
	

