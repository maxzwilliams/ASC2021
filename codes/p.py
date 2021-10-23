"""
Polar Convection code

Written By Maximilian Williams
"""

import matplotlib.pyplot as plt
import time
import copy
import math
import random
import numpy as np


def updateStreamFunctionFast(vorticity, streamfunction, simulationSettings):
	"""
	updates the streamfunction using the jacobi method with SOR.	
	"""

	vorticity = np.array(vorticity)
	streamfunction = np.array(streamfunction)
	dPhi = simulationSettings['dPhi']
	dRadius = simulationSettings['dRadius']
	iterations = simulationSettings['poissonIterations']
	epsilon = simulationSettings['poissonError']
	omega = simulationSettings['SORParam']
	
	
	
	for iteration in range(iterations):
		ogsf = np.copy(streamfunction)
		## below works
		streamfunctionij1 = np.array(streamfunction)
		streamfunctionij1 = np.roll(streamfunctionij1, -1,axis=1) ## index 0 is the first index of the original

		## below works		
		streamfunctionijm1 = np.array(streamfunction)
		streamfunctionijm1 = np.roll(streamfunctionijm1, 1,axis=1) ## index 0 is the last index of the original
		
		## works
		streamfunctioni1j = np.array(streamfunction)
		streamfunctioni1j =  np.roll(streamfunctioni1j, -1, axis=0) ## index 0 is the first index of the original
		
		## works
		streamfunctionim1j = np.array(streamfunction)
		streamfunctionim1j = np.roll(streamfunctionim1j, 1, axis=0) ## index 0 is the last index of the original
		
		
		
		
		
		## index 0 on this is the last entry on streamfunction
		rs = [ 1/rIndexToR(simulationSettings, rIndex) for rIndex in range(simulationSettings['radiusSteps']) ]
		rss = [1/rIndexToR(simulationSettings, rIndex)**2 for rIndex in range(simulationSettings['radiusSteps']) ]
		updateCoeffecient = [rIndexToR(simulationSettings, rIndex)**2/(rIndexToR(simulationSettings, rIndex)**2*dPhi**2 + dRadius**2) for rIndex in range(simulationSettings['radiusSteps']) ]
		
		## good
		Irs = np.array([rs for index in range(simulationSettings['phiSteps'])])
		Irs = np.transpose(Irs)

		## good
		Irss = np.array([rss for index in range(simulationSettings['phiSteps'])])
		Irss = np.transpose(Irss)
		
		## not checked but looks correct
		updateCoeffecientMatrix = np.array([updateCoeffecient for index in range(simulationSettings['phiSteps'])])
		updateCoeffecientMatrix = np.transpose(updateCoeffecientMatrix)
		
		
		term1 = np.multiply(Irs, (streamfunctioni1j - streamfunctionim1j)   *(1/(2*dRadius)))
		
		term2 = (streamfunctioni1j + streamfunctionim1j) * 1/(dRadius**2)
		term3 = vorticity
		term4 = np.multiply(Irss, streamfunctionij1 + streamfunctionijm1) * (1/(dPhi**2))
		streamfunction = (1-omega)*streamfunction + omega * ((dRadius*dPhi)**2)/2 * np.multiply(updateCoeffecientMatrix, term1+term2+term3+term4)
		
		for rIndex in range(simulationSettings['radiusSteps']):
			streamfunction[rIndex][0] = 0
			streamfunction[rIndex][simulationSettings['phiSteps']-1] = 0
		for phiIndex in range(simulationSettings['phiSteps']):
			streamfunction[0][phiIndex] = 0
			streamfunction[simulationSettings['radiusSteps']-1][phiIndex] = 0
			
		difference = np.subtract(ogsf,streamfunction)
		epsilonDash = 0
		for rIndex in range(simulationSettings['radiusSteps']):
			for phiIndex in range(simulationSettings['phiSteps']):
				epsilonDash += abs(difference[rIndex][phiIndex])
		##print("iteration", iteration, "error", epsilonDash)
		if (epsilonDash < epsilon):
			break;
			
			
		##print(streamfunction)
		##time.sleep(1)
	##print(streamfunction)
	##time.sleep(1)
	return streamfunction
		
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
		oldstreamfunction = copy.deepcopy(streamfunction)
		for rIndex in range(1,simulationSettings['radiusSteps']-1):
			for phiIndex in range(simulationSettings['phiSteps']):
				r = rIndexToR(simulationSettings, rIndex)
				term1 = 1/r * (streamfunction[rIndex+1][phiIndex] - streamfunction[rIndex-1][phiIndex])/(2*dRadius)
				term2 = (streamfunction[rIndex+1][phiIndex] + streamfunction[rIndex-1][phiIndex])/(dRadius**2)
				term3 = vorticity[rIndex][phiIndex]
				term4 = 1/r**2 *( streamfunction[rIndex][(phiIndex+1)%simulationSettings['phiSteps']]+streamfunction[rIndex][(phiIndex-1)%simulationSettings['phiSteps']] )/(dPhi**2)
				streamfunction[rIndex][phiIndex] = (1-omega)*streamfunction[rIndex][phiIndex] + omega * (dRadius * r * dPhi)**2/(2*(r**2 * dPhi**2 + dRadius**2)) *(term1 + term2 + term3 + term4)
				
				if (math.isnan(streamfunction[rIndex][phiIndex])):
					print("nan detected")
					print(OGsf[rIndex][phiIndex])
					print(oldstreamfunction[rIndex][phiIndex])
					print(term1)
					print(term2)
					print(term3)
					print(term4)
					
					time.sleep(10000)
				
		epsilonDash = 0
		for rIndex in range(simulationSettings['radiusSteps']):
			for phiIndex in range(simulationSettings['phiSteps']):
				epsilonDash += abs(oldstreamfunction[rIndex][phiIndex] - streamfunction[rIndex][phiIndex])
				
		##print("iteration:", interation, "error:", epsilonDash)
				
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
	rho0 = simulationSettings['rho0']
	
	##for phiIndex in range(simulationSettings['phiSteps']):
		##temperature[0][phiIndex] = temperature[1][phiIndex]
		##temperature[simulationSettings['radiusSteps']-1][phiIndex] = temperature[simulationSettings['radiusSteps']-2][phiIndex]

	oldTemperature = copy.deepcopy(temperature)
	for rIndex in range(0, simulationSettings['radiusSteps']):
		for phiIndex in range(simulationSettings['phiSteps']):
			r = rIndexToR(simulationSettings, rIndex)
			if (rIndex == 0):
				preTemp = oldTemperature[1][phiIndex]
				advTemp = oldTemperature[rIndex+1][phiIndex]
				u = 0
				v = 0
			
			
			elif (rIndex == simulationSettings['radiusSteps']-1):
				advTemp = oldTemperature[simulationSettings['radiusSteps']-2][phiIndex]
				preTemp = oldTemperature[rIndex-1][phiIndex]
				u = 0 
				v = 0
			else:
				advTemp = oldTemperature[rIndex+1][phiIndex]
				preTemp = oldTemperature[rIndex-1][phiIndex]
				u = 1/r * (streamfunction[rIndex][(phiIndex+1)%simulationSettings['phiSteps']] - streamfunction[rIndex][(phiIndex-1)%simulationSettings['phiSteps']])/(2*dPhi)
				v = - (streamfunction[rIndex+1][phiIndex] - streamfunction[rIndex-1][phiIndex])/(2*dRadius)
			
			
			term1 = ( advTemp - 2*oldTemperature[rIndex][phiIndex] + preTemp)/(dRadius**2)
			term2 = 1/r *(advTemp - preTemp)/(2*dRadius)
			term3 = 1/r**2 *(oldTemperature[rIndex][(phiIndex+1)%simulationSettings['phiSteps']] - 2*oldTemperature[rIndex][phiIndex] + oldTemperature[rIndex][ (phiIndex-1)%simulationSettings['phiSteps']]  )/(dPhi**2)
			term4 = heat[rIndex][phiIndex]
			
			
			if (abs(u)*dt/dRadius >= 1 or abs(v)*dt/(r*dPhi) >= 1):
				return temperature, "Time Fail"
			
			
			if (u>=0):
				term5 = -u * (oldTemperature[rIndex][phiIndex] - preTemp)/(dRadius)
			else:
				term5 = -u * (advTemp - oldTemperature[rIndex][phiIndex])/(dRadius)
			if (v >=0):
				term6 = -1/r * v * (oldTemperature[rIndex][(phiIndex)%simulationSettings['phiSteps']] - oldTemperature[rIndex][(phiIndex-1)%simulationSettings['phiSteps']])/(dPhi)
			else:
				term6 = -1/r * v * (oldTemperature[rIndex][(phiIndex+1)%simulationSettings['phiSteps']] - oldTemperature[rIndex][(phiIndex)%simulationSettings['phiSteps']])/(dPhi)
			##term5 = -1/r * (streamfunction[rIndex][(phiIndex+1)%simulationSettings['phiSteps']] - streamfunction[rIndex][(phiIndex-1)%simulationSettings['phiSteps']])/(2*dPhi) *(oldTemperature[rIndex+1][phiIndex] - oldTemperature[rIndex-1][phiIndex])/(2*dRadius)
			##term6 = 1/r * (streamfunction[rIndex+1][phiIndex] - streamfunction[rIndex-1][phiIndex])/(2*dRadius) *(oldTemperature[rIndex][(phiIndex+1)%simulationSettings['phiSteps']] - oldTemperature[rIndex][(phiIndex-1)%simulationSettings['phiSteps']])/(2*dPhi)
			update = kappa*(term1 + term2 + term3 ) + term4 + term5 + term6
			##print("update", update)
			temperature[rIndex][phiIndex] = oldTemperature[rIndex][phiIndex] + dt*update
			
	return temperature, None

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
		
	oldVorticity = copy.deepcopy(vorticity)
	
	for rIndex in range(1, simulationSettings['radiusSteps']-1):
		for phiIndex in range(simulationSettings['phiSteps']):
			
			r = rIndexToR(simulationSettings, rIndex)
			
			g = 4/3 *math.pi * G*  rho0*r
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
	
def myRound(number):
	if (number - math.floor(number)>=0.5 ):
		return int(math.ceil(number))
	else:
		return int(math.floor(number))
		
	
def plotFields(temperature, vorticity, streamfunction, rMesh, thetaMesh, t, index, analyticDynamics):
	rCenter = 0
	thetaCenter = 0
	
	totalTemperature = 0
	for rIndex in range(len(temperature)):
		for thetaIndex in range(len(temperature[0])):
			totalTemperature += temperature[rIndex][thetaIndex]
	
	xPosition = 0
	yPosition = 0
	for rIndex in range(len(temperature)):
		for thetaIndex in range(len(temperature[0])):

			r = rMesh[rIndex]
			theta = thetaMesh[thetaIndex]

			x = r*math.cos(theta)
			y = r*math.sin(theta)

			xPosition += x*temperature[rIndex][thetaIndex]/totalTemperature
			yPosition += y*temperature[rIndex][thetaIndex]/totalTemperature




			rCenter += rMesh[rIndex]*temperature[rIndex][thetaIndex]/totalTemperature
			##thetaCenter += (thetaMesh[thetaIndex] )*temperature[rIndex][thetaIndex]/totalTemperature

	##rCenter = (xPosition**2 + yPosition**2)**(0.5)
	if (xPosition > 0 and yPosition >= 0):
		thetaCenter = math.atan(abs(yPosition/xPosition))
	if (xPosition < 0 and yPosition >= 0):
		thetaCenter = math.pi - math.atan(abs(yPosition/xPosition))
	if (xPosition < 0 and yPosition < 0):
		thetaCenter = math.pi + math.atan(abs(yPosition/xPosition))
	if (xPosition > 0 and yPosition < 0):
		thetaCenter = 2*math.pi - math.atan(abs(yPosition/xPosition))

			
	##rCenter = myRound(rCenter)
	##thetaCenter = myRound(thetaCenter)
	fig, ax = plt.subplots(dpi=120, subplot_kw=dict(projection='polar'))
	pc = ax.pcolormesh(thetaMesh, rMesh, temperature, antialiased=False)
	theta = [thetaCenter%(2*math.pi)]
	r = [rCenter]
	##plt.plot(theta, r, 'or')
	ax.set_rorigin(-0.1)
	cb = plt.colorbar(pc)
	cb.set_label("Temperature (arb. units)")
	plt.title(t)
	plt.savefig('tpolar/'+str(index)+".png")
	plt.close()
	plt.clf()
	fig.clear()
	
	fig, ax = plt.subplots(dpi=120)
	pc = ax.pcolormesh(thetaMesh, rMesh, temperature, antialiased=False)
	theta = [thetaCenter%(2*math.pi)]
	r = [rCenter]
	##plt.plot(theta, r, 'or')
	##ax.set_rorigin(-0.1)
	cb = plt.colorbar(pc)
	cb.set_label("Temperature (arb. units)")
	plt.title(t)
	plt.xlabel("Azimuthal Angle")
	plt.ylabel("Radius (meters)")
	plt.savefig('tpolar/'+str(index)+"CProjection.png")
	plt.close()
	plt.clf()
	fig.clear()
	
	
	fig, ax = plt.subplots(dpi=120, subplot_kw=dict(projection='polar'))
	pc = ax.pcolormesh(thetaMesh, rMesh, streamfunction, antialiased=False)
	ax.set_rorigin(-0.1)
	plt.title(t)
	cb = plt.colorbar(pc)
	cb.set_label("Streamfunction (m^2/s)")
	plt.savefig('sfspolar/'+str(index)+".png")
	plt.close()
	plt.clf()
	fig.clear()
	
	fig, ax = plt.subplots(dpi=120, subplot_kw=dict(projection='polar'))
	pc = ax.pcolormesh(thetaMesh, rMesh, vorticity, antialiased=False)
	ax.set_rorigin(-0.1)
	plt.title(t)
	cb = plt.colorbar(pc)
	cb.set_label("Vorticity (1/s)")
	plt.savefig('vtspolar/'+str(index)+".png")
	plt.close()
	plt.clf()
	fig.clear()

	

def getTotalQ(T, simulationSettings):
	s = 0
	for rIndex in range(simulationSettings['radiusSteps']):
		for phiIndex in range(simulationSettings['phiSteps']):
			r = simulationSettings['innerRadius'] + rIndex*simulationSettings['dRadius']
			area = r * simulationSettings['dRadius']*simulationSettings['dPhi']
			s += area * T[rIndex][phiIndex]
			
	return s
	
def main():



	vel = 0.0005
	## settings for the simulation, how large our space is, how fine our grid is
	simulationSettings = dict()
	simulationSettings['outerRadius'] = 1
	simulationSettings['innerRadius'] = 0.1
	simulationSettings['radiusSteps'] = 100
	simulationSettings['phiSteps'] = 100
	simulationSettings['dRadius'] = (simulationSettings['outerRadius'] - simulationSettings['innerRadius'])/(simulationSettings['radiusSteps']-1)
	##simulationSettings['dRadius'] = simulationSettings['radiusSteps']/(simulationSettings['outerRadius'] - simulationSettings['innerRadius'])
	simulationSettings['dPhi'] = 2*math.pi/(simulationSettings['phiSteps']-1)
	simulationSettings['poissonIterations'] = 1000
	simulationSettings['poissonError'] = 0.001
	simulationSettings['SORParam'] = 0.5
	simulationSettings['dt'] = 10
	simulationSettings['Cp'] = 4000
	simulationSettings['kappa'] = 0.000001
	simulationSettings['rho0'] = 1000
	simulationSettings['nu'] = 1 * simulationSettings['kappa'] 
	simulationSettings['G'] = 3/(4 * math.pi * simulationSettings['rho0'])
	simulationSettings['Heat'] = 0.0001
	
	simulationSettings['alpha'] = 0.000001
	rs = np.arange(simulationSettings['innerRadius'], simulationSettings['outerRadius'], simulationSettings['dRadius'] )
	thetas = np.arange(0, 2*math.pi+simulationSettings['dPhi'], simulationSettings['dPhi'] )



	
	rmesh, thetaMesh = np.meshgrid(rs, thetas)
	
	totalQ = []
	time = []
	



	maxG = 4/3 *math.pi * simulationSettings['G'] * simulationSettings['rho0'] * simulationSettings['outerRadius']




	streamfunction, vorticity, temperature, heat = generateFields(simulationSettings)
	
	analyticDynamics = dict()


	Pr = simulationSettings['nu']/simulationSettings['kappa']
	Ra = abs(maxG) * simulationSettings['alpha'] * simulationSettings['Heat'] * simulationSettings['outerRadius']**5 /(simulationSettings['nu'] * simulationSettings['kappa']**2)

	print("stats about this simulation")
	print("Pr:", Pr, "Ra", Ra)

	t = 0
	
	##for phiIndex in range(int(0.4*simulationSettings['radiusSteps']),int(0.6*simulationSettings['radiusSteps']) ):
		##for rIndex in range(int(0.4*simulationSettings['phiSteps']),int(0.6*simulationSettings['phiSteps']) ):
			##temperature[rIndex][phiIndex] = 1
	

	rand1 = [[0 for _ in range(simulationSettings['phiSteps'])] for _ in range(simulationSettings['radiusSteps'])]
	rand2 = [[0 for _ in range(simulationSettings['phiSteps'])] for _ in range(simulationSettings['radiusSteps'])]
	for rIndex in range(simulationSettings['radiusSteps']):
		for phiIndex in range(simulationSettings['phiSteps']):
			rand1[rIndex][phiIndex] = (2.0*random.uniform(0,1) -1.0)
			rand2[rIndex][phiIndex] = (2.0*random.uniform(0,1) -1.0)


	index = 0
	while index<100000:	
		print("index", index, end="\r")
		totalQ.append(getTotalQ(temperature, simulationSettings))
		time.append(t)
		
			
		if (index % 100 == 0):
			plotFields(temperature, vorticity, streamfunction, rs, thetas, t, index, analyticDynamics)
			plt.plot(time, totalQ)
			plt.xlabel("time (arb. units)")
			plt.ylabel("Total Thermal Energy (arb. units)")
			plt.savefig("polarTotalQ.png")
			
		streamfunction = updateStreamFunction(vorticity, streamfunction, simulationSettings)
		
		oldstreamfunction = copy.deepcopy(streamfunction)
		oldtemperature = copy.deepcopy(temperature)
		oldvorticity = copy.deepcopy(vorticity)
		
		
		
		##streamfunction = updateStreamFunctionFast(vorticity, streamfunction, simulationSettings)
		
		temperature, state = updateTemperature(streamfunction, temperature, heat, simulationSettings)

		for phiIndex in range(simulationSettings['phiSteps']):
			for rIndex in range(simulationSettings['radiusSteps']):
				da = (simulationSettings['innerRadius'] + rIndex*simulationSettings['dRadius']) * simulationSettings['dRadius'] * simulationSettings['dPhi']
				if (rIndex < 45):
					temperature[rIndex][phiIndex] += (1 + 0.01 * random.uniform(-1,1) )*simulationSettings['Heat'] * simulationSettings['dt']
				else:
					temperature[rIndex][phiIndex] += -(1 + 0.01 * random.uniform(-1,1))*simulationSettings['Heat'] * simulationSettings['dt'] * (44**2/(50**2 - 44**2))
	
		
		if (state == "Time Fail"):
			
			temperature = oldtemperature
			streamfunction = oldstreamfunction
			simulationSettings['dt'] = simulationSettings['dt']/2
			print("time halving, now has value", simulationSettings['dt'])
			continue 
		
		vorticity = updateVorticity(vorticity, temperature, streamfunction, simulationSettings)
		
		
		t = t + simulationSettings['dt']
		
		index += 1

main()
		
		

	
	
	
			


	
			
	
	
	
	
	
	
	

	
	
			
			

				
				
	
		
	
	
	

