"""
Fluids confined to a box with a thermal source at the bottom
"""

import matplotlib.pyplot as plt
import time
import copy
import math
import random
import numpy as np

"""
Lets change the boundry conditions so that in the x direction we have periodic boundries

so streamfunction[i][j]=streamfunction[i+N][j]
"""


def updateStreamFunction(vorticity, streamfunction, n, iterations, epsilon, omega, h):
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
		for i in range(n):
			for j in range(1,n-1):
				streamfunction[i][j] = 0.25*omega*(streamfunction[(i+1)%n][j]+streamfunction[(i-1)%n][j]+streamfunction[i][j+1]+streamfunction[i][j-1]+h*h*vorticity[i][j])+ (1.0-omega)*streamfunction[i][j]
				
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
		for j in range(0, n):
			if (j == 0):
				u = 0
				v = -(streamfunction[(i+1)%n][j] - streamfunction[(i-1)%n][j])/(2*h)
				
				if (abs(u)*dt/h >= 1 or abs(v)*dt/h >= 1):
					return temperature, 1
				if (u >=0):
					Tx = (temperature[(i)%n][j] - temperature[(i-1)%n][j])/h
				else:
					Tx = (temperature[(i+1)%n][j] - temperature[(i)%n][j])/h
					
				if (v >=0):
					Ty = (temperature[i][j] - temperature[i][1])/h
				else:
					Ty = (temperature[i][j+1] - temperature[i][j])/h
				
				term1 = -u*Tx
				term2 = -v*Ty
				term3 = kappa * (  (temperature[(i+1)%n][j] - 2*temperature[i][j] + temperature[(i-1)%n][j])/(h**2) + (temperature[i][j+1] - 2*temperature[i][j] + temperature[i][1])/(h**2)  )
			
			
			if (j == n-1):
				u = 0
				v = -(streamfunction[(i+1)%n][j] - streamfunction[(i-1)%n][j])/(2*h)
				if (abs(u)*dt/h >= 1 or abs(v)*dt/h >= 1):
					return temperature, 1
				if (u >=0):
					Tx = (temperature[(i)%n][j] - temperature[(i-1)%n][j])/h
				else:
					Tx = (temperature[(i+1)%n][j] - temperature[(i)%n][j])/h
					
				if (v >=0):
					Ty = (temperature[i][j] - temperature[i][j-1])/h
				else:
					Ty = (temperature[i][n-2] - temperature[i][j])/h
				
				term1 = -u*Tx
				term2 = -v*Ty
				term3 = kappa * (  (temperature[(i+1)%n][j] - 2*temperature[i][j] + temperature[(i-1)%n][j])/(h**2) + (temperature[i][n-2] - 2*temperature[i][j] + temperature[i][j-1])/(h**2)  )
				
			if (j!= 0 and j!= n-1):
				u = (streamfunction[i][j+1] - streamfunction[i][j-1])/(2*h)
				v = -(streamfunction[(i+1)%n][j] - streamfunction[(i-1)%n][j])/(2*h)
				if (abs(u)*dt/h >= 1 or abs(v)*dt/h >= 1):
					return temperature, 1
				##if (u*0.1/h >= 1 or v*0.1/h >= 1):
					##print(u,v)
					##print("unstable now")
				if (u >=0):
					Tx = (temperature[(i)%n][j] - temperature[(i-1)%n][j])/h
				else:
					Tx = (temperature[(i+1)%n][j] - temperature[(i)%n][j])/h
					
				if (v >=0):
					Ty = (temperature[i][j] - temperature[i][j-1])/h
				else:
					Ty = (temperature[i][j+1] - temperature[i][j])/h
				
				term1 = -u*Tx
				term2 = -v*Ty
				term3 = kappa * (  (temperature[(i+1)%n][j] - 2*temperature[i][j] + temperature[(i-1)%n][j])/(h**2) + (temperature[i][j+1] - 2*temperature[i][j] + temperature[i][j-1])/(h**2)  )
				
				
				
				
			term4 = heat[i][j]/Cp
			
			update[i][j] = term1 + term2 + term3 + term4
			
			newTemperature[i][j] = temperature[i][j] + dt*update[i][j]
			
			
	## what sort of boundry conditions should I do??
	

	
	## now I need to do the surface:
	## I will hold the surface at a constant temperature		
	return newTemperature, None
	
## now we need to write something that updates the vorticity equation for creeping flow

def updateVorticity(vorticity, temperature, streamfunction, nu, alpha, g, rho0, h, n, dt):
	newVorticity = [[0 for i in range(n)] for j in range(n)]
	update = [[0 for i in range(n)] for j in range(n)]	
	## Now we need to look at eh bounrdy conditions.
	
	for i in range(n):
		vorticity[i][0] = -2*streamfunction[i][1]/(h**2) 
		
		vorticity[i][n-1] = -2*streamfunction[i][n-2]/(h**2)
		
		newVorticity[i][0] = -2*streamfunction[i][1]/(h**2)
		
		newVorticity[i][n-1] = -2*streamfunction[i][n-2]/(h**2)
		
		

	for i in range(n):
		for j in range(1, n-1):
			update[i][j] = -g*alpha/rho0 * ( temperature[(i+1)%n][j] - temperature[(i-1)%n][j])/(2*h) + nu*( (vorticity[(i+1)%n][j] - 2*vorticity[i][j] + vorticity[(i-1)%n][j]  )/(h**2) + (vorticity[i][j+1] -2*vorticity[i][j] + vorticity[i][j-1])/(h**2) ) 
			
	
	for i in range(n):
		for j in range(1, n-1):
			newVorticity[i][j] = vorticity[i][j] + dt*update[i][j]
		
	return  newVorticity


	
	
def plotFields(temperature, vorticity, streamfunction, xMesh, yMesh, t, index):

	temperature = np.rot90(temperature)
	temperature = np.rot90(temperature)
	temperature = np.rot90(temperature)

	totalTemperature = 0
	for rIndex in range(len(temperature)):
		for thetaIndex in range(len(temperature[0])):
			totalTemperature += temperature[rIndex][thetaIndex]
	

	xPosition = 0
	yPosition = 0

	for xIndex in range(len(temperature)):
		for yIndex in range(len(temperature[0])):

			xPosition += xMesh[xIndex][yIndex]*temperature[xIndex][yIndex]/totalTemperature
			yPosition += yMesh[xIndex][yIndex]*temperature[xIndex][yIndex]/totalTemperature

	fig, ax = plt.subplots(dpi=120)
	pc = ax.pcolormesh(xMesh, yMesh, temperature, antialiased=False)
	x = [xPosition]
	y = [yPosition]
	print(t, index, x, y)
	plt.plot(x, y, 'or')
	plt.colorbar(pc)
	plt.xlabel("x (horizontal)")
	plt.ylabel("y (vertical)")
	plt.title(str(t))
	plt.savefig('t/'+str(index)+".png")
	plt.close()
	plt.clf()
	fig.clear()

	
	fig, ax = plt.subplots(dpi=120)
	vorticity = np.rot90(vorticity)
	vorticity = np.rot90(vorticity)
	vorticity = np.rot90(vorticity)
	pc = ax.pcolormesh(xMesh, yMesh, vorticity, antialiased=False)
	plt.colorbar(pc)
	##ax.contourf(xMesh, yMesh, vorticity, 100)
	plt.xlabel("x (horizontal)")
	plt.ylabel("y (vertical)")
	plt.savefig('vts/'+str(index)+".png")
	plt.close()
	plt.close()
	plt.clf()
	fig.clear()
	
	fig, ax = plt.subplots(dpi=120)
	streamfunction = np.rot90(streamfunction)
	streamfunction = np.rot90(streamfunction)
	streamfunction = np.rot90(streamfunction)
	pc = ax.pcolormesh(xMesh, yMesh, streamfunction, antialiased=False)
	plt.colorbar(pc)
	##ax.contourf(xMesh, yMesh, vorticity, 100)
	plt.xlabel("x (horizontal)")
	plt.ylabel("y (vertical)")
	plt.savefig('sfs/'+str(index)+".png")
	plt.close()
	plt.close()
	plt.clf()

	fig.clear()
	plt.close()
	
	
	
	
	
def getTotalQ(T):
	s = 0
	for row in T:
		for el in row:
			s += el
	return s
	
	
	
def main():

	## function that ties everything together
	
	scale = 1 ## length scale of the box in meters
	n = 50
	timesteps = 10000000
	dt = 0.05
	
	iterations = 1000
	epsilon = 0.0001
	omega = 1.5
	
	h = 1/(n-1) * scale
	Cp=4000
	kappa=0
	
	Touter=0
	rho0=1000
	##alpha=10**(-7)
	alpha=0
	nu=10**(-10)
	g=-10
	t = 0
	
	
	xs = np.arange(0, scale+h, h)
	ys = np.arange(0, scale+h, h)
	
	xMesh, yMesh = np.meshgrid(xs, ys)
	
	TotalQ = []
	times = []
	
	
		
	
	temperature = [[0 for i in range(n)] for j in range(n)]
	##for i in range(int(0.4*n), int(0.6*n)):
		##for j in range(int(0.4*n), int(0.6*n)):
			##temperature[i][j] = 1

	
	
	streamfunction = [[0 for i in range(n)] for j in range(n)]
	
	vorticity = [[0 for i in range(n)] for j in range(n)]
	
	heat = [[0 for i in range(n)] for j in range(n)]
	
	
	streamfunction1 = [[-0.0005*(i+j)/math.sqrt(2) for i in range(n)] for j in range(n)]
	streamfunction2 = [[0.0005*(i+j)/math.sqrt(2) for i in range(n)] for j in range(n)]
	
	interval = 6

	for xIndex in range(20, 30):
		for yIndex in range(20, 30):
			temperature[xIndex][yIndex] = 1

	for timestep in range(timesteps):


		if (t < interval):
			streamfunction = streamfunction1
		elif (t < 3*interval):
			print("back it up")
			streamfunction = streamfunction2
		elif (t < 4*interval):
			streamfunction = streamfunction1
		else:
			streamfunction = [[0 for i in range(n)] for j in range(n) ]
	

	
		oldstreamfunction = copy.deepcopy(streamfunction)
		oldtemperature = copy.deepcopy(temperature)
		oldvorticity = copy.deepcopy(vorticity)
		TotalQ.append(getTotalQ(temperature))
		times.append(t)
		
		##print("loop: " + str(timestep))

		if ( timestep % 1 == 0):
			plotFields(temperature, vorticity, streamfunction, xMesh, yMesh, t, timestep)
			plt.plot(times, TotalQ)
			plt.ylabel("Total Energy (arb. units)")
			plt.xlabel("time (arb units)")
			plt.savefig("totalq.png")
	
		##streamfunction = updateStreamFunction(vorticity, streamfunction, n, iterations, epsilon, omega, h)
		
		temperature, error = updateTemperature(streamfunction, temperature, heat, h, n, Cp, kappa, Touter,dt)

		if (error != None):
			streamfunction = oldstreamfunction
			oldtemperature = temperature
			dt = dt/2
			print("new dt", dt)
			continue
		

		
		
		##for i in range(n):
			##heat[i][1] = 100 + random.uniform(0,1)
			##heat[i][n-2] = -100 - random.uniform(0,1)
			
		
		vorticity = updateVorticity(vorticity, temperature, streamfunction, nu, alpha, g, rho0, h, n, dt)
		
		##print(vorticity)
		
		
		

		t = t + dt


			
	
		
		
		
		
main()

		

	
	
	
			


	
			
	
	
	
	
	
	
	

	
	
			
			

				
				
	
		
	
	
	

