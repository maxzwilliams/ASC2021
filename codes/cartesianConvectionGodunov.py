"""
Bousinesque fluid confined between plates from the top and bottom with left and right boundaries periodic.
"""

import matplotlib.pyplot as plt
import time
import copy
import math
import random
import numpy as np


def updateStreamFunction(vorticity, streamfunction, n, iterations, epsilon, omega, h):
	"""
	Function that finds the streamfunction using the SOR method in cartesian coordiantes
	
	Parameters:
		vorticity (n lists of size n): vorticity at each point in the domein
		streamfunction (n lists of size n): streamfunction at each point in the domain
		n (int): size of the domain
		iterations (int): number of iterations used for the SOR solver
		epsilon (float): acceptable error in streamfunction solution
	 	omega (float): SOR parameter
	 	h (float): mesh size (distance between adjacent mesh points)
	 
	 Returns:
	 	streamfunction: new streamfunction
	"""
	
	## Begin to iteratively solve the poisson equation
	for iteration in range(iterations):
		oldstreamfunction = copy.deepcopy(streamfunction) ## save the current streamfunction solution

		## apply SOR formula for all interal points in the domain
		for i in range(n):
			for j in range(1,n-1):
				streamfunction[i][j] = 0.25*omega*(streamfunction[(i+1)%n][j]+streamfunction[(i-1)%n][j]+streamfunction[i][j+1]+streamfunction[i][j-1]+h*h*vorticity[i][j])+ (1.0-omega)*streamfunction[i][j]
		
		## get error (difference between current and previous iterations streamfunction)
		epsilonDash = 0
		for i in range(n):
			for j in range(n):
				epsilonDash += abs(oldstreamfunction[i][j] - streamfunction[i][j])
		
		## If error is suffeciently small, then we say the solution has converged and we break the loop and return
		if (epsilonDash < epsilon):
			break;
		## if we have run out of iterations and not converged, we panic and stop the program.
		if (iteration == iterations-1):
			raise Exception("when solving for streamfunction the method did not converge")

	## return our converged streamfunction
	return streamfunction
	
	
def updateTemperature(streamfunction, temperature, heat, h, n, kappa, dt):
	"""
	Function that updates the temperature field

	Parameters:
		streamfunction (n lists of size n): streamfunction at each point in the domain
		temperature (n lists of size n): temperature at each point in the domain
		heat (n lists of size n): heat source at each point
		h (float): mesh size (distance between adjacent mesh points)
		n (int): size of the domain
		kappa (float): thermal diffusivity
		dt (float): timestep 
	 	
	returns:
		newTemperature (n lists of size n): updated temperature at each location
		error (either None or 1): is 1 if conditions for an instability occured, else None
	"""
	

	update = [[0 for i in range(n)] for j in range(n)] ## term that will be added to temperature to update it
	newTemperature = [[0 for i in range(n)] for j in range(n)] ## the new temperature
	## for each position
	for i in range(n): 
		for j in range(0, n):

			## if at bottom boundary we have to handle velocities and boundary temperatures differently
			if (j == 0):
				## in this if statement, we have to use temperature[i][1] rather than temperature[i][j-1]
				## set velocities
				u = 0 
				v = -(streamfunction[(i+1)%n][j] - streamfunction[(i-1)%n][j])/(2*h)
				## check for stability conditions
				if (abs(u)*dt/h >= 1 or abs(v)*dt/h >= 1):
					return temperature, 1 ## if instabiltiy is detected return whatever the current temperature is and 1 denoting an error.

				## apply godanonv scheme
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
			
			## if at top boundary we have to handle velocities and boundary temperatures differently
			if (j == n-1):
				## in this if statement, we have to use temperature[i][n-2] rather than temperature[i][j+1]
				## set velocities
				u = 0
				v = -(streamfunction[(i+1)%n][j] - streamfunction[(i-1)%n][j])/(2*h)
				## check for stability conditions
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
			
			## if not at the top or bottom boundaries 
			if (j!= 0 and j!= n-1):
				u = (streamfunction[i][j+1] - streamfunction[i][j-1])/(2*h)
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
					Ty = (temperature[i][j+1] - temperature[i][j])/h
				
				term1 = -u*Tx
				term2 = -v*Ty
				term3 = kappa * (  (temperature[(i+1)%n][j] - 2*temperature[i][j] + temperature[(i-1)%n][j])/(h**2) + (temperature[i][j+1] - 2*temperature[i][j] + temperature[i][j-1])/(h**2)  )
				
				
				
				
			term4 = heat[i][j] 
			
			update[i][j] = term1 + term2 + term3 + term4 ## generate the update from the computations above
			
			newTemperature[i][j] = temperature[i][j] + dt*update[i][j] ## set the new temperature
			
	return newTemperature, None ## return the new temperature and no error as non occured.
	


def updateVorticity(vorticity, temperature, streamfunction, nu, alpha, g, rho0, h, n, dt):
	"""
	update the vorticity at every point

	Parameters:
		vorticity (n lists of size n): vorticity at each point in the domein
		temperature (n lists of size n): temperature at each point in the domain
		streamfunction (n lists of size n): streamfunction at each point in the domain
		nu (float): kinematic viscousity
		alpha (float): Thermal expansion coeffecient
		g (float): gravitational field
		rho0 (float64): reference density
		h (float): mesh size (distance between adjacent mesh points)


		n (int): size of the domain
		dt (float): timestep 

	Returns:
		newVorticity (n lists of size n): the updated vorticity

	"""
	newVorticity = [[0 for i in range(n)] for j in range(n)] ## make an array to store the new vorticity in
	update = [[0 for i in range(n)] for j in range(n)]	## update array to add to vortivity to make newVorticity
	
	## enforce boundary conditions on vorticity
	for i in range(n):
		vorticity[i][0] = -2*streamfunction[i][1]/(h**2) 
		
		vorticity[i][n-1] = -2*streamfunction[i][n-2]/(h**2)
		
		newVorticity[i][0] = -2*streamfunction[i][1]/(h**2)
		
		newVorticity[i][n-1] = -2*streamfunction[i][n-2]/(h**2)
		
		
	## for all internal points get the update part for the vortivity
	for i in range(n):
		for j in range(1, n-1):
			update[i][j] = -g*alpha/rho0 * ( temperature[(i+1)%n][j] - temperature[(i-1)%n][j])/(2*h) + nu*( (vorticity[(i+1)%n][j] - 2*vorticity[i][j] + vorticity[(i-1)%n][j]  )/(h**2) + (vorticity[i][j+1] -2*vorticity[i][j] + vorticity[i][j-1])/(h**2) ) 
			
	
	## sum the current vorticity and the update together
	for i in range(n):
		for j in range(1, n-1):
			newVorticity[i][j] = vorticity[i][j] + dt*update[i][j]
	
	## return the new vorticity
	return  newVorticity
	

"""
BEGIN HELPER FUNCTION SECTION
"""
def getTotalQ(T):
	s = 0
	for row in T:
		for el in row:
			s += el
	return s


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
	##print(t, index, x, y)
	##plt.plot(x, y, 'or')
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
		

"""
END HELPER FUNCTION SECTION
"""
	
def main():
	print("Program Running!!")
	## function that ties everything together
	
	scale = 1 ## length scale of the box in meters
	n = 50 ## number of compuational points along the horizontal and vertical directions
	timesteps = 10000000 ## total number of timesteps used
	dt = 0.01 ## length of each timestep

	## parameters for Jacobi method with SOR used for solving the streamfunction
	## these have nothing to do with the physics, they are purely mathematical	
	iterations = 1000 ## iterations for the poisson solver
	epsilon = 0.0001/35**2 * n**2 ## acceptable error for the poisson solver
	omega = 1.5 ## SOR parameter
	
	## distnace between adjacent compuational nodes
	h = 1/(n-1) * scale 
	##Cp=4000 
	kappa=0.0000000001 ## thermal diffusivity
	
	
	rho0=1000 ## reference density
	##alpha=10**(-7) 
	alpha=1 ## thermal expansion coeffecient
	nu=10**(-10) ## kinematic viscousity
	g=-100 ## graviational acceleration 
	t = 0 ## current time
	Heat = 1000 ## heating field 
	
	## the following three lines are useful only for faster plotting, not to do with the main algorithm 
	xs = np.arange(0, scale+h, h) ## all horizontal position values for compuational nodes 
	ys = np.arange(0, scale+h, h) ## all vertical position values for compuational nodes
	xMesh, yMesh = np.meshgrid(xs, ys) ## making meshes for fast plotting
	

	TotalQ = [] ## keeping track of the total energy
	times = [] ## times when we looked at the total energy 
	

	Pr = nu/kappa ## Prantl number
	Ra = g * alpha * scale**5 * Heat/(nu * kappa**2) ## Rayleigh number

	print("")
	print("Pr:",Pr, "Ra:", Ra)
	print("")
		
	## setup arrays for all the important properties
	temperature = [[0 for i in range(n)] for j in range(n)] 
	streamfunction = [[0 for i in range(n)] for j in range(n)]
	vorticity = [[0 for i in range(n)] for j in range(n)]
	heat = [[0 for i in range(n)] for j in range(n)]



	for xIndex in range(n):
		for yIndex in range(n):
			if (yIndex < 3):
				heat[xIndex][yIndex] = Heat + 0.00001*(2.0*random.uniform(0,1) -1.0)
			if (yIndex > n-3):
				heat[xIndex][yIndex] = -Heat + 0.00001*(2.0*random.uniform(0,1) -1.0)


	## iterate over time
	for timestep in range(timesteps):

		print("timestep:", timestep, end="\r")



	

		## keep track of the fundamental properties incase we need to restart a loop
		## due to instability conditions being met
		oldstreamfunction = copy.deepcopy(streamfunction)
		oldtemperature = copy.deepcopy(temperature)
		oldvorticity = copy.deepcopy(vorticity)

		## keep track of the total energy 
		TotalQ.append(getTotalQ(temperature))
		times.append(t)

		## every 100 timesteps we plot and save
		if ( timestep % 100 == 0):
			plotFields(temperature, vorticity, streamfunction, xMesh, yMesh, t, timestep)
			plt.plot(times, TotalQ)
			plt.ylabel("Total Energy (arb. units)")
			plt.xlabel("time (arb units)")
			plt.savefig("totalq.png")
		
		## updat the streamfunction
		streamfunction = updateStreamFunction(vorticity, streamfunction, n, iterations, epsilon, omega, h)
		## update the temperature
		temperature, error = updateTemperature(streamfunction, temperature, heat, h, n, kappa, dt)



		## apply heating conditions to the temperature
		for xIndex in range(n):
			for yIndex in range(n):
				temperature[xIndex][yIndex] += heat[xIndex][yIndex] * dt

		if (error != None): ## if we got an error due to instability conditions from updating temperature
			## reset everything back one timestep 
			streamfunction = oldstreamfunction 
			oldtemperature = temperature
			dt = dt/2 ## decrease the timestep length to get into a stable condition
			print("new dt", dt) ## tell the terminal we are doing this
			continue ## go back to the top of the loop and redo the previous timestep
			
		## update the vortcity
		vorticity = updateVorticity(vorticity, temperature, streamfunction, nu, alpha, g, rho0, h, n, dt)


		t = t + dt ## keep track of time
	
		
main() ## run the main loop when the program is run

		

	
	
	
			


	
			
	
	
	
	
	
	
	

	
	
			
			

				
				
	
		
	
	
	

