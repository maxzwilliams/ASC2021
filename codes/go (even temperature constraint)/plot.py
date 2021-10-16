"""
Python program for plotting the output of the go LB code 
"""

import math
import csv
import sys

import copy
import matplotlib.pyplot as plt
import csv


def convertToArray(fileName):
	file = open(fileName)
	csvr = csv.reader(file)
	rows = []
	for row in csvr:
		rows.append(row)
	return rows
	
def convertStringArrayToFloatArray(stringArray2D):

	rtn = copy.deepcopy(stringArray2D)
	
	for rowIndex in range(len(stringArray2D)):
		for colIndex in range(len(stringArray2D[rowIndex])):
			rtn[rowIndex][colIndex] = float(stringArray2D[rowIndex][colIndex])
			
	return rtn
	
	

	
	
def goPlot(data, vmin, vmax, name):
	print("plotting")
	fig, ax = plt.subplots(dpi=120)
	pc = ax.pcolorfast(data, vmin=vmin, vmax=vmax)
	plt.colorbar(pc)
	plt.xlabel("x")
	plt.ylabel("y")
	plt.title(name)
	plt.savefig(name+".png")
	print("done plotting")
	plt.close()
	
	

	

if __name__ == "__main__":
	plotsMax = int(sys.argv[2])
	plotsMin = int(sys.argv[1])
	stepSize = int(sys.argv[3])
	
	fileSets = []
	for index in range(plotsMin,plotsMax+1, stepSize):
		fileSets.append(["u//ux"+str(index)+".csv", "u//uy"+str(index)+".csv", "rho//rho"+str(index)+".csv", "ep//ep"+str(index)+".csv"  ])
		##fileSets.append(["u//ux"+str(index)+".csv", "u//uy"+str(index)+".csv", "rho//rho"+str(index)+".csv"])

		
	
	
	for el in fileSets:
		uxData = convertToArray(el[0])
		uxData = convertStringArrayToFloatArray(uxData)
		uyData = convertToArray(el[1])
		uyData = convertStringArrayToFloatArray(uyData)
		rhoData = convertToArray(el[2])
		rhoData = convertStringArrayToFloatArray(rhoData)


		
		goPlot(uxData, None, None, "p"+el[0] )
		goPlot(uyData, None, None, "p"+el[1])
		goPlot(rhoData, None, None,"p"+el[2])

		epData = convertToArray(el[3])
		epData = convertStringArrayToFloatArray(epData)
		goPlot(epData, None, None,"p"+el[3])
		
