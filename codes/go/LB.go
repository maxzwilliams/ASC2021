// LBM for natural convection in the earth

package main
import "fmt"
import "math"
import "math/rand"
import "os"
import "encoding/csv"
import "strconv"
import "time"




type D2Field struct {
	entries [][][]float64 // orientation, x coorindate, y coorindate
	name string // name of the field
}

type D3Field struct {
	entries [][][][]float64 // orientation, x coordinate, y coorindate, z coordinate
	name string // name of the field
}

type D2Property struct{
	entries [][]float64 // a number for each point in (x,y) space
	name string // name of the property
}


type lattice struct{ // defines the lattice we are solving the problem on
	weights []float64 // weights for different directions in the lattice
	directions [][]int // direction vectors in the lattice
	opposites []int // pointer array. At index i, directions[opposites[i]] points opposite to directions[i]
	vels []float64 // velocity constants for finding equilibrium
	eqWeights []float64 // 
	nx int // number of latice points in x direction
	ny int // number of lattice points in y direction
	Q int // directions in the lattice following the D2Q9 style convection
	D int // dimension of the lattice
	dt float64 // timestep, here it is always assumed 1 (dont change this, the code will break)
	dx float64 // space timestep
	S float64 // lattice speed
	tau_f float64 // density distribution function relaxation time
	tau_g float64 // energy distribution function relaxation time
	tau_grav float64 // relaxation time for gravity term
	rho float64 // reference density
	alpha float64 // thermal expansion coeffecient

	Cv float64 // specific heat at constant volume
	Heat float64 // heat generated per unit volume

}


/* A boundary function. It returns true if you are in a solid
and false if you are in a fluid. Takes arguments of the coordinates
you are currently at.
*/
type boundary func(cord1 int, cord2 int) bool 


func streamingStep2D(f1 D2Field, l lattice, b boundary) D2Field{
	/*
	Function that performs the streaming step on a field f1 on a lattice l given a boundary b

	Parameters:
		f1 (D2Field): Field that will be streamed
		l (lattice): lattice that will be used for the streaming step
		b (boundary): Boundary function that is used for simulating 

	Returns:
		(D2Field): updated field f1 streamed on timestep
	*/

	// generate an empty field to return
	var f2 D2Field
	entries := D3Constant(l.nx, l.ny, l.Q, 0.0)
	f2.entries = entries 
	f2.name = f1.name

	for dIndex, d := range l.directions{ // for each direction in the lattice
	for i:=0;i<l.nx;i++{ // for each location in the lattice
		for j:=0;j<l.ny;j++{

				// if you are streaming to a solid location
				if (b(i-d[0],j-d[1])){
					// reflect the streaming step backwards
					f2.entries[i][j][dIndex] = f1.entries[i][j][l.opposites[dIndex]]
				} else{ // if you are streaming within the fluid
					// stream each direction forward in time
					f2.entries[i][j][dIndex] = f1.entries[i-d[0]][j-d[1]][dIndex]
				}
			}
		}
	}
	return f2 // return our updated streamed field
}

// old function, dont use
func getMacro2D(f D2Field, l lattice, b boundary) (D2Property, D2Field){
	var rho D2Property
	var u D2Field
	
	rhoentries := D2Constant(l.nx, l.ny, 0.0)
	uentries := D3Constant(l.nx, l.ny, l.D, 0.0)
	
	rho.entries = rhoentries
	u.entries = uentries

	for xIndex:=0;xIndex<l.nx;xIndex++{
		for yIndex:=0;yIndex<l.ny;yIndex++{
			if (!b(xIndex, yIndex)){
				for dIndex, _ := range l.directions{
					rho.entries[xIndex][yIndex] += f.entries[xIndex][yIndex][dIndex]
				}
				for dIndex, d := range l.directions{

					u.entries[xIndex][yIndex][0] += float64(d[0]) * f.entries[xIndex][yIndex][dIndex]/rhoentries[xIndex][yIndex]
					u.entries[xIndex][yIndex][1] += float64(d[1]) * f.entries[xIndex][yIndex][dIndex]/rhoentries[xIndex][yIndex]
				}
			}
		}
	}
	return rho, u // rho and u are zero in the dead nodes
}

func getMacro2DWithEp(f D2Field, g D2Field, l lattice, b boundary) (D2Property, D2Property, D2Field) {
	/*
	This function finds the macroscopic properties (density, internal energy and velocity) 
	from the particle and energy distribution functions

	Parameters:
		f (D2Field): Particle distribution function
		g (D2Field): Energy distribution function
		l (lattice): lattice that we are simulating on
		b (boundary): liquid-solid boundary for this problem

	Returns:
		(D2Property): Macroscopic density
		(D2Property): Macroscopic internal energy density
		(D2Field): Macroscopic velocity field
	*/

	// define return variables
	var rho D2Property
	var ep D2Property
	var u D2Field
	
	// create blank entries for these return variables
	rhoentries := D2Constant(l.nx, l.ny, 0.0)
	epentries := D2Constant(l.nx, l.ny, 0.0)
	uentries := D3Constant(l.nx, l.ny, l.D, 0.0)
	// Assign blank entries
	rho.entries = rhoentries
	ep.entries = epentries
	u.entries = uentries
	
	// for each position in the lattice
	for xIndex:=0;xIndex<l.nx;xIndex++{
		for yIndex:=0;yIndex<l.ny;yIndex++{

			if (!b(xIndex, yIndex)){ // If the current position is a liquid

				// get density at this point
				for dIndex, _ := range l.directions{ 
					rho.entries[xIndex][yIndex] += f.entries[xIndex][yIndex][dIndex]
				}

				// perform similar sums to get the other macroscopic properties
				for dIndex, d := range l.directions{
					u.entries[xIndex][yIndex][0] += float64(d[0]) * f.entries[xIndex][yIndex][dIndex]/rho.entries[xIndex][yIndex]
					u.entries[xIndex][yIndex][1] += float64(d[1]) * f.entries[xIndex][yIndex][dIndex]/rho.entries[xIndex][yIndex]
					ep.entries[xIndex][yIndex] += g.entries[xIndex][yIndex][dIndex]/rho.entries[xIndex][yIndex]
				}
			}
		}
	}
	return rho, ep, u // return our macroscopic properties
}

//old function
func getEq2D(rho D2Property, u D2Field, l lattice) D2Field{
	// this looks good
	
	var f D2Field
	entries := D3Constant(l.nx, l.ny, l.Q, 0.0)
	f.entries = entries
	f.name = "f"

	var term1 float64
	var term2 float64
	var term3 float64
	var term4 float64
	
	for dIndex, d := range l.directions{
		for xIndex:=0;xIndex<l.nx;xIndex++{
			for yIndex:=0;yIndex<l.ny;yIndex++{
			
				term1 = l.vels[0]
				term2 = l.vels[1] * ( float64(d[0])*u.entries[xIndex][yIndex][0] + float64(d[1])*u.entries[xIndex][yIndex][1] )
				term3 = l.vels[2] * math.Pow(( float64(d[0])*u.entries[xIndex][yIndex][0] + float64(d[1])*u.entries[xIndex][yIndex][1] ), 2.0)
				term4 = l.vels[3] * (u.entries[xIndex][yIndex][0]*u.entries[xIndex][yIndex][0] + u.entries[xIndex][yIndex][1]*u.entries[xIndex][yIndex][1]   )
				
				f.entries[xIndex][yIndex][dIndex] = rho.entries[xIndex][yIndex] * l.weights[dIndex]*(term1 + term2 + term3 + term4)
			}
		}
	}
	return f
}


func getEq2Df(rho D2Property, u D2Field, l lattice) D2Field{
	/*
	Function that gets equilibrium particle distribution function

	Parameters:
		rho (D2Property): Macroscopic density 
		u (D2Field): Macroscopic velocity
		l (lattice): lattice we are performing our simulation on

	Returns:
		(D2Field): the equilibrium particle distribution function
	*/
	
	// generate an empty field to return
	var f D2Field
	entries := D3Constant(l.nx, l.ny, l.Q, 0.0)
	f.entries = entries
	f.name = "f"
	
	
	// define variables once rather in loop for faster compute
	var term1 float64
	var term2 float64
	var term3 float64
	var term4 float64
	
	// for each lattice direction for all lattice locations
	for dIndex, d := range l.directions{
		for xIndex:=0;xIndex<l.nx;xIndex++{
			for yIndex:=0;yIndex<l.ny;yIndex++{
				
				// calculate terms in f^{eq}
				term1 = l.vels[0]
				term2 = l.vels[1] * ( float64(d[0])*u.entries[xIndex][yIndex][0] + float64(d[1])*u.entries[xIndex][yIndex][1] )
				term3 = l.vels[2] * math.Pow(( float64(d[0])*u.entries[xIndex][yIndex][0] + float64(d[1])*u.entries[xIndex][yIndex][1] ), 2.0)
				term4 = l.vels[3] * (u.entries[xIndex][yIndex][0]*u.entries[xIndex][yIndex][0] + u.entries[xIndex][yIndex][1]*u.entries[xIndex][yIndex][1]   )
				
				// assign f^{eq} according to LB theory
				f.entries[xIndex][yIndex][dIndex] = rho.entries[xIndex][yIndex] * l.weights[dIndex] * (term1 + term2 + term3 + term4)
			}
		}
	}
	return f // return the equilibrium particle distribution function
}

func getEq2Dg(rho D2Property, ep D2Property, u D2Field, l lattice) D2Field{
	/*
	Function that gets equilibrium energy distribution function

	Parameters:
		rho (D2Property): Macroscopic density 
		ep (D2Property): Macroscopic internal energy density
		u (D2Field): Macroscopic velocity
		l (lattice): lattice we are performing our simulation on

	Returns:
		(D2Field): the equilibrium energy distribution function
	*/
	
	// same algorithm is used to getEqDf
	var g D2Field
	entries := D3Constant(l.nx, l.ny, l.Q, 0.0)
	g.entries = entries
	g.name = "g"
	
	
	
	var term1 float64
	var term2 float64
	var term3 float64
	var term4 float64
	
	for dIndex, d := range l.directions{
		for xIndex:=0;xIndex<l.nx;xIndex++{
			for yIndex:=0;yIndex<l.ny;yIndex++{
			
				term1 = l.vels[0]
				term2 = l.vels[1] * ( float64(d[0])*u.entries[xIndex][yIndex][0] + float64(d[1])*u.entries[xIndex][yIndex][1] )
				term3 = l.vels[2] * math.Pow(( float64(d[0])*u.entries[xIndex][yIndex][0] + float64(d[1])*u.entries[xIndex][yIndex][1] ), 2.0)
				term4 = l.vels[3] * (u.entries[xIndex][yIndex][0]*u.entries[xIndex][yIndex][0] + u.entries[xIndex][yIndex][1]*u.entries[xIndex][yIndex][1]   )
				
				g.entries[xIndex][yIndex][dIndex] = ep.entries[xIndex][yIndex]* rho.entries[xIndex][yIndex] * l.weights[dIndex]*(term1 + term2 + term3 + term4)
			}
		}
	}
	return g 
}


func collisionStep(f D2Field, feq D2Field, l lattice) D2Field{
	/*
	Performs the collision step on a field given its current value, its equilibrium value and the lattice

	Parameters:
		f (D2Field): field we are performing the collision step on
		feq (D2Field): equilibrium value for the field
		l (lattice): lattice we are performing the simulation on

	Returns:
		(D2Field): f after collision. 
	*/
	// empty field to return
	var rtn D2Field
	entries := D3Constant(l.nx, l.ny, l.Q, 0.0)
	rtn.entries = entries
	rtn.name = f.name
	
	var tau float64 // relaxation time for the field
	// change relaxation time depending on the type of field
	if ( f.name == "f"){
		tau = l.tau_f
	} else if (f.name == "g"){
		tau = l.tau_g
	} else {
		panic("unknown field being used")
	}
	// for each direction and position in the lattice
	for dIndex, _ := range l.directions{
		for xIndex:=0;xIndex<l.nx;xIndex++{
			for yIndex:=0;yIndex<l.ny;yIndex++{
				// apply the collision term
				rtn.entries[xIndex][yIndex][dIndex] = f.entries[xIndex][yIndex][dIndex] +  1.0/(tau) *  ( feq.entries[xIndex][yIndex][dIndex] - f.entries[xIndex][yIndex][dIndex] )
			}
		}
	}
	return rtn // return updated field with collision step
}


func gravityAffect(f D2Field, grav D2Field, rho D2Property, ep D2Property, l lattice) D2Field{
	/*
	Applys gravity to a field f. It assumes that f is a particle distribution function

	Parameters:
		f (D2Field): particle distribution function
		grav (D2FIeld): gravity field 
		rho (D2Property): macroscopic density
		ep (D2Property): macroscopic internal energy density
		l (lattice): lattice which we are solving the problem on

	Returns:
		(D2Field): updated field f with gravity term added.	
	*/
	
	// make a return field
	var rtn D2Field
	entries := D3Constant(l.nx, l.ny, l.Q, 0.0)
	rtn.entries = entries
	rtn.name = "f"
	
	// for each direction and position in the lattice
	for dIndex, d := range l.directions{
		for xIndex:=0;xIndex<l.nx;xIndex++{
			for yIndex:=0;yIndex<l.ny;yIndex++{
				// if null direction, then no gravity affect.
				if (dIndex == 0){
					rtn.entries[xIndex][yIndex][dIndex] = f.entries[xIndex][yIndex][dIndex] + 0.0
				
				} else{
					// update f according to theory from LBM
					term := -1.0/(l.tau_grav) * l.weights[dIndex] * rho.entries[xIndex][yIndex] * ( ep.entries[xIndex][yIndex]*2.0/float64(l.D) ) * l.alpha * math.Pow( math.Pow(float64(d[0]), 2.0) + math.Pow(float64(d[1]), 2.0), -0.5) * (grav.entries[xIndex][yIndex][0]*float64(d[0]) + grav.entries[xIndex][yIndex][1]*float64(d[1]))
					rtn.entries[xIndex][yIndex][dIndex] = f.entries[xIndex][yIndex][dIndex] + term
				}
			}
		}
	}
	return rtn // return our updated particle distribution function
}


// Radius 50 boundary function for a 102 by 102 lattice
func circleBoundaryR50(xIndex int, yIndex int) bool{	
	if (  math.Pow( float64(xIndex - 51), 2.0) + math.Pow(float64(yIndex - 51), 2.0) < math.Pow(float64(50), 2.0) ){
		return false
	}
	return true
}

/*
BEGIN HELPER FUNCTION SECTION
*/
func D3Constant(size1 int, size2 int, size3 int, c float64) [][][]float64{
	f := make([][][]float64,0)
	
	for i:=0; i<size1; i++{
		finner := make([][]float64,0)
		for j:=0; j<size2; j++{
			
			finnerinner := make([]float64, 0)
			
			for k:=0;k<size3;k++{
				finnerinner = append(finnerinner, c)
			}
			
			finner = append(finner, finnerinner)
		}
		
		f = append(f, finner)
	}		
	return f
}


func D2Constant(size1 int, size2 int, c float64) [][]float64{
	rtn := make([][]float64, 0)
	for i:=0;i<size1;i++{
		inner := make([]float64, 0)
		for j:=0;j<size2;j++{
			inner = append(inner, c)
		}
		rtn = append(rtn, inner)
	}
	return rtn 
}


func convertFloatArrayToStringArray(floatArray []float64) []string{
	rtn := make([]string, 0)
	for _, el := range(floatArray){
		rtn = append(rtn, fmt.Sprintf("%f", el))	
	}
	return rtn
}
    
func writeArrayCSV(floatArray [][]float64, fileName string){
    file, _ := os.Create(fileName)
    defer file.Close()

    writer := csv.NewWriter(file)
    defer writer.Flush()

    for _, value := range floatArray {
    	stringValue := convertFloatArrayToStringArray(value)
        writer.Write(stringValue)
    }
}


func testArrayNaN2D(array [][]float64) bool{
	for _, row := range array{
		for _, el := range row{
			if (math.IsNaN(el)){
				return true
			}
		}	
	}
	return false
}

func testArrayNaN3D(array [][][]float64) bool{
	for _, outer := range array{
		for _, row := range outer{
			for _, el := range row{
				if (math.IsNaN(el)){
					return true
				}
			}
		}	
	}
	return false
}

/*
END HELPER FUNCTION SECTION
*/

// main function which is run when program is run
func main(){
	fmt.Println("Program Running")

	// generate a seed for random numbers
	s1 := rand.NewSource(time.Now().UnixNano())
	rand.New(s1)

	var D2Q9 lattice // define a lattice
	
	var dt, dx float64 = 1, 1 // time and space discritzation
	var S float64 = dx/dt // lattice speed
	
	Q := 9 // directions in the lattice
	D := 2 // dimensions of the lattice
	
	
	
	// direction weights
	ws := []float64{4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0}
	

	v1 := []int{0,0} 
	v2 := []int{1,0} 
	v3 := []int{-1,0} 
	v4 := []int{0,1} 
	v5 := []int{0,-1} 
	v6 := []int{1,1} 
	v7 := []int{-1,-1} 
	v8 := []int{1,-1} 
	v9 := []int{-1,1} 
	// directions within the lattice
	cs := [][]int{v1, v2, v3, v4, v5, v6, v7, v8, v9}
	
	// pointer array. At index i, directions[opposites[i]] points opposite to directions[i]. This is the opposites field for D2Q9
	op := []int{0, 2, 1, 4, 3, 6, 5, 8, 7}
	
	// constants used in finding equilibrium distributions
	vels := []float64{ 1.0, 3.0/(math.Pow(S, 2.0)), 9.0/(2.0*math.Pow(S,4.0)), -3.0/(2.0*math.Pow(S, 2.0))}

	var vf float64 = 0.125 // kinematic viscosity
	var tau_f float64 = vf/(math.Pow(S/math.Sqrt(3), 2.0) * dt) + 0.5 // relaxation time for particle distribution 
	
	var vg float64 = 0.000025 // a measure of thermal diffusivity
	var tau_g float64 = vg/(math.Pow(S/math.Sqrt(3), 2.0) * dt  ) + 0.5 // relaxation time associated with thermal diffusion (associated with energy distribution function)

	var vgrav float64 = 0.6 
	var tau_grav float64 = vgrav/(math.Pow(S/math.Sqrt(3), 2.0) * dt  ) + 0.5 // relaxation time for gravity
	
	
	var rhoBack float64 = 1000.0 // background density
	var alpha float64 = 2.5*math.Pow(10.0, -11.0) // thermal expansion coeffecient
	
	// define the discritization
	nt := 10000000 // timesteps
	nx := 103 // steps horizontally
	ny := 103 // steps vertically
	
	// Assign all values to the lattice
	D2Q9.weights = ws
	D2Q9.directions = cs
	D2Q9.opposites = op
	D2Q9.vels = vels
	D2Q9.nx = nx
	D2Q9.ny = ny
	D2Q9.Q = Q
	D2Q9.D = D
	D2Q9.dt = dt
	D2Q9.dx = dx
	D2Q9.S = S
	//D2Q9.Cv = 400
	D2Q9.Heat = 1 
	D2Q9.rho = rhoBack
	D2Q9.alpha = alpha
	D2Q9.tau_f = tau_f
	D2Q9.tau_g = tau_g
	D2Q9.tau_grav = tau_grav



	maxG := -1.0 // define maximum gravity (gravity at exterior of boundary)
	// Get important dimensionless parameters and print them out
	Ra := math.Abs(maxG*D2Q9.alpha*D2Q9.Heat*math.Pow(50.0, 5.0)/(vf * math.Pow(vg, 2.0) ))
	Pr := vf/vg
	timeScaleDiff := math.Pow(50.0, 2.0)/vg
	timeScaleAdv := (vg*vf)/(math.Abs(maxG) * D2Q9.alpha * D2Q9.Heat * math.Pow(50.0, 3.0))
	fmt.Println("Ra: ", Ra)
	fmt.Println("Pr:", Pr)
	fmt.Println("Characteristic Diffusion TimeScale:", timeScaleDiff)
	fmt.Println("Characteristic Advection TimeScale:", timeScaleAdv)
	

	
	// then we define our field and properties that we will update as we evolve
	var f D2Field
	var g D2Field
	var fstream D2Field
	var gstream D2Field
	var feq D2Field
	var geq D2Field
	var rho D2Property
	var ep D2Property
	var u D2Field
	var gravity D2Field

	// populate these
	f.entries = D3Constant(nx, ny, Q, 0.0)
	f.name = "f" 
	g.entries = D3Constant(nx, ny, Q, 0.0)
	g.name = "g"
	fstream.entries =  D3Constant(nx, ny, Q, 0.0)
	fstream.name = "f"
	gstream.entries =  D3Constant(nx, ny, Q, 0.0)
	gstream.name = "g"
	feq.entries = D3Constant(nx, ny, Q, 0.0)
	feq.name = "f"
	geq.entries = D3Constant(nx, ny, Q, 0.0)
	geq.name = "g"
	rho.entries = D2Constant(nx, ny, 1.0)
	ep.entries = D2Constant(nx, ny, 0.0)
	u.entries = D3Constant(nx, ny, D, 0.0)
	
	// Define a gravity field
	gravity.name = "gravity"
	gravity.entries = D3Constant(nx, ny, 2, 0.0)
	gravityx := D2Constant(nx, ny, 0)
	gravityy := D2Constant(nx, ny, 0)
	gravityStrength := D2Constant(nx, ny, 0)
	// Set a gravity field for a self gravitating fluid
	for xIndex:=0;xIndex<D2Q9.nx;xIndex++{
		for yIndex:=0;yIndex<D2Q9.ny;yIndex++{
			xCord := float64(xIndex - 51)
			yCord := float64(yIndex - 51)

			theta := 0.0
			if (xCord == 0 && yCord > 0){
				theta = math.Pi/2
			} else if (xCord == 0 && yCord <= 0){
				theta = 3*math.Pi/2
			} else if (xCord > 0 && yCord >=0){
				theta = math.Atan(math.Abs(yCord/xCord))
			} else if (xCord <= 0 && yCord >= 0){
				theta = math.Pi - math.Atan(math.Abs(yCord/xCord))
			} else if (xCord <= 0 && yCord < 0){
				theta = math.Pi + math.Atan(math.Abs(yCord/xCord))
			} else if (xCord > 0 && yCord < 0){
				theta = math.Pi*2 - math.Atan(math.Abs(yCord/xCord))
			}

			gravity.entries[xIndex][yIndex][0] = maxG*math.Sqrt( math.Pow(xCord, 2.0) + math.Pow(yCord, 2.0))/50.0 * math.Cos(theta)
			gravity.entries[xIndex][yIndex][1] = maxG*math.Sqrt( math.Pow(xCord, 2.0) + math.Pow(yCord, 2.0))/50.0 * math.Sin(theta)
			gravityx[xIndex][yIndex] = gravity.entries[xIndex][yIndex][0]
			gravityy[xIndex][yIndex] = gravity.entries[xIndex][yIndex][1]
			gravityStrength[xIndex][yIndex] = math.Pow(  math.Pow(gravity.entries[xIndex][yIndex][0], 2.0) +  math.Pow(gravity.entries[xIndex][yIndex][1], 2.0), 0.5 )
		}
	}
	
	
	
	
	// set particle distribution function uniformly accross the fluid region
	for xIndex:=0;xIndex<D2Q9.nx;xIndex++{
		for yIndex:=0;yIndex<D2Q9.ny;yIndex++{
			if (!circleBoundaryR50(xIndex, yIndex)){
				for dIndex, _ := range D2Q9.directions{
					f.entries[xIndex][yIndex][dIndex] = D2Q9.weights[dIndex] * D2Q9.rho
				}				
			}
		}
	}
	
	// some arrays for keeping track of horizontal and vertical velocities 
	ux := D2Constant(nx, ny, 0.0)
	uy := D2Constant(nx, ny, 0.0)

	// arrays of constant random numbers for use in perterbations
	randomNumbers1 := D2Constant(nx, ny, 0.0)
	randomNumbers2 := D2Constant(nx, ny, 0.0)
	for index:=0;index<D2Q9.nx;index++{
		for innerIndex :=0;innerIndex<D2Q9.ny;innerIndex++{
			randomNumbers1[index][innerIndex] = (2.0*rand.Float64() - 1.0)
			randomNumbers2[index][innerIndex] = (2.0*rand.Float64() - 1.0)
		}
	}

	fmt.Println("starting computation") 
	displayTimes := false // display computation times for each step

	innerRadius := float64(40)
	innerCounter := 0.0
	outerCounter := 0.0
	for xIndex:=0;xIndex<D2Q9.nx;xIndex++{
		for yIndex:=0;yIndex<D2Q9.ny;yIndex++{
			if (!circleBoundaryR50(xIndex, yIndex)){
				if (math.Pow(float64(xIndex - 51), 2.0) + math.Pow(  float64(yIndex - 51), 2.0 ) < math.Pow(float64(innerRadius) + 3*randomNumbers1[xIndex][yIndex], 2.0)) {
					innerCounter += 1.0
				}else {
					outerCounter += 1.0
				}
			}
		}
	}

	// main loop
	for n:=0;n<nt;n++{
		fmt.Printf("\rstarting loop" + strconv.Itoa(n))
		// set energy distribution field according to heat source
		for xIndex:=0;xIndex<D2Q9.nx;xIndex++{
			for yIndex:=0;yIndex<D2Q9.ny;yIndex++{
				for dIndex, _ := range D2Q9.directions{

					if (!circleBoundaryR50(xIndex, yIndex)){		
						if (  math.Pow(float64(xIndex - 51), 2.0) + math.Pow(  float64(yIndex - 51), 2.0 ) < math.Pow(float64(innerRadius) + 3*randomNumbers1[xIndex][yIndex], 2.0) )   {
							g.entries[xIndex][yIndex][dIndex] += rho.entries[xIndex][yIndex] * D2Q9.Heat * D2Q9.weights[dIndex] * (1 + 0.1*randomNumbers1[xIndex][yIndex])
						} else {
							g.entries[xIndex][yIndex][dIndex] += -rho.entries[xIndex][yIndex] * D2Q9.Heat * D2Q9.weights[dIndex] * (innerCounter/outerCounter) * (1 + 0.1*randomNumbers1[xIndex][yIndex])
						}
						
						
					}
				}

			}	
		}


		
		start := time.Now() // ignore lines that have to do with time, they are used to test the performance of the code

		// perform the streaming step for both fields
		fstream = streamingStep2D(f, D2Q9, circleBoundaryR50)
		gstream = streamingStep2D(g, D2Q9, circleBoundaryR50) 
		if (displayTimes){
			fmt.Println("streaming took", time.Since(start))
		}

		f = fstream
		g = gstream
		
		
		
		start = time.Now()

		// get the macroscopic properties
		rho, ep, u = getMacro2DWithEp(f, g, D2Q9, circleBoundaryR50)

		if (displayTimes){
			fmt.Println("getting macros took", time.Since(start))
		}
		

		// every 100 steps we record the macroscopic properties
		if ( n % 100  == 0){
			name := "rho//rho"+strconv.Itoa(n)+".csv"
			writeArrayCSV(rho.entries, name)
			name = "ep//ep"+strconv.Itoa(n)+".csv"
			writeArrayCSV(ep.entries, name)
			for xIndex:=0;xIndex<D2Q9.nx;xIndex++{
				for yIndex:=0;yIndex<D2Q9.ny;yIndex++{
					ux[xIndex][yIndex] = u.entries[xIndex][yIndex][0]
					uy[xIndex][yIndex] = u.entries[xIndex][yIndex][1]
				}
			}
			name = "u//ux"+strconv.Itoa(n)+".csv"
			writeArrayCSV(gravityx, name)
			name = "u//uy"+strconv.Itoa(n)+".csv"
			writeArrayCSV(gravityy, name)
		}
		
		// get equilibrium quantities
		start = time.Now()
		geq = getEq2Dg(rho, ep, u, D2Q9)
		feq = getEq2Df(rho, u, D2Q9)
		if (displayTimes){
			fmt.Println("getting eq took", time.Since(start))
		}
		

		// perform the collision step
		start = time.Now()
		f = collisionStep(f, feq, D2Q9)
		g = collisionStep(g, geq, D2Q9)

		if (displayTimes){
			fmt.Println("collisin took", time.Since(start))
		}
		
		// Add the gravity term 
		start = time.Now()
		f = gravityAffect(f, gravity, rho, ep, D2Q9)

		if (displayTimes){
			fmt.Println("gravity took", time.Since(start))
		}	

		// the loop repeats, we do this many times
	}
	fmt.Println("Program terminated, done all timesteps requested. Try increasing nt to run for longer")
}
