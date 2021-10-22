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

type boundary func(cord1 int, cord2 int) bool 

func streamingStep2D(f1 D2Field, l lattice, b boundary) D2Field{
	var f2 D2Field
	entries := D3Constant(l.nx, l.ny, l.Q, 0.0)
	f2.entries = entries 
	f2.name = f1.name
	//copy(f2.entries, f1.entries)




	for dIndex, d := range l.directions{
	for i:=0;i<l.nx;i++{		
		for j:=0;j<l.ny;j++{
				if (b(i-d[0],j-d[1])){
					f2.entries[i][j][dIndex] = f1.entries[i][j][l.opposites[dIndex]]
				} else{
					f2.entries[i][j][dIndex] = f1.entries[i-d[0]][j-d[1]][dIndex]
				}
			}
		}
	}
	return f2
}

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
	var rho D2Property
	var ep D2Property
	var u D2Field
	
	rhoentries := D2Constant(l.nx, l.ny, 0.0)
	epentries := D2Constant(l.nx, l.ny, 0.0)
	uentries := D3Constant(l.nx, l.ny, l.D, 0.0)
	
	rho.entries = rhoentries
	ep.entries = epentries
	u.entries = uentries
	
	for xIndex:=0;xIndex<l.nx;xIndex++{
		for yIndex:=0;yIndex<l.ny;yIndex++{

			if (!b(xIndex, yIndex)){

				for dIndex, _ := range l.directions{
					rho.entries[xIndex][yIndex] += f.entries[xIndex][yIndex][dIndex]
				}

				for dIndex, d := range l.directions{
					u.entries[xIndex][yIndex][0] += float64(d[0]) * f.entries[xIndex][yIndex][dIndex]/rho.entries[xIndex][yIndex]
					u.entries[xIndex][yIndex][1] += float64(d[1]) * f.entries[xIndex][yIndex][dIndex]/rho.entries[xIndex][yIndex]
					ep.entries[xIndex][yIndex] += g.entries[xIndex][yIndex][dIndex]/rho.entries[xIndex][yIndex]
				}
			}
		}
	}
	return rho, ep, u // rho and u are zero in the dead nodes
}

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
				
				f.entries[xIndex][yIndex][dIndex] = rho.entries[xIndex][yIndex] * l.weights[dIndex] * (term1 + term2 + term3 + term4)
			}
		}
	}
	return f
}

func getEq2Dg(rho D2Property, ep D2Property, u D2Field, l lattice) D2Field{
	// this looks good
	
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
	var rtn D2Field
	entries := D3Constant(l.nx, l.ny, l.Q, 0.0)
	rtn.entries = entries
	rtn.name = f.name
	
	var tau float64
	if ( f.name == "f"){
		tau = l.tau_f
	} else if (f.name == "g"){
		tau = l.tau_g
	} else {
		panic("unknown field being used")
	}
	
	for dIndex, _ := range l.directions{
		for xIndex:=0;xIndex<l.nx;xIndex++{
			for yIndex:=0;yIndex<l.ny;yIndex++{
				rtn.entries[xIndex][yIndex][dIndex] = f.entries[xIndex][yIndex][dIndex] +  1.0/(tau) *  ( feq.entries[xIndex][yIndex][dIndex] - f.entries[xIndex][yIndex][dIndex] )
			}
		}
	}
	return rtn
}


func gravityAffect(f D2Field, grav D2Field, rho D2Property, ep D2Property, l lattice) D2Field{
	// pingback
	
	var rtn D2Field
	entries := D3Constant(l.nx, l.ny, l.Q, 0.0)
	rtn.entries = entries
	rtn.name = "f"
	
	for dIndex, d := range l.directions{
		for xIndex:=0;xIndex<l.nx;xIndex++{
			for yIndex:=0;yIndex<l.ny;yIndex++{



				//rtn.entries[xIndex][yIndex][dIndex] += -1.0/(l.tau_grav) * 1.0 * l.weights[dIndex] * l.rho * l.alpha * 1.0 * float64(d[1]) * ( ep.entries[xIndex][yIndex]*2.0/float64(l.D) )

				

				if (dIndex == 0){
					rtn.entries[xIndex][yIndex][dIndex] = f.entries[xIndex][yIndex][dIndex] + 0.0
				} else{

					//term := -1/(l.tau_grav) * 3.0 * l.weights[dIndex] * float64(d[1]) * (-rho.entries[xIndex][yIndex]) * l.alpha * 1.0 * ( ep.entries[xIndex][yIndex]*2.0/float64(l.D) )
					//term = 0.0
					term := -1.0/(l.tau_grav) * l.weights[dIndex] * rho.entries[xIndex][yIndex] * ( ep.entries[xIndex][yIndex]*2.0/float64(l.D) ) * l.alpha * math.Pow( math.Pow(float64(d[0]), 2.0) + math.Pow(float64(d[1]), 2.0), -0.5) * (grav.entries[xIndex][yIndex][0]*float64(d[0]) + grav.entries[xIndex][yIndex][1]*float64(d[1]))
					//fmt.Println(term)
					rtn.entries[xIndex][yIndex][dIndex] = f.entries[xIndex][yIndex][dIndex] + term
				}
				
				/*
				if (dIndex == 3 || dIndex == l.opposites[3] ){
					term := -1/(l.tau_grav) * 3.0 * l.weights[dIndex] * float64(d[1]) * (-rho.entries[xIndex][yIndex]) * l.alpha * 1.0 * ( ep.entries[xIndex][yIndex]*2.0/float64(l.D) )
					rtn.entries[xIndex][yIndex][dIndex] += term
				}
				*/
				

				
			}
		}
	}
	return rtn
}


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


func circleBoundaryR50(xIndex int, yIndex int) bool{	
	if (  math.Pow( float64(xIndex - 51), 2.0) + math.Pow(float64(yIndex - 51), 2.0) <= math.Pow(float64(50), 2.0) ){
		return false
	}
	return true
	
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


func main(){
	s1 := rand.NewSource(time.Now().UnixNano())
	rand.New(s1)
	fmt.Println("Running Main")
	// lets define the D2Q9 lattice
	
	var dt, dx float64 = 1, 1
	var S float64 = dx/dt
	
	Q := 9
	D := 2
	
	
	var D2Q9 lattice

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
	
	cs := [][]int{v1, v2, v3, v4, v5, v6, v7, v8, v9}
	// here I write the opposites
	
	op := []int{0, 2, 1, 4, 3, 6, 5, 8, 7}
	

	vels := []float64{ 1.0, 3.0/(math.Pow(S, 2.0)), 9.0/(2.0*math.Pow(S,4.0)), -3.0/(2.0*math.Pow(S, 2.0))}
	fmt.Println(vels)

	// trying a different eq method
	//vels = []float64{ 1.0, 0.0, 0.0, 0.0}
	
	// define the discritization
	nt := 10000000
	nx := 103
	ny := 103
		
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
	D2Q9.Cv = 400
	D2Q9.Heat = 0.01



	
	
	
	// lets define our constants 
	var vf float64 = 0.1 // kinematic viscosity
	var tau_f float64 = vf/(math.Pow(S/math.Sqrt(3), 2.0) * dt) + 0.5 // check this
	
	var vg float64 = 0.00002 // a measure of thermal diffusivity
	//var vg float64 = 10000 // a measure of thermal diffusivity
	var tau_g float64 = vg/(math.Pow(S/math.Sqrt(3), 2.0) * dt  ) + 0.5 // relaxation time associated with thermal diffusion
	//tau_g = 0.0 // this should throw a million errors // so clearly this isnt being used properly by the program.

	var vgrav float64 = 0.6
	var tau_grav float64 = vgrav/(math.Pow(S/math.Sqrt(3), 2.0) * dt  ) + 0.5
	
	
	var rhoBack float64 = 1000.0
	var alpha float64 = 0.1
	
	// set some fluid properties
	D2Q9.rho = rhoBack
	D2Q9.alpha = alpha
	
	
	
	
	// set the relaxation times into our lattice
	D2Q9.tau_f = tau_f
	D2Q9.tau_g = tau_g
	D2Q9.tau_grav = tau_grav

	maxG := -1.0
	Ra := maxG*D2Q9.alpha*D2Q9.Heat*math.Pow(50.0, 5.0)/(D2Q9.rho * vf * math.Pow(vg, 2.0)*D2Q9.Cv )
	Pr := vf/vg
	fmt.Println("Ra: ", Ra)
	fmt.Println("Pr:", Pr)
	

	
	// then we define our arrays
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
	
	// populate gravity Here I will assume a very simple gravitational field which we can change later
	gravity.name = "gravity"
	gravity.entries = D3Constant(nx, ny, 2, 0.0)
	gravityx := D2Constant(nx, ny, 0)
	gravityy := D2Constant(nx, ny, 0)
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
		}
	}
	
	
	
	// lets set a simple flow in the x direction
	fmt.Println("starting initial conditions")
	for xIndex:=0;xIndex<D2Q9.nx;xIndex++{
		for yIndex:=0;yIndex<D2Q9.ny;yIndex++{
			if (!circleBoundaryR50(xIndex, yIndex)){
				for dIndex, _ := range D2Q9.directions{
					f.entries[xIndex][yIndex][dIndex] = D2Q9.weights[dIndex] * D2Q9.rho
				}				
			}
		}
	}
	
	// We will also need to set the conditions for the q field initially.
	ux := D2Constant(nx, ny, 0.0)
	uy := D2Constant(nx, ny, 0.0)

	randomNumbers1 := D2Constant(nx, ny, 0.0)
	randomNumbers2 := D2Constant(nx, ny, 0.0)

	for index:=0;index<D2Q9.nx;index++{
		for innerIndex :=0;innerIndex<D2Q9.ny;innerIndex++{
			randomNumbers1[index][innerIndex] = (2.0*rand.Float64() - 1.0)
			randomNumbers2[index][innerIndex] = (2.0*rand.Float64() - 1.0)
		}
	}




	
	fmt.Println("starting computation")
	displayTimes := false
	for n:=0;n<nt;n++{

		
		fmt.Printf("\rstarting loop" + strconv.Itoa(n))

		// fix the ep value (fix the temperature)
		// fix g at the bonudry aswell
		for xIndex:=0;xIndex<D2Q9.nx;xIndex++{
			for yIndex:=0;yIndex<D2Q9.ny;yIndex++{


				if ( math.Pow(float64(xIndex - 51),2.0) + math.Pow(float64(yIndex-51),2.0) <= math.Pow(25,2.0) ){
					for dIndex, _ := range D2Q9.directions{
						g.entries[xIndex][yIndex][dIndex] += rho.entries[xIndex][yIndex] * (2.0/float64(D2Q9.D) * D2Q9.Heat/D2Q9.Cv)*(1+0.8*randomNumbers1[xIndex][yIndex]) * D2Q9.weights[dIndex]
					}
				} else{
					for dIndex, _ := range D2Q9.directions{
						g.entries[xIndex][yIndex][dIndex] -= rho.entries[xIndex][yIndex] * 1.0/3.0 * (2.0/float64(D2Q9.D) * D2Q9.Heat/D2Q9.Cv)*(1+0.8*randomNumbers2[xIndex][yIndex])  * D2Q9.weights[dIndex]
					}

				}


				/*
				if (yIndex<3){
					//ep.entries[xIndex][yIndex] = -0.001 + 0.00025*(math.Sin( 1.0/10.0 * float64(xIndex)))

					// fix g here
					for dIndex, _ := range D2Q9.directions{
						g.entries[xIndex][yIndex][dIndex] = rho.entries[xIndex][yIndex] * (2.0/float64(D2Q9.D) * (-0.001 + 0.00025*randomNumbers1[xIndex]  )) * D2Q9.weights[dIndex]
					}


				}
				if (yIndex > 96){
					//ep.entries[xIndex][yIndex] = +0.001 + 0.00025*(math.Sin(1.0/10.0 * float64(xIndex)))

					// fix g here
					for dIndex, _ := range D2Q9.directions{
						g.entries[xIndex][yIndex][dIndex] = rho.entries[xIndex][yIndex] * (2.0/float64(D2Q9.D) * (0.001 + 0.00025*randomNumbers2[xIndex]  )) * D2Q9.weights[dIndex]
					}

				}
				*/

			}	
		}


		start := time.Now()
		fstream = streamingStep2D(f, D2Q9, circleBoundaryR50)
		gstream = streamingStep2D(g, D2Q9, circleBoundaryR50) // there is a problem here, somehow g gets to be (after this) such that ep grows very big/wrong
		if (displayTimes){
			fmt.Println("streaming took", time.Since(start))
		}
		f = fstream
		g = gstream
		
		
		
		
		//rho, u = getMacro2D(f, D2Q9, circleBoundaryThousand) // this is used for the code without thermally driven flows
		start = time.Now()
		rho, ep, u = getMacro2DWithEp(f, g, D2Q9, circleBoundaryR50)

		if (displayTimes){
			fmt.Println("getting macros took", time.Since(start))
		}
		


		if ( n % 100  == 0){
			//fmt.Println("writting rho to csv")
			name := "rho//rho"+strconv.Itoa(n)+".csv"
			writeArrayCSV(rho.entries, name)
			name = "ep//ep"+strconv.Itoa(n)+".csv"
			writeArrayCSV(ep.entries, name)
			//fmt.Println("done writting rho")
			
			
			for xIndex:=0;xIndex<D2Q9.nx;xIndex++{
				for yIndex:=0;yIndex<D2Q9.ny;yIndex++{
					ux[xIndex][yIndex] = u.entries[xIndex][yIndex][0]
					uy[xIndex][yIndex] = u.entries[xIndex][yIndex][1]
				}
			}
			//fmt.Println("writting u to csv")
			name = "u//ux"+strconv.Itoa(n)+".csv"
			writeArrayCSV(ux, name)
			name = "u//uy"+strconv.Itoa(n)+".csv"
			writeArrayCSV(uy, name)
			//fmt.Println("done writting u")
		}




	
		
		//feq = getEq2D(rho, u, D2Q9) // this is used for the code without thermally driven flows
		
		start = time.Now()
		geq = getEq2Dg(rho, ep, u, D2Q9)
		feq = getEq2Df(rho, u, D2Q9)
		if (displayTimes){
			fmt.Println("getting eq took", time.Since(start))
		}
		


		start = time.Now()
		f = collisionStep(f, feq, D2Q9)
		g = collisionStep(g, geq, D2Q9)

		if (displayTimes){
			fmt.Println("collisin took", time.Since(start))
		}
		

		start = time.Now()
		f = gravityAffect(f, gravity, rho, ep, D2Q9)

		if (displayTimes){
			fmt.Println("gravity took", time.Since(start))
		}	
		
		// theres another term isnt there?

		
		// record the results



	
	}
	
	
	
	
	
	

		
	
	
	
	
}
