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
	entries [][]float64
	name string
}


type lattice struct{
	weights []float64
	directions [][]int
	opposites []int
	vels []float64
	eqWeights []float64
	nx int
	ny int
	Q int
	D int
	dt float64
	dx float64
	S float64
	tau_f float64
	tau_g float64
	
}

type boundary func(cord1 int, cord2 int) bool

func streamingStep2D(f1 D2Field, l lattice, b boundary) D2Field{
// somehow this step induces asymetry
	var f2 D2Field
	entries := D3Constant(l.nx, l.ny, l.Q, 0.0)
	f2.entries = entries


	// okay so we take a field and we stream it. Thats all that we need to do.
	for dIndex, d := range l.directions{
		for i:=0;i<l.nx;i++{		
			for j:=0;j<l.ny;j++{
				if (b(i-d[0],j-d[1])){
					f2.entries[i][j][dIndex] = f1.entries[i][j][l.opposites[dIndex]]
				} else{
					f2.entries[i][j][dIndex] = f1.entries[i-d[0]][j-d[1]][dIndex]
				}
				
				
				/*
				// try to write this properly
				if (b(i+d[0], j+d[1])) {
					f2.entries[i][j][l.opposites[dIndex]] = f1.entries[i][j][dIndex]
					
				} else {
					f2.entries[i+d[0]][j+d[1]][dIndex] = f1.entries[i][j][dIndex]
				}
				
				*/
			}
		}
	}
	
	
	
	return f2
}

func getMacro2D(f D2Field, l lattice, b boundary) (D2Property, D2Field){
	var rho D2Property
	var u D2Field
	
	rhoentries := D2Constant(l.nx, l.ny, 0.0)
	uentries := D3Constant(l.nx, l.ny, l.D  , 0.0)

	for xIndex:=0;xIndex<l.nx;xIndex++{
		for yIndex:=0;yIndex<l.ny;yIndex++{
		
		
		
			if (!b(xIndex, yIndex)){
				for dIndex, _ := range l.directions{
					rhoentries[xIndex][yIndex] += f.entries[xIndex][yIndex][dIndex]
				}
				for dIndex, d := range l.directions{

					uentries[xIndex][yIndex][0] += float64(d[0]) * f.entries[xIndex][yIndex][dIndex]/rhoentries[xIndex][yIndex]
					uentries[xIndex][yIndex][1] += float64(d[1]) * f.entries[xIndex][yIndex][dIndex]/rhoentries[xIndex][yIndex]
				}
			}
		}
	}
	rho.entries = rhoentries
	u.entries = uentries
	return rho, u // rho and u are zero in the dead nodes
}


func getEq2D(rho D2Property, u D2Field, l lattice) D2Field{
	// this looks good
	
	var f D2Field
	entries := D3Constant(l.nx, l.ny, l.Q, 0.0)
	f.entries = entries
	
	
	
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


func collisionStep(f D2Field, feq D2Field, l lattice) D2Field{
	var rtn D2Field
	entries := D3Constant(l.nx, l.ny, l.Q, 0.0)
	rtn.entries = entries
	
	for dIndex, _ := range l.directions{
		for xIndex:=0;xIndex<l.nx;xIndex++{
			for yIndex:=0;yIndex<l.ny;yIndex++{
				rtn.entries[xIndex][yIndex][dIndex] = math.Abs(f.entries[xIndex][yIndex][dIndex] +  1.0/(l.tau_f) *  ( feq.entries[xIndex][yIndex][dIndex] - f.entries[xIndex][yIndex][dIndex] ))
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


func circleBoundaryThousand(xIndex int, yIndex int) bool{

	if (45-30<xIndex && 55+30>xIndex && 45-30<yIndex && 55+30>yIndex){
		return false
	}
	return true
	/*
	if (  math.Pow( float64(xIndex - 50), 2.0) + math.Pow(float64(yIndex - 50), 2.0) < math.Pow(float64(45), 2.0) ){
		return false
	}
	return true
	*/
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
	fmt.Println(ws[0])
	
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
	
	// define the discritization
	nt := 100
	nx := 100
	ny := 100
		
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
	
	
	
	// lets define our constants 
	var vf float64 = 0.00001 // kinematic viscosity
	var tau_f float64 = vf/(math.Pow(S/math.Sqrt(3), 2.0)*dt) + 0.5 // check this
	
	D2Q9.tau_f = tau_f
	

	
	// then we define our arrays
	var f D2Field
	var fstream D2Field
	var feq D2Field
	var rho D2Property
	var u D2Field
	
	f.entries = D3Constant(nx, ny, Q, 0.0)
	fstream.entries =  D3Constant(nx, ny, Q, 0.0)
	feq.entries = D3Constant(nx, ny, Q, 0.0)
	rho.entries = D2Constant(nx, ny, 1.0)
	u.entries = D3Constant(nx, ny, D, 0.0)
	
	
	// lets set a simple flow in the x direction
	fmt.Println("starting initial conditions")
	for xIndex:=0;xIndex<D2Q9.nx;xIndex++{
		for yIndex:=0;yIndex<D2Q9.ny;yIndex++{
			if (!circleBoundaryThousand(xIndex,yIndex)){
				for dIndex, _ := range D2Q9.directions{
				
					if (40 < xIndex && 60 > xIndex && 40 < yIndex && 60 > yIndex){
						f.entries[xIndex][yIndex][dIndex] = 0.001*rand.Float64()
					} else {
						f.entries[xIndex][yIndex][dIndex] = 0.001*rand.Float64()
					}
					
				}
				/*
				if (xIndex > 40 && xIndex < 60 && yIndex > 40 && yIndex < 60){
					f.entries[xIndex][yIndex][1] = 0.001*rand.Float64()
				} else{
					f.entries[xIndex][yIndex][dIndex] = 0.0 + 0.0001*rand.Float64()
				}
				*/
				
					
				
				
			}
		}
	}
	
	//feq := D3Zeros(nx, ny, Q)
	//Deltaf := D3Zeros(nx, ny, Q)
	
	ux := D2Constant(nx, ny, 0.0)
	uy := D2Constant(nx, ny, 0.0)
	
	for n:=0;n<nt;n++{
		// The only initial condition that we have is for f everywhere in the domainf
		
		fstream = streamingStep2D(f, D2Q9, circleBoundaryThousand)
		f = fstream
		rho, u = getMacro2D(f, D2Q9, circleBoundaryThousand) // looks good
		feq = getEq2D(rho, u, D2Q9)


	
		
		// record the results
		fmt.Println("writting rho to csv")
		name := "rho//rho"+strconv.Itoa(n)+".csv"
		writeArrayCSV(rho.entries, name)
		fmt.Println("done writting rho")
		
		
		for xIndex:=0;xIndex<D2Q9.nx;xIndex++{
			for yIndex:=0;yIndex<D2Q9.ny;yIndex++{
				ux[xIndex][yIndex] = u.entries[xIndex][yIndex][0]
				uy[xIndex][yIndex] = u.entries[xIndex][yIndex][1]
			}
		}
		fmt.Println("writting u to csv")
		name = "u//ux"+strconv.Itoa(n)+".csv"
		writeArrayCSV(ux, name)
		name = "u//uy"+strconv.Itoa(n)+".csv"
		writeArrayCSV(uy, name)
		fmt.Println("done writting u")
		
		f = collisionStep(f, feq, D2Q9)
		
		
		
		
		fmt.Println("step: " + strconv.Itoa(n))
	
	}
	
	
	
	
	
	

		
	
	
	
	
}
