package newton

import (
	"fmt"
	"math"
	"time"

	"github.com/almostinf/computational_mathematics/task2/matrix"
)

var Eps float64 = 1e-14

type NewtonMethod struct {
	F    []float64
	J    *matrix.Matrix
	X    []float64
	size int
}

func New() *NewtonMethod {
	n := &NewtonMethod{
		size: 10,
	}

	n.InitX()
	n.InitF()
	n.InitJ()

	return n
}

func (n *NewtonMethod) InitF() {
	n.F = make([]float64, 10)

	n.F[0] = math.Cos(n.X[0]*n.X[1]) - math.Exp(-3.0*n.X[2]) + n.X[3]*n.X[4]*n.X[4] - n.X[5] - math.Sinh(2.0*n.X[7])*n.X[8] + 2.0*n.X[9] + 2.0004339741653854440
	n.F[1] = math.Sin(n.X[0]*n.X[1]) + n.X[2]*n.X[8]*n.X[6] - math.Exp(-n.X[9]+n.X[5]) + 3.0*n.X[4]*n.X[4] - n.X[5]*(n.X[7]+1.0) + 10.886272036407019994
	n.F[2] = n.X[0] - n.X[1] + n.X[2] - n.X[3] + n.X[4] - n.X[5] + n.X[6] - n.X[7] + n.X[8] - n.X[9] - 3.1361904761904761904
	n.F[3] = 2.0*math.Cos(-n.X[8]+n.X[3]) + n.X[4]/(n.X[2]+n.X[0]) - math.Sin(n.X[1]*n.X[1]) + math.Cos(n.X[6]*n.X[9])*math.Cos(n.X[6]*n.X[9]) - n.X[7] - 0.1707472705022304757
	n.F[4] = math.Sin(n.X[4]) + 2.0*n.X[7]*(n.X[2]+n.X[0]) - math.Exp(-n.X[6]*(-n.X[9]+n.X[5])) + 2.0*math.Cos(n.X[1]) - 1.00/(n.X[3]-n.X[8]) - 0.3685896273101277862
	n.F[5] = math.Exp(n.X[0]-n.X[3]-n.X[8]) + n.X[4]*n.X[4]/n.X[7] + math.Cos(3.0*n.X[9]*n.X[1])/2.0 - n.X[5]*n.X[2] + 2.0491086016771875115
	n.F[6] = n.X[1]*n.X[1]*n.X[1]*n.X[6] - math.Sin(n.X[9]/n.X[4]+n.X[7]) + (n.X[0]-n.X[5])*math.Cos(n.X[3]) + n.X[2] - 0.7380430076202798014
	n.F[7] = n.X[4]*(n.X[0]-2.0*n.X[5])*(n.X[0]-2.0*n.X[5]) - 2.0*math.Sin(-n.X[8]+n.X[2]) + 1.5*n.X[3] - math.Exp(n.X[1]*n.X[6]+n.X[9]) + 3.5668321989693809040
	n.F[8] = 7.0/n.X[5] + math.Exp(n.X[4]+n.X[3]) - 2.0*n.X[1]*n.X[7]*n.X[9]*n.X[6] + 3.0*n.X[8] - 3.0*n.X[0] - 8.4394734508383257499
	n.F[9] = n.X[9]*n.X[0] + n.X[8]*n.X[1] - n.X[7]*n.X[2] + math.Sin(n.X[3]+n.X[4]+n.X[5])*n.X[6] - 0.78238095238095238096
}

func (n *NewtonMethod) InitJ() {
	n.J = matrix.New(10, 10)

	n.J.A[0][0] = -math.Sin(n.X[0]*n.X[1]) * n.X[1]
	n.J.A[0][1] = -math.Sin(n.X[0]*n.X[1]) * n.X[0]
	n.J.A[0][2] = 3 * math.Exp(-(3.0 * n.X[2]))
	n.J.A[0][3] = n.X[4] * n.X[4]
	n.J.A[0][4] = 2 * n.X[3] * n.X[4]
	n.J.A[0][5] = -1
	n.J.A[0][6] = 0
	n.J.A[0][7] = -2 * math.Cosh((2 * n.X[7])) * n.X[8]
	n.J.A[0][8] = -math.Sinh((2 * n.X[7]))
	n.J.A[0][9] = 2
	n.J.A[1][0] = math.Cos(n.X[0]*n.X[1]) * n.X[1]
	n.J.A[1][1] = math.Cos(n.X[0]*n.X[1]) * n.X[0]
	n.J.A[1][2] = n.X[8] * n.X[6]
	n.J.A[1][3] = 0
	n.J.A[1][4] = 6 * n.X[4]
	n.J.A[1][5] = -math.Exp(-n.X[9]+n.X[5]) - n.X[7] - 0.1e1
	n.J.A[1][6] = n.X[2] * n.X[8]
	n.J.A[1][7] = -n.X[5]
	n.J.A[1][8] = n.X[2] * n.X[6]
	n.J.A[1][9] = math.Exp(-n.X[9] + n.X[5])
	n.J.A[2][0] = 1
	n.J.A[2][1] = -1
	n.J.A[2][2] = 1
	n.J.A[2][3] = -1
	n.J.A[2][4] = 1
	n.J.A[2][5] = -1
	n.J.A[2][6] = 1
	n.J.A[2][7] = -1
	n.J.A[2][8] = 1
	n.J.A[2][9] = -1
	n.J.A[3][0] = -n.X[4] * math.Pow(n.X[2]+n.X[0], -2)
	n.J.A[3][1] = -2 * math.Cos(n.X[1]*n.X[1]) * n.X[1]
	n.J.A[3][2] = -n.X[4] * math.Pow(n.X[2]+n.X[0], -2)
	n.J.A[3][3] = -2 * math.Sin(-n.X[8]+n.X[3])
	n.J.A[3][4] = 1 / (n.X[2] + n.X[0])
	n.J.A[3][5] = 0
	n.J.A[3][6] = -2 * math.Cos(n.X[6]*n.X[9]) * math.Sin(n.X[6]*n.X[9]) * n.X[9]
	n.J.A[3][7] = -1
	n.J.A[3][8] = 2 * math.Sin(-n.X[8]+n.X[3])
	n.J.A[3][9] = -2 * math.Cos(n.X[6]*n.X[9]) * math.Sin(n.X[6]*n.X[9]) * n.X[6]
	n.J.A[4][0] = 2 * n.X[7]
	n.J.A[4][1] = -2 * math.Sin(n.X[1])
	n.J.A[4][2] = 2 * n.X[7]
	n.J.A[4][3] = math.Pow(-n.X[8]+n.X[3], -2)
	n.J.A[4][4] = math.Cos(n.X[4])
	n.J.A[4][5] = n.X[6] * math.Exp(-n.X[6]*(-n.X[9]+n.X[5]))
	n.J.A[4][6] = -(n.X[9] - n.X[5]) * math.Exp(-n.X[6]*(-n.X[9]+n.X[5]))
	n.J.A[4][7] = (2 * n.X[2]) + 2*n.X[0]
	n.J.A[4][8] = -math.Pow(-n.X[8]+n.X[3], -2)
	n.J.A[4][9] = -n.X[6] * math.Exp(-n.X[6]*(-n.X[9]+n.X[5]))
	n.J.A[5][0] = math.Exp(n.X[0] - n.X[3] - n.X[8])
	n.J.A[5][1] = -3.0 / 2.0 * math.Sin(3*n.X[9]*n.X[1]) * n.X[9]
	n.J.A[5][2] = -n.X[5]
	n.J.A[5][3] = -math.Exp(n.X[0] - n.X[3] - n.X[8])
	n.J.A[5][4] = 2 * n.X[4] / n.X[7]
	n.J.A[5][5] = -n.X[2]
	n.J.A[5][6] = 0
	n.J.A[5][7] = -n.X[4] * n.X[4] * math.Pow(n.X[7], (-2))
	n.J.A[5][8] = -math.Exp(n.X[0] - n.X[3] - n.X[8])
	n.J.A[5][9] = -3.0 / 2.0 * math.Sin(3*n.X[9]*n.X[1]) * n.X[1]
	n.J.A[6][0] = math.Cos(n.X[3])
	n.J.A[6][1] = 3 * n.X[1] * n.X[1] * n.X[6]
	n.J.A[6][2] = 1
	n.J.A[6][3] = -(n.X[0] - n.X[5]) * math.Sin(n.X[3])
	n.J.A[6][4] = math.Cos(n.X[9]/n.X[4]+n.X[7]) * n.X[9] * math.Pow(n.X[4], (-2))
	n.J.A[6][5] = -math.Cos(n.X[3])
	n.J.A[6][6] = math.Pow(n.X[1], 3)
	n.J.A[6][7] = -math.Cos(n.X[9]/n.X[4] + n.X[7])
	n.J.A[6][8] = 0
	n.J.A[6][9] = -math.Cos(n.X[9]/n.X[4]+n.X[7]) / n.X[4]
	n.J.A[7][0] = 2 * n.X[4] * (n.X[0] - 2*n.X[5])
	n.J.A[7][1] = -n.X[6] * math.Exp(n.X[1]*n.X[6]+n.X[9])
	n.J.A[7][2] = -2 * math.Cos(-n.X[8]+n.X[2])
	n.J.A[7][3] = 0.15e1
	n.J.A[7][4] = math.Pow(n.X[0]-2*n.X[5], 2)
	n.J.A[7][5] = -4 * n.X[4] * (n.X[0] - 2*n.X[5])
	n.J.A[7][6] = -n.X[1] * math.Exp(n.X[1]*n.X[6]+n.X[9])
	n.J.A[7][7] = 0
	n.J.A[7][8] = 2 * math.Cos(-n.X[8]+n.X[2])
	n.J.A[7][9] = -math.Exp(n.X[1]*n.X[6] + n.X[9])
	n.J.A[8][0] = -3
	n.J.A[8][1] = -2 * n.X[7] * n.X[9] * n.X[6]
	n.J.A[8][2] = 0
	n.J.A[8][3] = math.Exp((n.X[4] + n.X[3]))
	n.J.A[8][4] = math.Exp((n.X[4] + n.X[3]))
	n.J.A[8][5] = -0.7e1 * math.Pow(n.X[5], -2)
	n.J.A[8][6] = -2 * n.X[1] * n.X[7] * n.X[9]
	n.J.A[8][7] = -2 * n.X[1] * n.X[9] * n.X[6]
	n.J.A[8][8] = 3
	n.J.A[8][9] = -2 * n.X[1] * n.X[7] * n.X[6]
	n.J.A[9][0] = n.X[9]
	n.J.A[9][1] = n.X[8]
	n.J.A[9][2] = -n.X[7]
	n.J.A[9][3] = math.Cos(n.X[3]+n.X[4]+n.X[5]) * n.X[6]
	n.J.A[9][4] = math.Cos(n.X[3]+n.X[4]+n.X[5]) * n.X[6]
	n.J.A[9][5] = math.Cos(n.X[3]+n.X[4]+n.X[5]) * n.X[6]
	n.J.A[9][6] = math.Sin(n.X[3] + n.X[4] + n.X[5])
	n.J.A[9][7] = -n.X[2]
	n.J.A[9][8] = n.X[1]
	n.J.A[9][9] = n.X[0]
}

func (n *NewtonMethod) InitX() {
	n.X = make([]float64, 10)

	n.X[0] = 0.5
	n.X[1] = 0.5
	n.X[2] = 1.5
	n.X[3] = -1
	n.X[4] = -0.5
	n.X[5] = 1.5 
	n.X[6] = 0.5
	n.X[7] = -0.5
	n.X[8] = 1.5
	n.X[9] = -1.5
}

// https://slemeshevsky.github.io/num-mmf/snes/html/._snes-FlatUI001.html
func (n *NewtonMethod) SolveSystem() {
	fmt.Println("--------------------------------------------------")
	fmt.Println("NEWTON METHOD")
	operations := 0
	iterations := 0
	now := time.Now()
	n.InitX()
	for {
		n.InitF()
		n.InitJ()

		err := n.J.LUDecomposition()
		if err != nil {
			fmt.Println(err)
			return
		}

		fixedF, err := matrix.MultOnVecRight(n.J.P, n.F)
		if err != nil {
			fmt.Println(err)
			return
		}
		operations += n.size

		dx, err := n.J.SLAESolution(fixedF)
		operations += n.J.Operations

		if err != nil {
			fmt.Println("Error: matrix and vector dimensions don't match")
			return
		}

		for i := range n.X {
			n.X[i] -= dx[i]
			operations++
		}

		if matrix.MaxInVec(dx) < Eps {
			break
		}
		iterations++
	}

	dur := time.Since(now).Seconds()
	fmt.Println("Solution: ", n.X)
	fmt.Println("Number of operations: ", operations)
	fmt.Println("Number of iterations: ", iterations)
	fmt.Println("Elapsed time in second: ", dur)
	fmt.Println("--------------------------------------------------")
	fmt.Println()
}

func (n *NewtonMethod) ModifiedSolveSystem() {
	fmt.Println("--------------------------------------------------")
	fmt.Println("MODIFIED NEWTON METHOD")
	operations := 0
	iterations := 0
	now := time.Now()

	n.InitX()
	n.InitJ()
	err := n.J.LUDecomposition()
	if err != nil {
		fmt.Println(err)
		return
	}

	for iterations < 1000 {
		n.InitF()

		fixedF, err := matrix.MultOnVecRight(n.J.P, n.F)
		if err != nil {
			fmt.Println(err)
			return
		}
		operations += n.size * n.size

		dx, err := n.J.SLAESolution(fixedF)

		if err != nil {
			fmt.Println("Error: matrix and vector dimensions don't match")
			return
		}

		for i := range n.X {
			n.X[i] -= dx[i]
			operations++
		}

		if matrix.MaxInVec(dx) > Eps {
			iterations++
			continue
		}
		break
	}

	operations += n.J.Operations
	dur := time.Since(now).Seconds()
	fmt.Println("Solution: ", n.X)
	fmt.Println("Number of operations: ", operations)
	fmt.Println("Number of iterations: ", iterations)
	fmt.Println("Elapsed time in second: ", dur)
	fmt.Println("--------------------------------------------------")
	fmt.Println()
}

func (n *NewtonMethod) ModifiedSolveSystemOnlyKIterations(k int) {
	fmt.Println("--------------------------------------------------")
	fmt.Println("MODIFIED METHOD OF NEWTON AFTER K ITERATIONS")
	operations := 0
	iterations := 0
	now := time.Now()
	n.InitX()
	n.InitJ()

	err := n.J.LUDecomposition()
	if err != nil {
		fmt.Println(err)
		return
	}

	for iterations < 1000 {
		n.InitF()

		if k >= 0 {
			operations += n.J.Operations
			n.InitJ()

			err := n.J.LUDecomposition()
			if err != nil {
				fmt.Println(err)
				return
			}
		}

		fixedF, err := matrix.MultOnVecRight(n.J.P, n.F)
		if err != nil {
			fmt.Println(err)
			return
		}
		operations += n.size * n.size

		dx, err := n.J.SLAESolution(fixedF)

		if err != nil {
			fmt.Println("Error: matrix and vector dimensions don't match")
			return
		}

		for i := range n.X {
			n.X[i] -= dx[i]
			operations++
		}

		if matrix.MaxInVec(dx) > Eps {
			k--
			iterations++
			continue
		}
		break
	}

	operations += n.J.Operations
	dur := time.Since(now).Seconds()
	fmt.Println("Solution: ", n.X)
	fmt.Println("Number of operations: ", operations)
	fmt.Println("Number of iterations: ", iterations)
	fmt.Println("Elapsed time in second: ", dur)
	fmt.Println("--------------------------------------------------")
	fmt.Println()
}

func (n *NewtonMethod) SolveSystemWithMIterations(m int) {
	fmt.Println("--------------------------------------------------")
	fmt.Println("NEWTON METHOD WITH RECOUNTING J EVERY M ITERATIONS")
	operations := 0
	iterations := 0
	now := time.Now()
	n.InitX()
	n.InitJ()

	err := n.J.LUDecomposition()
	if err != nil {
		fmt.Println(err)
		return
	}

	for iterations < 1000 {
		n.InitF()

		if iterations%m == 0 {
			operations += n.J.Operations
			n.InitJ()

			err := n.J.LUDecomposition()
			if err != nil {
				fmt.Println(err)
				return
			}
		}

		fixedF, err := matrix.MultOnVecRight(n.J.P, n.F)
		if err != nil {
			fmt.Println(err)
			return
		}
		operations += n.size * n.size

		dx, err := n.J.SLAESolution(fixedF)

		if err != nil {
			fmt.Println("Error: matrix and vector dimensions don't match")
			return
		}

		for i := range n.X {
			n.X[i] -= dx[i]
			operations++
		}

		if matrix.MaxInVec(dx) > Eps {
			iterations++
			continue
		}
		break
	}

	operations += n.J.Operations
	dur := time.Since(now).Seconds()
	fmt.Println("Solution: ", n.X)
	fmt.Println("Number of operations: ", operations)
	fmt.Println("Number of iterations: ", iterations)
	fmt.Println("Elapsed time in second: ", dur)
	fmt.Println("--------------------------------------------------")
	fmt.Println()
}

func (n *NewtonMethod) MethodsTransition(m, k int) {
	fmt.Println("--------------------------------------------------")
	fmt.Println("TRANSITION FROM NEWTON METHOD TO MODIFIED NEWTON METHOD AFTER K ITERATIONS")
	operations := 0
	iterations := 0
	now := time.Now()
	n.InitX()
	n.InitJ()

	err := n.J.LUDecomposition()
	if err != nil {
		fmt.Println(err)
		return
	}

	for iterations < 1000 {
		n.InitF()

		if m >= 0 || iterations%k == 0 {
			operations += n.J.Operations
			n.InitJ()

			err := n.J.LUDecomposition()
			if err != nil {
				fmt.Println(err)
				return
			}
		}

		fixedF, err := matrix.MultOnVecRight(n.J.P, n.F)
		if err != nil {
			fmt.Println(err)
			return
		}
		operations += n.size * n.size

		dx, err := n.J.SLAESolution(fixedF)

		if err != nil {
			fmt.Println("Error: matrix and vector dimensions don't match")
			return
		}

		for i := range n.X {
			n.X[i] -= dx[i]
			operations++
		}

		if matrix.MaxInVec(dx) > Eps {
			m--
			iterations++
			continue
		}
		break
	}

	operations += n.J.Operations
	dur := time.Since(now).Seconds()
	fmt.Println("Solution: ", n.X)
	fmt.Println("Number of operations: ", operations)
	fmt.Println("Number of iterations: ", iterations)
	fmt.Println("Elapsed time in second: ", dur)
	fmt.Println("--------------------------------------------------")
	fmt.Println()
}
