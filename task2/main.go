package main

import (
	"fmt"

	"github.com/almostinf/computational_mathematics/task2/matrix"
)

const N int = 5

func main() {
	// m := matrix.New(N, N)
	// if err := m.LUDecomposition(); err != nil {
	// 	fmt.Println(err)
	// 	return
	// }
	// matrix.Print(m.A, "A")
	// matrix.Print(m.L, "L")
	// matrix.Print(m.U, "U")
	// matrix.Print(m.P, "P")
	// matrix.Print(m.Q, "Q")

	// res, err := m.Determinant()
	// if err != nil {
	// 	fmt.Println(err)
	// 	return
	// }

	// fmt.Println("Determinant: ", res)
	// fmt.Println()

	// inv, _ := m.Inverse()
	// matrix.Print(inv, "Inverse")

	// matrix.Print(matrix.Mult(inv, m.A), "Eye")

	// degen_m := matrix.New(N, N)

	// for i := 0; i < N; i++ {
	// 	degen_m.A[i][2] = degen_m.A[i][0] + degen_m.A[i][1]
	// 	degen_m.A[i][3] = degen_m.A[i][0] - degen_m.A[i][1]
	// }

	// if err := degen_m.LUDecomposition(); err != nil {
	// 	fmt.Println(err)
	// 	return
	// }

	// matrix.Print(degen_m.A, "A")
	// matrix.Print(degen_m.L, "L")
	// matrix.Print(degen_m.U, "U")
	// matrix.Print(degen_m.P, "P")
	// matrix.Print(degen_m.Q, "Q")

	// fmt.Println("Rank:", degen_m.GetRank())
	// fmt.Println()

	// degen_b := matrix.GetRandomVector(N)
	// m_b,_ := matrix.MultOnVecRight(degen_m.A, degen_b)
	// fmt.Println("degen_b: ", m_b)
	// fmt.Println()

	// fixed_b, _ := matrix.MultOnVecRight(degen_m.P, m_b)
	// fmt.Println("fixed_b: ", fixed_b)
	// fmt.Println()

	// fmt.Println("SLAE Solution of degen: ")
	// solution, err := degen_m.SLAESolution(fixed_b)
	// if err != nil {
	// 	fmt.Println(err)
	// 	return
	// }
	// fmt.Println(solution)
	// fmt.Println()

	// ax, _ := matrix.MultOnVecRight(degen_m.A, solution)
	// fmt.Println("ax: ")
	// fmt.Println(ax)
	// fmt.Println()

	// matrix.Print(degen_m.A, "degen_A")
	// matrix.Print(degen_m.L, "degen_L")
	// matrix.Print(degen_m.U, "degen_U")
	// matrix.Print(degen_m.P, "degen_P")
	// matrix.Print(degen_m.Q, "degen_Q")

	// --------------------------------------------------
	m := matrix.New(N, N)
	m.QR()
	matrix.Print(m.A, "A")
	matrix.Print(m.Q, "Q")
	matrix.Print(m.R, "R")

	qr := matrix.Mult(m.Q, m.R)
	matrix.Print(qr, "QR")

	b := matrix.GetRandomVector(N)
	fmt.Println("b: ", b)
	fmt.Println()

	x, _ := m.SolveQR(b)
	Ax, _ := matrix.MultOnVecRight(m.A, x)
	fmt.Println("Ax: ", Ax)

	m = matrix.GetMatrixDiagonalDominance(N)

	x, iterations, prioriEstimate := m.JacobiMethod(b)
	fmt.Println("Jacobi")
	fmt.Println("x: ", x)
	fmt.Println("iterations: ", iterations)
	fmt.Println("priori estimate: ", prioriEstimate)

	Ax, _ = matrix.MultOnVecRight(m.A, x)
	fmt.Println("Ax: ", Ax)

	fmt.Println()

	x, iterations, prioriEstimate = m.SeidelMethod(b)
	fmt.Println("Seidel")
	fmt.Println("x: ", x)
	fmt.Println("iterations: ", iterations)
	fmt.Println("priori estimate: ", prioriEstimate)

	Ax, _ = matrix.MultOnVecRight(m.A, x)
	fmt.Println("Ax: ", Ax)

	fmt.Println()
	fmt.Println("-------------------------")
	fmt.Println()

	pm := matrix.GetPositivelyDefiniteMatrix(N)
	matrix.Print(pm.A, "Positive matrix")
	fmt.Println()

	x, iterations, prioriEstimate = pm.JacobiMethod(b)
	fmt.Println("Jacobi")
	fmt.Println("x: ", x)
	fmt.Println("iterations: ", iterations)
	fmt.Println("priori estimate: ", prioriEstimate)

	Ax, _ = matrix.MultOnVecRight(pm.A, x)
	fmt.Println("Ax: ", Ax)

	fmt.Println()

	x, iterations, prioriEstimate = pm.SeidelMethod(b)
	fmt.Println("Seidel")
	fmt.Println("x: ", x)
	fmt.Println("iterations: ", iterations)
	fmt.Println("priori estimate: ", prioriEstimate)

	Ax, _ = matrix.MultOnVecRight(pm.A, x)
	fmt.Println("Ax: ", Ax)
}
