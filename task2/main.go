package main

import (
	"fmt"

	"github.com/almostinf/computational_mathematics/task2/matrix"
)

const N int = 3

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
	// matrix.Print(matrix.Mult(m.L, m.U), "Mult")

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

	degen_m := matrix.New(N, N)
	matr := [][]float64{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}
	degen_m.Set(matr)

	if err := degen_m.LUDecomposition(); err != nil {
		fmt.Println(err)
		return
	}

	fmt.Println("Rank:", degen_m.GetRank())
	fmt.Println()

	degen_b := []float64{2.,10.,10.}
	fmt.Println("degen_b: ", degen_b)
	fmt.Println()

	fixed_b, _ := matrix.MultOnVecRight(degen_m.P, degen_b)
	fmt.Println("fixed_b: ", fixed_b)
	fmt.Println()

	fmt.Println("SLAE Solution of degen: ")
	solution, err := degen_m.SLAESolution(fixed_b)
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Println(solution)
	fmt.Println()

	ax, _ := matrix.MultOnVecRight(degen_m.A, solution)
	fmt.Println("ax: ")
	fmt.Println(ax)
	fmt.Println()

	matrix.Print(degen_m.A, "degen_A")
	matrix.Print(degen_m.L, "degen_L")
	matrix.Print(degen_m.U, "degen_U")
	matrix.Print(degen_m.P, "degen_P")
	matrix.Print(degen_m.Q, "degen_Q")
}
