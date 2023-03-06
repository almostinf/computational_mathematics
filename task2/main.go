package main

import (
	"fmt"

	"github.com/almostinf/computational_mathematics/task2/matrix"
)

const N int = 3

func main() {
	m := matrix.New(N, N)
	if err := m.LUDecomposition(); err != nil {
		fmt.Println(err)
		return
	}
	matrix.Print(m.A, "A")
	matrix.Print(m.L, "L")
	matrix.Print(m.U, "U")
	matrix.Print(m.P, "P")
	matrix.Print(m.Q, "Q")
	matrix.Print(matrix.Mult(m.L, m.U), "Mult")

	res, err := m.Determinant()
	if err != nil {
		fmt.Println(err)
		return
	}

	fmt.Println("Determinant: ", res)
	fmt.Println()

	inv, _ := m.Inverse()
	matrix.Print(inv, "Inverse")

	matrix.Print(matrix.Mult(inv, m.A), "Eye")

	b := matrix.GetRandomVector(N)
	matrix.PrintVec(b, "b")

	x, _ := m.SLAESolution(b)
	matrix.PrintVec(x, "x")

	fixed_A := matrix.Mult(m.P, m.A)

	new_b, _ := matrix.MultOnVecRight(fixed_A, x)
	matrix.PrintVec(new_b, "new_b")

	fixed_b, _ := matrix.MultOnVecRight(m.P, new_b)
	matrix.PrintVec(fixed_b, "fixed_b")

}
