package main

import (
	"fmt"

	"github.com/almostinf/computational_mathematics/task3/newton"
)

func main() {
	n := newton.New()

	fmt.Println("FIRST TASK")
	n.SolveSystem()

	fmt.Println("SECOND TASK")
	n.ModifiedSolveSystem()

	fmt.Println("THIRD TASK")
	n.ModifiedSolveSystemOnlyKIterations(3)

	fmt.Println("FOURTH TASK")
	n.SolveSystemWithMIterations(7)

	fmt.Println("FIFTH TASK")
	n.MethodsTransition(7, 7)

	fmt.Println("SIXTH TASK")
	fmt.Println("k < 7")
	n.ModifiedSolveSystemOnlyKIterations(4)
	fmt.Println("k == 7")
	n.ModifiedSolveSystemOnlyKIterations(6)
	fmt.Println("k > 7")
	n.ModifiedSolveSystemOnlyKIterations(7)
}
