package main

import (
	"github.com/almostinf/computational_mathematics/task3/newton"
)

func main() {
	n := newton.New()
	n.SolveSystem()
	n.ModifiedSolveSystem()
	n.SolveSystemWithKIterations(7)
	n.MethodsTransition(7, 7)
}
