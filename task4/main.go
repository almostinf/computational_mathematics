package main

import (
	"fmt"
	"math"

	"github.com/almostinf/computational_mathematics/task4/integral"
)

const (
	eps           = 1e-6
	expectedValue = 20.73027110955
)

func getAnswerWithAccuracy(a, b, alpha, step float64, err *float64, f func(x float64) float64, method func(low, high, b, step, alpha float64, f func(x float64) float64) float64) float64 {
	high := 0.
	ans := 0.
	l := 2.
	deg := 4.

	for *err > eps {
		high = a + step
		ans = method(a, high, b, step, alpha, f)

		step /= l
		high = a + step
		ans2 := method(a, high, b, step, alpha, f)

		step /= l
		high = a + step
		ans3 := method(a, high, b, step, alpha, f)
		
		speed := -math.Log(math.Abs((ans3-ans2)/(ans2-ans))) / math.Log(l)

		*err = math.Abs((ans2 - ans) / (math.Pow(l, deg) - 1))
		ans += (ans2 - ans) / (1 - math.Pow(l, -deg))
		fmt.Println("Answer: ", ans)
		fmt.Println("Expected: ", expectedValue)
		fmt.Println("Mistake: ", math.Abs(ans-expectedValue))
		fmt.Println("Error: ", *err)
		fmt.Println("Speed: ", speed)
		fmt.Println()
	}

	return ans
}

func main() {
	// 9 variant
	a := 2.5
	b := 3.3
	alpha := 2. / 3.

	f := func(x float64) float64 {
		return 3*math.Cos(1.5*x)*math.Exp(x/4) + 4*math.Sin(3.5*x)*math.Exp(-3*x) + 4*x
	}

	fmt.Println("-----------------------------------------------------")
	fmt.Println("Newton Cots: ")
	ans := integral.NewtonCotes(a, b, b, 1, alpha, f)
	fmt.Println("Ans: ", ans)
	fmt.Println("Expected: ", expectedValue)
	fmt.Println("Mistake: ", math.Abs(ans-expectedValue))
	fmt.Println("-----------------------------------------------------")
	fmt.Println()
	fmt.Println("-----------------------------------------------------")
	fmt.Println("Newton Cots modified: ")
	step := (b - a) / 10
	l := 2.
	deg := 4.
	err := 1.

	ans = getAnswerWithAccuracy(a, b, alpha, step, &err, f, integral.NewtonCotes)

	fmt.Println("Answer: ", ans)
	fmt.Println("Expected: ", expectedValue)
	fmt.Println("Mistake: ", math.Abs(ans-expectedValue))
	fmt.Println("-----------------------------------------------------")
	fmt.Println()
	fmt.Println("-----------------------------------------------------")
	fmt.Println("Finding hOpt for Newton Cots: ")
	hOpt := (b - a) / math.Ceil((b-a)/(step*l*math.Pow((eps/err), 1./deg)))

	ans = getAnswerWithAccuracy(a, b, alpha, hOpt, &err, f, integral.NewtonCotes)

	fmt.Println("hOpt: ", hOpt)
	fmt.Println("Answer: ", ans)
	fmt.Println("Expected: ", expectedValue)
	fmt.Println("Mistake: ", math.Abs(ans-expectedValue))
	fmt.Println("-----------------------------------------------------")
	fmt.Println()
	fmt.Println("-----------------------------------------------------")
	fmt.Println("Gauss: ")

	ans = integral.Gauss(a, b, b, 1, alpha, f)

	fmt.Println("Ans: ", ans)
	fmt.Println("Expected: ", expectedValue)
	fmt.Println("Mistake: ", math.Abs(ans-expectedValue))
	fmt.Println("-----------------------------------------------------")
	fmt.Println()
	fmt.Println("-----------------------------------------------------")
	fmt.Println("Gauss modified: ")
	step = (b - a) / 10
	err = 1.

	ans = getAnswerWithAccuracy(a, b, alpha, step, &err, f, integral.Gauss)

	fmt.Println("Answer: ", ans)
	fmt.Println("Expected: ", expectedValue)
	fmt.Println("Mistake: ", math.Abs(ans-expectedValue))
	fmt.Println("-----------------------------------------------------")
	fmt.Println()
	fmt.Println("-----------------------------------------------------")
	fmt.Println("Finding hOpt for Gauss Cots: ")
	hOpt = (b - a) / math.Ceil((b-a)/(step*l*math.Pow((eps/err), 1./deg)))
	err = 1.

	ans = getAnswerWithAccuracy(a, b, alpha, hOpt, &err, f, integral.Gauss)

	fmt.Println("hOpt: ", hOpt)
	fmt.Println("Answer: ", ans)
	fmt.Println("Expected: ", expectedValue)
	fmt.Println("Mistake: ", math.Abs(ans-expectedValue))
	fmt.Println("-----------------------------------------------------")
}
