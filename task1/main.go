package main

import (
	"fmt"
	"math"
)

var (
	eps1 = math.Pow10(-6) / 0.18
	eps2 = math.Pow10(-6) / 2.7
	eps3 = math.Pow10(-6) / 3
)

func sin(x float64, accuracy float64) float64 {
	// Start with the first term of the power series
	term := x
	sum := term

	// Keep calculating terms until the accuracy is reached
	for i := 1; math.Abs(term) >= accuracy; i++ {
		// Calculate the next term of the power series
		term = (-1) * term * x * x / float64((2*i)*(2*i+1))

		// Add the term to the sum
		sum += term
	}

	return sum
}

func sh(x float64, accuracy float64) float64 {
	// Start with the first term of the power series
	term := x
	sum := term

	// Keep calculating terms until the accuracy is reached
	for i := 1; math.Abs(term) >= accuracy; i++ {
		// Calculate the next term of the power series
		term = term * x * x / float64((2*i)*(2*i+1))

		// Add the term to the sum
		sum += term
	}

	return sum
}

func sqrt(x float64, accuracy float64) float64 {
	z := math.Max(x, 1)

	// Keep refining the guess until the desired accuracy is reached
	for math.Abs(z*z-x) >= accuracy {
		// Use the Geron formula to refine the guess
		z = (z + x/z) / 2
	}

	return z
}

func calculate(start, end, step float64) {
	for i := start; i <= end+step; i += step {
		my_calc := sqrt(sin(i+0.74, eps1), eps3) * sh(0.8*math.Pow(i, 2)+0.1, eps2)
		stand_calc := math.Sqrt(math.Sin(i+0.74)) * math.Sinh(0.8*math.Pow(i, 2)+0.1)
		err := math.Abs(my_calc - stand_calc)
		fmt.Printf("step = %.2g\t| my_calc = %.10g \t| stand_calc = %.10g\t| error = %.10g\n", i, my_calc, stand_calc, err)
	}
}

func main() {
	calculate(0.1, 0.2, 0.01)
}
