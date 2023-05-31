package integral

import (
	"log"
	"math"
	"sort"

	"github.com/almostinf/computational_mathematics/task2/matrix"
)

func getMoments(a, alpha, z1, z2 float64) []float64 {
	moments := make([]float64, 6)
	moments[0] = (math.Pow(z1-a, 1.-alpha) - math.Pow(z2-a, 1.-alpha)) / (1. - alpha)
	moments[1] = (math.Pow(z1-a, 2.-alpha)-math.Pow(z2-a, 2.-alpha))/(2.-alpha) + a*moments[0]
	moments[2] = (math.Pow(z1-a, 3.-alpha)-math.Pow(z2-a, 3.-alpha))/(3.-alpha) + 2.*a*moments[1] - a*a*moments[0]
	moments[3] = (math.Pow(z1-a, 4.-alpha)-math.Pow(z2-a, 4.-alpha))/(4.-alpha) + 3.*a*moments[2] - 3*a*a*moments[1] + a*a*a*moments[0]
	moments[4] = (math.Pow(z1-a, 5.-alpha)-math.Pow(z2-a, 5.-alpha))/(5.-alpha) + 4.*a*moments[3] - 6*a*a*moments[2] + 4*a*a*a*moments[1] - a*a*a*a*moments[0]
	moments[5] = (math.Pow(z1-a, 6.-alpha)-math.Pow(z2-a, 6.-alpha))/(6.-alpha) + 5.*a*moments[4] - 10*a*a*moments[3] + 10*a*a*a*moments[2] - 5*a*a*a*a*moments[1] + a*a*a*a*a*moments[0]

	return moments
}

func NewtonCotes(low, high, b, step, alpha float64, f func(x float64) float64) float64 {
	var ans float64
	z1, z2, z3 := low, (low+high)/2, high
	n := int(math.Ceil((b-high)/step) + 1)

	for i := 1; i <= n; i++ {
		moments := getMoments(low, alpha, z3, z1)
		neededMoments := []float64{moments[0], moments[1], moments[2]}

		A := matrix.New(3, 3)
		A.A[0][0] = 1
		A.A[0][1] = 1
		A.A[0][2] = 1
		A.A[1][0] = z1
		A.A[1][1] = z2
		A.A[1][2] = z3
		A.A[2][0] = z1 * z1
		A.A[2][1] = z2 * z2
		A.A[2][2] = z3 * z3

		err := A.LUDecomposition()
		if err != nil {
			log.Fatal("Newton Cotes LU: ", err)
		}

		fixedMoments, err := matrix.MultOnVecRight(A.P, neededMoments)
		if err != nil {
			log.Fatal("Newton Cots mult on matrix P: ", err)
		}

		s, err := A.SLAESolution(fixedMoments)
		if err != nil {
			log.Fatal("Newton Cots SLAE Solution: ", err)
		}

		b := []float64{z1, z2, z3}
		for j := 0; j < 3; j++ {
			ans += s[j] * f(b[j])
		}

		z1 = low + float64(i)*step
		z3 = low + float64(i+1)*step
		z2 = (z1 + z3) / 2
	}

	return math.Abs(ans)
}

func NewtonCotesLim(a, b, step, alpha float64, f func(x float64) float64) float64 {
	var ans float64
	lim1 := a
	lim2 := lim1 + step
	n := int(math.Round((b - lim2) / step))

	for i := 0; i < n; i++ {
		moments := getMoments(a, alpha, lim2, lim1)
		neededMoments := []float64{moments[0], moments[1], moments[2]}

		nodes := []float64{lim1, (lim1 + lim2) / 2, lim2}
		A := matrix.New(3, 3)
		A.A[0][0] = 1
		A.A[0][1] = 1
		A.A[0][2] = 1
		A.A[1][0] = nodes[0]
		A.A[1][1] = nodes[1]
		A.A[1][2] = nodes[2]
		A.A[2][0] = nodes[0] * nodes[0]
		A.A[2][1] = nodes[1] * nodes[1]
		A.A[2][2] = nodes[2] * nodes[2]

		err := A.LUDecomposition()
		if err != nil {
			log.Fatal("Newton Cotes LU: ", err)
		}

		fixedMoments, err := matrix.MultOnVecRight(A.P, neededMoments)
		if err != nil {
			log.Fatal("Newton Cots mult on matrix P: ", err)
		}

		s, err := A.SLAESolution(fixedMoments)
		if err != nil {
			log.Fatal("Newton Cots SLAE Solution: ", err)
		}

		for j := 0; j < 3; j++ {
			ans += s[j] * f(nodes[j])
		}

		lim1 += step
		lim2 += step
	}

	return ans
}

func kardano(a, x []float64) []float64 {
	a[0], a[2] = a[2], a[0]

	p := a[1] - a[0]*a[0]/3.0
	q := a[2] + 2.0*a[0]*a[0]*a[0]/27.0 - a[0]*a[1]/3.0
	determinant := q*q/4.0 + p*p*p/27.0

	if determinant < 0 {
		fi := 0.0
		if q < 0 {
			fi = math.Atan(2.0 * math.Sqrt(-determinant) / (-q))
		}
		if q > 0 {
			fi = math.Atan(2.0*math.Sqrt(-determinant)/(-q) + math.Pi)
		}
		if q == 0 {
			fi = math.Pi / 2.0
		}

		x[0] = 2.0*math.Sqrt(-p/3.0)*math.Cos(fi/3.0) - a[0]/3.0
		x[1] = 2.0*math.Sqrt(-p/3.0)*math.Cos(fi/3.0+2.0*math.Pi/3.0) - a[0]/3.0
		x[2] = 2.0*math.Sqrt(-p/3.0)*math.Cos(fi/3.0+4.0*math.Pi/3.0) - a[0]/3.0
	} else if determinant > 0 {
		x[1] = 0.0
		if (-q)/2.0+math.Pow(determinant, 1.0/2.0) < 0 {
			x[1] += -math.Pow((-q)/2.0-math.Pow(determinant, 1.0/2.0), 1.0/3.0)
		} else {
			x[1] += math.Pow((-q)/2.0+math.Pow(determinant, 1.0/2.0), 1.0/3.0)
		}
		if -q/2.0-math.Pow(determinant, 1.0/2.0) < 0 {
			x[1] += -math.Pow(q/2.0+math.Pow(determinant, 1.0/2.0), 1.0/3.0) - a[0]/3.0
		} else {
			x[1] += math.Pow(-q/2.0-math.Pow(determinant, 1.0/2.0), 1.0/3.0) - a[0]/3.0
		}
	} else if determinant == 0 {
		x[0] = 2*math.Pow(-q/2.0, 1.0/3.0) - a[0]/3.0
		x[1] = -math.Pow(-q/2.0, 1.0/3.0) - a[0]/3.0
	}

	return x
}

func Gauss(low, high, b, step, alpha float64, f func(x float64) float64) float64 {
	var ans float64
	z1, z2, z3 := low, (low+high)/2, high
	n := int(math.Ceil((b-high)/step) + 1)

	for i := 1; i <= n; i++ {
		m := matrix.New(3, 3)

		moments := getMoments(low, alpha, z3, z1)

		m.A[0][0] = moments[0]
		m.A[0][1] = moments[1]
		m.A[0][2] = moments[2]
		m.A[1][0] = moments[1]
		m.A[1][1] = moments[2]
		m.A[1][2] = moments[3]
		m.A[2][0] = moments[2]
		m.A[2][1] = moments[3]
		m.A[2][2] = moments[4]

		b := []float64{-moments[3], -moments[4], -moments[5]}

		err := m.LUDecomposition()
		if err != nil {
			log.Fatal("Gauss LU: ", err)
		}

		fixedB, err := matrix.MultOnVecRight(m.P, b)
		if err != nil {
			log.Fatal("Gauss mult on matrix P: ", err)
		}

		a, err := m.SLAESolution(fixedB)
		if err != nil {
			log.Fatal("Gauss: ", err)
		}

		coef := []float64{z1, z2, z3}

		sol := kardano(a, coef)

		sort.Float64s(sol)

		A := matrix.New(3, 3)
		A.A[0][0] = 1
		A.A[0][1] = 1
		A.A[0][2] = 1
		A.A[1][0] = sol[0]
		A.A[1][1] = sol[1]
		A.A[1][2] = sol[2]
		A.A[2][0] = sol[0] * sol[0]
		A.A[2][1] = sol[1] * sol[1]
		A.A[2][2] = sol[2] * sol[2]

		err = A.LUDecomposition()
		if err != nil {
			log.Fatal("Gauss LU: ", err)
		}

		neededMoments := []float64{moments[0], moments[1], moments[2]}

		fixedMoments, err := matrix.MultOnVecRight(A.P, neededMoments)
		if err != nil {
			log.Fatal("Gauss mult on matrix P: ", err)
		}

		s, err := A.SLAESolution(fixedMoments)
		if err != nil {
			log.Fatal("Gauss SLAE: ", err)
		}

		for j := 0; j < 3; j++ {
			ans += s[j] * f(coef[j])
		}

		z1 = low + float64(i)*step
		z3 = low + float64(i+1)*step
		z2 = (z1 + z3) / 2
	}

	return ans
}
