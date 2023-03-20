package matrix

import "math"

func GetBc(A [][]float64, b []float64) ([][]float64, []float64) {
	n := len(A)
	l := Zeros(n)
	invD := Zeros(n)
	r := Zeros(n)

	for i := range invD {
		invD[i][i] = 1.0 / A[i][i]
	}

	for i := 0; i < n; i++ {
		for j := 0; j < i; j++ {
			l[i][j] = A[i][j]
		}
	}

	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			r[i][j] = A[i][j]
		}
	}

	B := Negative(Mult(invD, Sum(l, r)))
	c, _ := MultOnVecRight(invD, b)
	return B, c
}

func (m *Matrix) JacobiMethod(b []float64) ([]float64, int, int) {
	B, c := GetBc(m.A, b)
	prevX := make([]float64, m.rows)
	copy(prevX, c)

	Bx, _ := MultOnVecRight(B, prevX)
	x := SumVectors(Bx, c)

	q := Norm(B)
	var cond float64
	if q > 1 {
		cond = 0.01
	} else {
		cond = (1 - q) / q
	}

	prioriEstimate := (math.Log10(Eps) - math.Log10(NormVector(c)) + math.Log10((1 - q))) / math.Log10(q)

	iterations := 0

	for NormVector(DiffVectors(x, prevX)) > cond*Eps*Norm(m.A) {
		if iterations == 1000 {
			return x, iterations, int(prioriEstimate)
		}

		copy(prevX, x)
		Bx, _ = MultOnVecRight(B, x)
		x = SumVectors(Bx, c)
		iterations++
	}

	return x, iterations, int(prioriEstimate)
}

func (m *Matrix) SeidelMethod(b []float64) ([]float64, int, int) {
	B, c := GetBc(m.A, b)
	prevX := make([]float64, m.rows)
	copy(prevX, c)

	Bx, _ := MultOnVecRight(B, prevX)
	x := SumVectors(Bx, c)

	q := Norm(B)
	var cond float64
	if q > 1 {
		cond = 0.01
	} else {
		cond = (1 - q) / q
	}

	prioriEstimate := (math.Log10(Eps) - math.Log10(NormVector(c)) + math.Log10((1 - q))) / math.Log10(q)

	iterations := 0

	for NormVector(DiffVectors(x, prevX)) > cond*Eps*Norm(m.A) {
		if iterations == 100000 {
			return x, iterations, int(prioriEstimate)
		}

		copy(prevX, x)

		for i := 0; i < m.rows; i++ {
			xk1 := c[i]
			for j := 0; j < m.cols; j++ {
				xk1 += x[j] * B[i][j]
			}
			x[i] = xk1
		}

		iterations++
	}

	return x, iterations, int(prioriEstimate)
}
