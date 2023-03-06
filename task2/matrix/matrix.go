package matrix

import (
	"errors"
	"fmt"
	"math"
	"math/rand"
	"time"
)

var Eps float64 = math.Pow10(-14)

type Matrix struct {
	rows, cols    int
	A, L, U, P, Q [][]float64
	rowExchanges  int
	isLU bool
}

func New(rows, cols int) *Matrix {
	matrix := &Matrix{
		rows: rows,
		cols: cols,
		A:    make([][]float64, rows),
		L:    make([][]float64, rows),
		U:    make([][]float64, rows),
		P:    make([][]float64, rows),
		Q:    make([][]float64, rows),
		isLU: false,
	}

	for i := range matrix.A {
		matrix.A[i] = make([]float64, cols)
		matrix.L[i] = make([]float64, cols)
		matrix.U[i] = make([]float64, cols)
		matrix.P[i] = make([]float64, cols)
		matrix.Q[i] = make([]float64, cols)
	}

	rand.Seed(time.Now().UnixNano())

	// Initialize A matrix
	for i, r := range matrix.A {
		for j := range r {
			matrix.A[i][j] = rand.Float64() * 100
		}
	}

	if rows == cols {
		// Initialize P, Q, L matrix
		for i := range matrix.P {
			matrix.P[i][i] = 1.
			matrix.Q[i][i] = 1.
		}
	}

	return matrix
}

func Copy(m [][]float64) [][]float64 {
	res := NewMatrix(len(m), len(m[0]))

	for i, r := range m {
		for j := range r {
			res[i][j] = m[i][j]
		}
	}

	return res
}

func Mult(A, B [][]float64) [][]float64 {
	m, n, p := len(A), len(A[0]), len(B[0])
	C := NewMatrix(m, p)

	for i := 0; i < m; i++ {
		for j := 0; j < p; j++ {
			for k := 0; k < n; k++ {
				C[i][j] += A[i][k] * B[k][j]
			}
		}
	}

	return C
}

func Print(m [][]float64, name string) {
	fmt.Printf("matrix %s: \n", name)
	for i, r := range m {
		for j := range r {
			fmt.Printf("%.10f ", m[i][j])
		}
		fmt.Println()
	}
	fmt.Println()
}

func NewMatrix(rows, cols int) [][]float64 {
	res := make([][]float64, rows)
	for i := range res {
		res[i] = make([]float64, cols)
	}
	return res
}

func Norm(m [][]float64) float64 {
	norm := 0.

	for i := range m {
		rowSum := 0.
		for j := range m[i] {
			rowSum += math.Abs(m[i][j])
		}
		if rowSum > norm {
			norm = rowSum
		}
	}

	return norm
}

func GetRandomVector(size int) []float64 {
	res := make([]float64, size)

	// Fill with random numbers
	for i := range res {
		res[i] = rand.Float64() * 100
	}

	return res
}

func MultOnVecRight(matrix [][]float64, vec []float64) ([]float64, error) {
	rows := len(matrix)
	cols := len(matrix[0])
	if len(vec) != cols {
		return nil, errors.New("vector length does not match matrix columns")
	}

	res := make([]float64, rows)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			res[i] += vec[j] * matrix[i][j]
		}
	}

	return res, nil
}

func MultOnVecLeft(vec []float64, matrix [][]float64) ([]float64, error) {
	rows := len(matrix)
	cols := len(matrix[0])
	if len(vec) != rows {
		return nil, errors.New("vector length does not match matrix rows")
	}

	res := make([]float64, cols)
	for j := 0; j < cols; j++ {
		for i := 0; i < rows; i++ {
			res[j] += vec[i] * matrix[i][j]
		}
	}

	return res, nil
}

func PrintVec(vec []float64, name string) {
	fmt.Printf("%s: \n", name)
	for i := range vec {
		fmt.Printf("%.2f ", vec[i])
	}
	fmt.Println()
	fmt.Println()
}

func EqualToA(A, matrix [][]float64, eps float64) bool {
	if len(A) != len(matrix) || len(A[0]) != len(matrix[0]) {
		return false
	}

	for i, r := range A {
		for j := range r {
			if math.Abs(A[i][j]-matrix[i][j]) > eps {
				return false
			}
		}
	}

	return true
}

func EqualVecs(vec1, vec2 []float64, eps float64) bool {
	if len(vec1) != len(vec2) {
		return false
	}

	for i := range vec1 {
		if math.Abs(vec1[i]-vec2[i]) > eps {
			return false
		}
	}

	return true
}

func MultVecs(vec1, vec2 []float64) ([]float64, error) {
	if len(vec1) != len(vec2) {
		return nil, errors.New("lengths don't equal")
	}

	res := make([]float64, len(vec1))
	for i := range vec1 {
		res[i] = vec1[i] * vec2[i]
	}

	return res, nil
}

func IsEye(matrix [][]float64, eps float64) bool {
	if len(matrix) != len(matrix[0]) {
		return false
	}

	for i := range matrix {
		for j := range matrix[0] {
			num := 0.0
			if i == j {
				num = 1.0
			}
			if math.Abs(matrix[i][j] - num) > eps {
				return false
			}
		}
	}
	return true
}
