package matrix

import (
	"errors"
	"math"
)

func (m *Matrix) LUDecomposition() error {
	if m.isLU {
		return nil
	}
	m.isLU = true

	// Check if A is square
	if m.rows != m.cols {
		return errors.New("matrix A is not square")
	}

	// Initialize U matrix
	for i, r := range m.U {
		for j := range r {
			m.U[i][j] = m.A[i][j]
		}
	}

	for i := 0; i < m.rows; i++ {
		max := m.U[i][i]
		pivot := i
		for j := i + 1; j < m.rows; j++ {
			if math.Abs(m.U[j][i]) > max && m.U[j][i] != 0 {
				max = m.U[j][i]
				pivot = j
			}
		}

		if pivot != i {
			m.U[i], m.U[pivot] = m.U[pivot], m.U[i]
			m.P[i], m.P[pivot] = m.P[pivot], m.P[i]
			m.L[i], m.L[pivot] = m.L[pivot], m.L[i]
			m.rowExchanges++
		}

		for k := i + 1; k < m.rows; k++ {
			if m.U[i][i] != 0 {
				coef := m.U[k][i] / m.U[i][i]
				for j := i; j < m.cols; j++ {
					m.U[k][j] -= m.U[i][j] * coef
				}
				m.L[k][i] = coef
			}
		}

		m.L[i][i] = 1.0
	}

	return nil
}

func (m *Matrix) Determinant() (float64, error) {
	// Check if A is square
	if m.rows != m.cols {
		return 0, errors.New("matrix A is not square")
	}

	// Perform LU decomposition
	err := m.LUDecomposition()
	if err != nil {
		return 0, err
	}

	// Compute determinant
	det := 1.0
	for i := 0; i < m.rows; i++ {
		det *= m.U[i][i]
	}

	if m.rowExchanges%2 != 0 { // TODO with sign
		det *= -1.
	}

	return det, nil
}

func (m *Matrix) Inverse() ([][]float64, error) {
	// Perform LU decomposition
	if err := m.LUDecomposition(); err != nil {
		return nil, err
	}

	inverse := NewMatrix(m.rows, m.cols)

	// Solve for inverse using forward/back substitution
	for i := 0; i < m.rows; i++ {
		// Solve L*y = I
		y := make([]float64, m.rows)
		// y[i] = 1

		for j := 0; j < m.rows; j++ {
			y[j] = m.P[j][i]
		}

		// y_fixed, err := MultOnVecRight(m.P, y)
		// if err != nil {
		// 	return nil, err
		// }

		x, err := m.SLAESolution(y)
		if err != nil {
			return nil, err
		}

		// Set the i-th column of the inverse matrix
		for j := 0; j < m.rows; j++ {
			inverse[j][i] = x[j]
		}
	}

	return inverse, nil
}

func (m *Matrix) ConditionNumber() (float64, error) {
	if err := m.LUDecomposition(); err != nil {
		return -1., err
	}

	inverse, err := m.Inverse()
	if err != nil {
		return -1., err
	}

	return Norm(inverse) * Norm(m.A), nil
}

func (m *Matrix) SLAESolution(b []float64) ([]float64, error) {
	if err := m.LUDecomposition(); err != nil {
		return nil, err
	}

	// Solve for y in Ly = b using forward substitution
	y := make([]float64, m.rows)
	for i := 0; i < m.rows; i++ {
		y[i] = b[i]
		for j := 0; j < i; j++ {
			y[i] -= m.L[i][j] * y[j]
		}
	}

	// Forward substitution formula:
	// $$ x_m = \frac{b_m - \sum_{i=1}^{m-1}{l_{m,i} \cdot x_i}}{l_{m,m}} $$
	// wiki: https://en.wikipedia.org/wiki/Triangular_matrix#Forward_and_back_substitution

	// Solve for x in Ux = y using backward substitution
	x := make([]float64, m.rows)
	for i := m.rows - 1; i >= 0; i-- {
		x[i] = y[i]
		for j := i + 1; j < m.cols; j++ {
			x[i] -= m.U[i][j] * x[j]
		}
		x[i] /= m.U[i][i]
	}

	return x, nil
}