package matrix

import (
	"errors"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

const N = 5

func TestLUOnEqualPAQ(t *testing.T) {
	m := New(N, N)

	err := m.LUDecomposition()
	require.NoError(t, err)

	lu := Mult(m.L, m.U)
	pa := Mult(m.P, m.A) // removing permutation effect

	assert.True(t, EqualMatrixes(pa, lu, Eps*Norm(m.A)))
}

func TestSLAE(t *testing.T) {
	m := New(N, N)

	err := m.LUDecomposition()
	require.NoError(t, err)

	b := GetRandomVector(N)
	b_fixed, err := MultOnVecRight(m.P, b) // Lx = Pb
	require.NoError(t, err)

	x, err := m.SLAESolution(b_fixed)
	require.NoError(t, err)

	Ax, err := MultOnVecRight(m.A, x)
	require.NoError(t, err)

	assert.True(t, EqualVecs(Ax, b, Eps*Norm(m.A)))
}

func TestInverseLeft(t *testing.T) {
	m := New(N, N)

	inv, err := m.Inverse()
	require.NoError(t, err)

	eye := Mult(inv, m.A)

	assert.True(t, IsEye(eye, Eps*Norm(m.A)))
}

func TestInverseRight(t *testing.T) {
	m := New(N, N)

	inv, err := m.Inverse()
	require.NoError(t, err)

	eye := Mult(m.A, inv)

	assert.True(t, IsEye(eye, Eps*Norm(m.A)))
}

// TESTS FOR DEGEN MATRIX
func TestSLAEonDegenMatrixNoSolutions(t *testing.T) {
	degenM := New(3, 3)
	matr := [][]float64{
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
		{7.0, 8.0, 9.0},
	}
	degenM.Set(matr)

	err := degenM.LUDecomposition()
	assert.NoError(t, err)

	b := []float64{1.0, 2.0, 5.0}
	b_fixed, err := MultOnVecRight(degenM.P, b)
	assert.NoError(t, err)

	expectedError := errors.New("system has no solutions")

	_, err = degenM.SLAESolution(b_fixed)
	assert.EqualError(t, err, expectedError.Error())
}

func TestSLAEonDegenMatrixPartialSolution(t *testing.T) {
	degenM := New(3, 3)
	matr := [][]float64{
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
		{7.0, 8.0, 9.0},
	}
	degenM.Set(matr)

	err := degenM.LUDecomposition()
	assert.NoError(t, err)

	b := []float64{1.0, 2.0, 3.0}
	b_fixed, err := MultOnVecRight(degenM.P, b)
	assert.NoError(t, err)

	x, err := degenM.SLAESolution(b_fixed)
	assert.NoError(t, err)

	Ax, err := MultOnVecRight(degenM.A, x)
	require.NoError(t, err)

	assert.True(t, EqualVecs(Ax, b, Eps*Norm(degenM.A)))
}
