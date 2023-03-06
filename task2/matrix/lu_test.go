package matrix

import (
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
	A_fixed := Mult(m.P, m.A) // removing permutation effect

	assert.True(t,EqualToA(A_fixed, lu, Eps * Norm(m.A)))
}

func TestSLAE(t *testing.T) {
	m := New(N, N)

	err := m.LUDecomposition();
	require.NoError(t, err)

	b := GetRandomVector(N)
	b_fixed, err := MultOnVecRight(m.P, b) // Lx = Pb
	require.NoError(t, err)

	x, err := m.SLAESolution(b_fixed)
	require.NoError(t, err)

	Ax, err := MultOnVecRight(m.A, x)
	require.NoError(t, err)

	assert.True(t, EqualVecs(Ax, b, Eps * Norm(m.A)))
}

func TestInverseLeft(t *testing.T) {
	m := New(N, N)

	inv, err := m.Inverse()
	require.NoError(t, err)

	eye := Mult(inv, m.A)

	assert.True(t, IsEye(eye, Eps * Norm(m.A)))
}

func TestInverseRight(t *testing.T) {
	m := New(N, N)

	inv, err := m.Inverse()
	require.NoError(t, err)

	eye := Mult(m.A, inv)

	assert.True(t, IsEye(eye, Eps * Norm(m.A)))
}
