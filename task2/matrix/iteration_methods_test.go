package matrix

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestJacobiOnSLAE(t *testing.T) {
	m := GetMatrixDiagonalDominance(N)

	b := GetRandomVector(N)

	x, _, _ := m.JacobiMethod(b)
	Ax, err := MultOnVecRight(m.A, x)
	require.NoError(t, err)

	assert.True(t, EqualVecs(Ax, b, Eps*Norm(m.A)))
}

func TestSeidelOnSLAE(t *testing.T) {
	m := GetMatrixDiagonalDominance(N)

	b := GetRandomVector(N)

	x, _, _ := m.SeidelMethod(b)
	Ax, err := MultOnVecRight(m.A, x)
	require.NoError(t, err)

	assert.True(t, EqualVecs(Ax, b, Eps*Norm(m.A)))
}
