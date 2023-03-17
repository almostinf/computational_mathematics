package matrix

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestQRDecomposition(t *testing.T) {
	m := New(N, N)
	err := m.QR()
	require.NoError(t, err)

	qr := Mult(m.Q, m.R)
	assert.True(t, EqualMatrixes(m.A, qr, Eps * Norm(m.A)))
}

func TestQRSLAESolution(t *testing.T) {
	m := New(N, N)
	b := GetRandomVector(N)

	x, err := m.SolveQR(b)
	require.NoError(t, err)

	Ax, err := MultOnVecRight(m.A, x)
	require.NoError(t, err)

	assert.True(t, EqualVecs(Ax, b, Eps * Norm(m.A)))
}