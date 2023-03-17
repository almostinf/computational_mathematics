package matrix

import (
	"errors"
	"math"
)

// Use Givens method
// Here is explanation: https://www.tspi.at/2021/12/08/qrgivens.html
func (m *Matrix) QR() error {
    if m.isQR {
        return errors.New("QR decomposition is also done")
    }
    m.isQR = true

	m.Q = Eye(m.rows)
	G := Eye(m.rows)
	for j := 0; j < m.cols; j++ {
		for i := m.rows - 1; i > j; i-- {
			if math.Abs(m.R[i][j]) < Eps*Norm(m.A) {
				continue
			}

			c := m.R[j][j] / math.Sqrt(m.R[j][j]*m.R[j][j]+m.R[i][j]*m.R[i][j])
			s := m.R[i][j] / math.Sqrt(m.R[j][j]*m.R[j][j]+m.R[i][j]*m.R[i][j])

			G[j][j], G[i][i] = c, c
			G[i][j], G[j][i] = -s, s

			m.Q = Mult(G, m.Q)
			m.R = Mult(G, m.R)

			G[j][j], G[i][i] = 1, 1
			G[i][j], G[j][i] = 0, 0
		}
	}

    m.Q = Transpose(m.Q)

    return nil
}

func (m *Matrix) SolveQR(b []float64) ([]float64, error) {
    if !m.isQR {
        m.QR()
    }

    QT := Transpose(m.Q)

    Qb, err := MultOnVecRight(QT, b)
    if err != nil {
        return nil, err
    }

    x := make([]float64, m.cols)
    for i := m.rows - 1; i >= 0; i-- {
        x[i] = Qb[i]
        for j := i + 1; j < m.cols; j++ {
            x[i] -= m.R[i][j] * x[j]
        }
        x[i] /= m.R[i][i]
    }
    return x, nil
}
