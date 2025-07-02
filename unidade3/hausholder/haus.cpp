#include "../utilitarios.h"
#include <cmath>

Matrix hausholder(const Matrix& A, int j, int n) {
    Vector v(n, 0.0), v_linha(n, 0.0);
    for (int i = j + 1; i < n; i++) v[i] = A[i][j];

    double Lv = norma(v);
    v_linha[j + 1] = Lv;

    Vector N(n);
    for (int i = 0; i < n; i++) N[i] = v[i] - v_linha[i];

    Vector n_normalizado = normalize(N);

    Matrix H = identidade(n);
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            H[i][k] -= 2.0 * n_normalizado[i] * n_normalizado[k];

    return H;
}

pair<Matrix, Matrix> mhausholder(Matrix A, int n) {
    Matrix H = identidade(n);
    for (int j = 1; j <= n - 2; j++) {
        Matrix Hj = hausholder(A, j, n);
        A = multiplicarMatrizes(multiplicarMatrizes(transpor(Hj), A), Hj);
        H = multiplicarMatrizes(H, Hj);
    }
    return {A, H};
}
