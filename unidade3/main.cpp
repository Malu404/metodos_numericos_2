#include <iostream>
#include <iomanip>
#include "utilitarios.h"
using namespace std;

// Declaracao das funcoes definidas nos outros arquivos
pair<Matrix, Matrix> mhausholder(Matrix A, int n);
pair<Matrix, Vector> metodoQRcode(const Matrix& A, int n, double epsilon);

// Funcao para imprimir autovalores e autovetores
void imprimirResultado(const Matrix& A_original, const Matrix& P, const Vector& autovalores) {
    int n = A_original.size();
    cout << fixed << setprecision(6);

    cout << "\nAutovalores encontrados: ";
    for (double val : autovalores) cout << val << " ";
    cout << "\n\nAutovetores (colunas de P):" << endl;
    for (int i = 0; i < n; i++) {
        cout << "v" << i+1 << ": ";
        for (int j = 0; j < n; j++) cout << P[j][i] << " ";
        cout << endl;
    }

    cout << "\nVerificacao A*v = lambda*v:\n";
    for (int i = 0; i < n; i++) {
        Vector vi(n);
        for (int j = 0; j < n; j++) vi[j] = P[j][i];
        Vector Avi = multiplicarMatrizVetor(A_original, vi);

        cout << "lambda = " << autovalores[i] << endl;
        cout << "A*v: ";
        for (double val : Avi) cout << val << " ";
        cout << "\nlambda*v: ";
        for (int j = 0; j < n; j++) cout << autovalores[i] * vi[j] << " ";

        double erro = 0.0;
        for (int j = 0; j < n; j++) erro += pow(Avi[j] - autovalores[i] * vi[j], 2);
        cout << "\nErro ||A*v - lambda*v||: " << sqrt(erro) << "\n" << endl;
    }
}

int main() {
    Matrix A = {
        {5.0, 2.0, 1.0},
        {2.0, 3.0, 1.0},
        {1.0, 1.0, 1.0}
    };
    int n = A.size();
    double epsilon = 1e-6;

    cout << "\nEtapa 1: Tridiagonalizacao via Householder\n";
    auto [T, H] = mhausholder(A, n);

    cout << "\nEtapa 2: Aplicando QR na matriz tridiagonal\n";
    auto [P, autovalores] = metodoQRcode(T, n, epsilon);

    imprimirResultado(A, multiplicarMatrizes(H, P), autovalores);

    return 0;
}
