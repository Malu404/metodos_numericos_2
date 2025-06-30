#ifndef UTILITARIOS_H
#define UTILITARIOS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <utility>
using namespace std;

// Tipo para matrizes e vetores
typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

// Função para calcular a norma de um vetor
double norma(const Vector& v) {
    double sum = 0.0;
    for (double val : v) {
        sum += val * val;
    }
    return sqrt(sum);
}

// Função para normalizar um vetor
Vector normalize(const Vector& v) {
    double norm = norma(v);
    Vector normalized(v.size());
    if (norm != 0) {
        for (int i = 0; i < v.size(); i++) {
            normalized[i] = v[i] / norm;
        }
    }
    return normalized;
}

// Função para criar matriz identidade
Matrix identidade(int n) {
    Matrix I(n, Vector(n, 0.0));
    for (int i = 0; i < n; i++) {
        I[i][i] = 1.0;
    }
    return I;
}

// Função para multiplicar matrizes
Matrix multiplicarMatrizes(const Matrix& A, const Matrix& B) {
    int n = A.size();
    int m = B[0].size();
    int p = B.size();
    
    Matrix C(n, Vector(m, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < p; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// Função para multiplicar matriz por vetor
Vector multiplicarMatrizVetor(const Matrix& A, const Vector& v) {
    int n = A.size();
    Vector resultado(n, 0.0);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < v.size(); j++) {
            resultado[i] += A[i][j] * v[j];
        }
    }
    return resultado;
}

// Função para transpor uma matriz
Matrix transpor(const Matrix& A) {
    int n = A.size();
    int m = A[0].size();
    Matrix At(m, Vector(n));
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            At[j][i] = A[i][j];
        }
    }
    return At;
}

double produtoescalar(const Vector& a, const Vector& b){
    double resultado = 0.0;
    for (size_t i = 0; i<a.size(); i++){
        resultado += a[i] * b[i];
    }
    return resultado;
}

// Função para inverter uma matriz (Gauss-Jordan)
Matrix inversa(const Matrix& A) {
    int n = A.size();
    Matrix I = identidade(n);
    Matrix B = A;

    for (int i = 0; i < n; ++i) {
        double diag = B[i][i];
        if (fabs(diag) < 1e-10) {
            cerr << "Erro: pivo zero na inversão" << endl;
            exit(1);
        }

        // Normaliza a linha
        for (int j = 0; j < n; ++j) {
            B[i][j] /= diag;
            I[i][j] /= diag;
        }

        // Elimina os outros elementos da coluna
        for (int k = 0; k < n; ++k) {
            if (k == i) continue;
            double fator = B[k][i];
            for (int j = 0; j < n; ++j) {
                B[k][j] -= fator * B[i][j];
                I[k][j] -= fator * I[i][j];
            }
        }
    }

    return I;
}

// Função para calcular o método de potência regular
pair<double, Vector> Regularparainversa(const Matrix& A, const Vector& v0, double epsilon) {
    Vector vk_new = normalize(v0);
    Vector vk_old = v0;
    double lambda_new = 0.0;
    double lambda_old = 0.0;
    Vector x1;
    
    do {
        vk_old = vk_new;
        x1 = multiplicarMatrizVetor(A, vk_old);
        lambda_old = lambda_new;
        lambda_new = produtoescalar(x1, vk_old);
        vk_new = normalize(x1);
    } while (fabs(lambda_new - lambda_old) > epsilon);
    
    return{lambda_new,vk_new};
}

#endif // UTILITARIOS_H