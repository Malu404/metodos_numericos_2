#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "../utilitarios.h"
using namespace std;

// Tipo para matrizes e vetores
typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

// Função para criar matriz de Householder H(j)
Matrix hausholder(const Matrix& A, int j, int n) {
    // comeco com meu vetor de zeros
    Vector v(n, 0.0);
    Vector v_linha(n, 0.0);
    
    // v((j+1):n) <- A((j+1):n , j)
    for (int i = j + 1; i < n; i++) {
        v[i] = A[i][j];
    }
    
    // normalizo meu vetor
    double Lv = norma(v);
    
    // v'(j+1) <- Lv
    v_linha[j + 1] = Lv;
    
    // N <- v - v'
    Vector N(n);
    for (int i = 0; i < n; i++) {
        N[i] = v[i] - v_linha[i];
    }
    
    // n <- normalize(N)
    Vector n_normalizado = normalize(N);
    
    //Tendo em mente que a matriz HH e :  H = I - 2*n*n^T , inicializo com identidade
    Matrix H = identidade(n);
    //depois tiro o 2*n*n^T de H
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            H[i][k] -= 2.0 * n_normalizado[i] * n_normalizado[k];
        }
    }
    
    return H;
}

//funcao para as transformacoes de hausholder
pair<Matrix, Matrix> mhausholder(Matrix A, int n) {
    Matrix H = identidade(n);
    
    // loop (j = 1:(n-2))
    for (int j = 1; j <= n - 2; j++) {
        // H(j) <- hausholder(A,j,n)
        Matrix Hj = hausholder(A, j, n);
        
        // A <- H(j)^T * A * H(j)
        Matrix HjT = transpor(Hj);
        Matrix temp = multiplicarMatrizes(HjT, A);
        A = multiplicarMatrizes(temp, Hj);
        
        // H <- H * H(j)
        H = multiplicarMatrizes(H, Hj);
    }
    
    return make_pair(A, H);
}

// print
void imprimirMatriz(const Matrix& mat, const string& nome) {
    cout << nome << ":" << endl;
    for (const auto& linha : mat) {
        for (double val : linha) {
            cout << setw(10) << setprecision(4) << fixed << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int main() {// inserir uma matriz aqui, de preferencia as que eu ja usei nos metodos anteriores pq ja tenho elas calculadas
    // Exemplo de uso
    int n;
    cout << "Digite o tamanho da matriz: ";
    cin >> n;
    
    Matrix A(n, Vector(n));
    
    cout << "Digite os elementos da matriz A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << "A[" << i << "][" << j << "]: ";
            cin >> A[i][j];
        }
    }
    
    cout << "\nMatriz original:" << endl;
    imprimirMatriz(A, "A");
    
    // aplicar transformacoes HH
    pair<Matrix, Matrix> resultado = mhausholder(A, n);
    Matrix T = resultado.first;  // Matriz tridiagonalizada
    Matrix H = resultado.second; // Matriz de transformação acumulada
    
    cout << "Resultado das transformações de Householder:" << endl;
    imprimirMatriz(T, "Matriz Tridiagonalizada (T)");
    imprimirMatriz(H, "Matriz de Transformação (H)");
    
    return 0;
}