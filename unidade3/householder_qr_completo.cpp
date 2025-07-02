#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "utilitarios.h"

using namespace std;

// Funcao de Householder para transformar matriz em forma tridiagonal
Matrix hausholder(const Matrix& A, int j, int n) {
    Vector v(n, 0.0), v_linha(n, 0.0);
    
    // Coletar elementos abaixo da diagonal
    for (int i = j + 1; i < n; i++) {
        v[i] = A[i][j];
    }
    
    double Lv = norma(v);
    if (Lv < 1e-10) return identidade(n);
    
    v_linha[j + 1] = (A[j+1][j] >= 0) ? Lv : -Lv;
    
    Vector N(n);
    for (int i = 0; i < n; i++) {
        N[i] = v[i] - v_linha[i];
    }
    
    double normaN = norma(N);
    if (normaN < 1e-10) return identidade(n);
    
    Vector n_normalizado = normalize(N);
    
    Matrix H = identidade(n);
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            H[i][k] -= 2.0 * n_normalizado[i] * n_normalizado[k];
        }
    }
    
    return H;
}

// Transformacao completa de Householder
pair<Matrix, Matrix> householderTridiagonal(const Matrix& A, int n) {
    Matrix H_acumulada = identidade(n);
    Matrix T = A;
    
    cout << "~~ Transformacoes de Householder ~~" << endl;
    
    for (int j = 0; j < n - 2; j++) {
        cout << "\nPasso " << j+1 << ":" << endl;
        
        Matrix Hj = hausholder(T, j, n);
        
        // T = Hj * T * Hj (transformacao de similaridade)
        Matrix temp = multiplicarMatrizes(Hj, T);
        T = multiplicarMatrizes(temp, Hj);
        
        // Acumular transformacoes
        H_acumulada = multiplicarMatrizes(H_acumulada, Hj);
        
        cout << "Matriz T apos passo " << j+1 << ":" << endl;
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                cout << fixed << setprecision(4) << T[i][k] << "\t";
            }
            cout << endl;
        }
    }
    
    return {H_acumulada, T};
}

// QR do meu arquivo .cpp. tive que trazer pra ca
double somaTermosQuadradosAbaixoDiagonal(const Matrix& A) {
    int n = A.size();
    double soma = 0.0;
    for (int i = 1; i < n; ++i)
        for (int j = 0; j < i; ++j)
            soma += A[i][j] * A[i][j];
    return soma;
}

Matrix matrizGivens(const Matrix& A, int i, int j, int n) {
    Matrix G = identidade(n);
    double epsilon = 1e-10;
    
    if (fabs(A[i][j]) <= epsilon) return G;

    double r = sqrt(A[j][j] * A[j][j] + A[i][j] * A[i][j]);
    double c = A[j][j] / r;
    double s = -A[i][j] / r;
    
    G[j][j] = c;
    G[i][i] = c;
    G[j][i] = s;
    G[i][j] = -s;
    
    return G;
}

pair<Matrix, Matrix> decomposicaoQR(const Matrix& A, int n) {
    Matrix Q = identidade(n);
    Matrix R = A;

    for (int j = 0; j < n - 1; ++j) {
        for (int i = j + 1; i < n; ++i) {
            if (fabs(R[i][j]) > 1e-10) {
                Matrix G = matrizGivens(R, i, j, n);
                R = multiplicarMatrizes(G, R);
                Q = multiplicarMatrizes(Q, transpor(G));
            }
        }
    }
    return {Q, R};
}

pair<Matrix, Vector> metodoQR(const Matrix& T, int n, double epsilon) {
    Matrix X = identidade(n);
    Matrix Ak = T;
    Vector autovalores(n);
    
    cout << "\n~~ Metodo QR ~~" << endl;
    
    for (int iter = 0; iter < 1000; iter++) {
        auto [Q, R] = decomposicaoQR(Ak, n);
        Ak = multiplicarMatrizes(R, Q);
        X = multiplicarMatrizes(X, Q);

        double convergencia = somaTermosQuadradosAbaixoDiagonal(Ak);
        
        if (iter % 100 == 0) {
            cout << "Iteracao " << iter << ": convergencia = " << scientific << setprecision(4) << convergencia << endl;
        }
        
        if (convergencia < epsilon) {
            cout << "QR convergiu em " << iter << " iteracoes!" << endl;
            break;
        }

        // Deflacao
        for (int i = 1; i < n; i++)
            for (int j = 0; j < i; j++)
                if (fabs(Ak[i][j]) < 1e-8) Ak[i][j] = 0.0;
    }

    for (int i = 0; i < n; i++) {
        autovalores[i] = Ak[i][i];
    }

    return {X, autovalores};
}

int main() {
    // Questao 2: Problema generalizado Ax = lambda*Bx
    Matrix B = {
        {2, 0, 0, 0, 0},
        {0, 3, 0, 0, 0},
        {0, 0, 4, 0, 0},
        {0, 0, 0, 5, 0},
        {0, 0, 0, 0, 6}
    };
    
    Matrix A = {
        {8, 4, 2, 6, 10},
        {6, 27, 9, 21, 24},
        {4, 12, 44, 24, 16},
        {15, 35, 30, 70, 35},
        {30, 48, 24, 42, 60}
    };
    
    int n = A.size();
    double epsilon = 1e-8;
    
    cout << "~~ QUESTAO 3: HOUSEHOLDER + QR ~~" << endl;
    
    // a) Transformacao para forma padrao: C = B^(-1) * A
    cout << "\na) Transformacao para forma padrao:" << endl;
    Matrix Binv = inversa(B);
    Matrix C = multiplicarMatrizes(Binv, A);
    
    cout << "Matriz C = B^(-1) * A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(4) << C[i][j] << "\t";
        }
        cout << endl;
    }
    
    // b) Metodo de Householder
    cout << "\nb) Metodo de Householder:" << endl;
    auto [H, T] = householderTridiagonal(C, n);
    
    cout << "\nMatriz H acumulada:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(4) << H[i][j] << "\t";
        }
        cout << endl;
    }
    
    cout << "\nMatriz T final (tridiagonal):" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(4) << T[i][j] << "\t";
        }
        cout << endl;
    }
    
    // c) Metodo QR
    cout << "\nc) Metodo QR:" << endl;
    auto [X, Lambda] = metodoQR(T, n, epsilon);
    
    cout << "\nMatriz X (autovetores de T):" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(4) << X[i][j] << "\t";
        }
        cout << endl;
    }
    
    cout << "\nMatriz Lambda (autovalores):" << endl;
    Matrix Diag = identidade(n);
    for (int i = 0; i < n; i++) {
        Diag[i][i] = Lambda[i];
        cout << "lambda" << i+1 << " = " << fixed << setprecision(6) << Lambda[i] << endl;
    }
    
    // d) Autovetores de C
    cout << "\nd) Autovetores de C:" << endl;
    Matrix Y = multiplicarMatrizes(H, X);
    
    cout << "Matriz Y = H * X (autovetores de C):" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(4) << Y[i][j] << "\t";
        }
        cout << endl;
    }
    
    // e) Verificacao da decomposicao espectral C = Y * Lambda * Y^T
    cout << "\ne) Verificacao da decomposicao espectral:" << endl;
    Matrix YDiag = multiplicarMatrizes(Y, Diag);
    Matrix YDiagYT = multiplicarMatrizes(YDiag, transpor(Y));
    
    cout << "C reconstruida = Y * Lambda * Y^T:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << fixed << setprecision(4) << YDiagYT[i][j] << "\t";
        }
        cout << endl;
    }
    
    // Calcular erro
    double erro = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double diff = C[i][j] - YDiagYT[i][j];
            erro += diff * diff;
        }
    }
    erro = sqrt(erro);
    
    cout << "\nErro ||C - Y*Lambda*Y^T||_F = " << scientific << setprecision(4) << erro << endl;
    
    // f) Verificacao do problema generalizado
    cout << "\nf) Verificacao A*x = lambda*B*x para o autovalor dominante:" << endl;
    int idx_dom = 0;
    for (int i = 1; i < n; i++) {
        if (fabs(Lambda[i]) > fabs(Lambda[idx_dom])) {
            idx_dom = i;
        }
    }
    
    Vector x_dom(n);
    for (int i = 0; i < n; i++) {
        x_dom[i] = Y[i][idx_dom];
    }
    
    Vector Ax = multiplicarMatrizVetor(A, x_dom);
    Vector Bx = multiplicarMatrizVetor(B, x_dom);
    Vector lambdaBx(n);
    for (int i = 0; i < n; i++) {
        lambdaBx[i] = Lambda[idx_dom] * Bx[i];
    }
    
    cout << "Autovalor dominante: lambda = " << Lambda[idx_dom] << endl;
    cout << "A*x = ";
    for (double val : Ax) cout << fixed << setprecision(4) << val << " ";
    cout << "\nlambda*B*x = ";
    for (double val : lambdaBx) cout << fixed << setprecision(4) << val << " ";
    
    double erro_gen = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = Ax[i] - lambdaBx[i];
        erro_gen += diff * diff;
    }
    erro_gen = sqrt(erro_gen);
    
    cout << "\nErro ||A*x - lambda*B*x|| = " << scientific << setprecision(4) << erro_gen << endl;
    
    return 0;
}
