#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "../utilitarios.h" 
//pra nao ter que escrever a funcao identidade, produto escalar, produto matrizes pela trigesima vez
//vo ficar usando de novo mesmo, mais facil jogar num header e so ficar pedindo dentro do cpp mesmo. mais facil incluir metodos/funcoes novas

using namespace std;
//infelizmente nomes cumpridos de funcao me ajudam a lembrar doq ela deve fazer, feio mas funciona
double somaTermosQuadradosAbaixoDiagonal(const Matrix& A) {
    int n = A.size();
    double soma = 0.0;
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            soma += A[i][j] * A[i][j];
        }
    }
    return soma;
}
//peguei o nome logo todo do pseudocodigo pra eu nao me perder
Matrix matrizJacobiBaseadaNoElemento_ij_DeRvelha(const Matrix& A, int i, int j, int n) {
    Matrix J = identidade(n);
    double epsilon = 1e-10;
    
    if (fabs(A[i][j]) <= epsilon) return J;

    // rotacao de givens pra zerar A[i][j]
    double r = sqrt(A[j][j] * A[j][j] + A[i][j] * A[i][j]);
    double cos = A[j][j] / r;  // coseno
    double sen = -A[i][j] / r; // seno
    
    J[j][j] = cos;
    J[i][i] = cos;
    J[j][i] = sen;
    J[i][j] = -sen;
    
    return J;
}
//a tal da decomposicao qr
pair<Matrix, Matrix> decomposicaoQR(const Matrix& A, int n) {
    Matrix Q = identidade(n);
    Matrix R = A;

    for (int j = 0; j < n - 1; ++j) {
        for (int i = j + 1; i < n; ++i) {
            if (fabs(R[i][j]) > 1e-10) {
                Matrix G = matrizJacobiBaseadaNoElemento_ij_DeRvelha(R, i, j, n);
                R = multiplicarMatrizes(G, R);
                Q = multiplicarMatrizes(Q, transpor(G));
            }
        }
    }
    return {Q, R};
}

// Qr agora com deflacao
pair<Matrix, Vector> metodoQRcode(const Matrix& A, int n, double epsilon) {
    Matrix P = identidade(n);
    Matrix Ak = A;
    Vector autovalores(n);
    
    cout << "~ Iniciando QR ~" << endl;
    
    // QR com criterio pesado de convergencia
    for (int iter = 0; iter < 500; iter++) {
        auto [Q, R] = decomposicaoQR(Ak, n);
        Ak = multiplicarMatrizes(R, Q);
        P = multiplicarMatrizes(P, Q);
        
        double convergencia = somaTermosQuadradosAbaixoDiagonal(Ak);
        
        if (iter % 50 == 0) {
            cout << "Iteração " << iter << ": convergência = " << scientific << setprecision(4) << convergencia << endl;
        }
        //Zerando os elementos da diagonal inferior, tive que pesquisar oq era a deflaca pra esse troco funcionar
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (fabs(Ak[i][j]) < 1e-8) {
                    Ak[i][j] = 0.0;
                }
            }
        }
        
        if (convergencia < epsilon) {
            cout << "Convergiu em " << iter << " iterações :D (tomara nao papoque meu pc)" << endl;
            break;
        }
    }
    
    // Coleta os autovalores da diagonal principal
    for (int i = 0; i < n; i++) {
        autovalores[i] = Ak[i][i];
    }
    
    //Bota os autovalores em ordem decrescente pra ficar organizado
    for (int i = 0; i < n-1; i++) {
        for (int j = i+1; j < n; j++) {
            if (fabs(autovalores[i]) < fabs(autovalores[j])) {
                swap(autovalores[i], autovalores[j]);
                // Trocar também as colunas de P
                for (int k = 0; k < n; k++) {
                    swap(P[k][i], P[k][j]);
                }
            }
        }
    }
    
    return {P, autovalores};
}

int main() { //tirei logo os valores de λ atraves de um calculador de matrizes online, pra conferir com os resultados
    Matrix A = {
        {5.0, 2.0, 1.0},
        {2.0, 3.0, 1.0},
        {1.0, 1.0, 1.0}
    };//λ1: 0.57637\dots ,\:λ2:1.84653 ,\:λ3:6.57708
    Matrix B = {
        {-2.7083, -2.6824, 0.4543},
        {0.1913, 0.7269, 0.1007},
        {-0.3235, -0.4052, 5.0453}
    };
    Matrix C = {
    {40, 8, 4, 2, 1},
    {8, 30, 12, 6, 2},
    {4, 12, 20, 1, 2},
    {2, 6, 1, 25, 4},
    {1, 2, 2, 4, 5}
    };//λ\approx \:4.01488\dots ,\:λ\approx \:11.64242\dots ,\:λ\approx \:23.64807\dots ,\:λ\approx \:31.31146\dots ,\:λ\approx \:49.38314\dots 
    Vector v0 = {1, 1, 1}; // Vetor inicial pra usar em A e B
    Vector v1 = {1, 1, 1, 1, 1}; // Vetor inicial para C
    double epsilon = 1e-6;
    int n = A.size();
    int m = B.size();
    int p = C.size();

    cout << "\n ~ Testando QR ~ " << endl;
    auto [P, autovalores] = metodoQRcode(A, n, epsilon);

    cout << fixed << setprecision(4);
    cout << "\n~ Resultados do Metodo QR  ~" << endl;
    cout << "Autovalores encontrados: ";
    for (double val : autovalores) cout << val << " ";
    cout << "\nAutovalores esperados: 6.5771 1.8465 0.5764 (calculei pelo matrixcalc.org)" << endl;
    
    cout << "\nAutovetores (colunas de P):\n";
    for (int i = 0; i < n; i++) {
        cout << "v" << i+1 << ": ";
        for (int j = 0; j < n; j++) {
            cout << P[j][i] << " ";
        }
        cout << "\n";
    }
    //Confere A*v = λ*v nos autovetores, pra ver se os valores batem
    cout << "\n~ Verificação A*v = λ*v ~" << endl;
    for (int i = 0; i < n; i++) {
        Vector vi(n);
        for (int j = 0; j < n; j++) {
            vi[j] = P[j][i];
        }
        Vector Avi = multiplicarMatrizVetor(A, vi);
        
        cout << "Autovalor " << i+1 << " (λ=" << autovalores[i] << "):" << endl;
        cout << "  A*v: ";
        for (double val : Avi) cout << val << " ";
        cout << "\n  λ*v: ";
        for (int j = 0; j < n; j++) {
            cout << autovalores[i] * vi[j] << " ";
        }
        
        // Calcula o erro 
        double erro = 0.0;
        for (int j = 0; j < n; j++) {
            double diff = Avi[j] - autovalores[i] * vi[j];
            erro += diff * diff;
        }
        cout << "\n  Erro é ||A*v - λ*v||: " << sqrt(erro) << endl;
    }
    
    /*// Teste de comparacao pq nao tava dando certo
    cout << "\n~ Comparação com Método da Potência ~" << endl;
    Vector v_inicial = {1, 1, 1};
    auto [lambda_dom, v_dom] = Regularparainversa(A, v_inicial, epsilon);
    cout << "Autovalor dominante (Potência): " << lambda_dom << endl;
    
    // Verificar autovalor menor usando potência inversa
    Matrix A_inv = inversa(A);
    auto [lambda_inv, v_inv] = Regularparainversa(A_inv, v_inicial, epsilon);
    double lambda_menor = 1.0 / lambda_inv;
    cout << "Autovalor menor (Potência Inversa): " << lambda_menor << endl;*/

    
    return 0;
}
