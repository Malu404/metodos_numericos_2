#ifndef UTILITARIOS_H
#define UTILITARIOS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <utility>
using namespace std;

// o tipo pra matriz e vetor,ai so chama na funcao que tiver precisando
typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

//  diretamente dos confins de cg, normalize.
double norma(const Vector& v) {
    double sum = 0.0;
    for (double val : v) {
        sum += val * val;
    }
    return sqrt(sum);
}

// diretamente dos confins de cg, normalize.
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

// gera uma identidade do tamanho que eu precisar(o n na chamada da funcao)
Matrix identidade(int n) {
    Matrix I(n, Vector(n, 0.0));
    for (int i = 0; i < n; i++) {
        I[i][i] = 1.0;
    }
    return I;
}

// matrix x matrix, usa especialmente quando comeca com uma identidade e depois altera ela
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

// matriz x vetor, precisa disso na regular e na inversa
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

// transposta da matriz
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
//vetor x vetor
double produtoescalar(const Vector& a, const Vector& b){
    double resultado = 0.0;
    for (size_t i = 0; i<a.size(); i++){
        resultado += a[i] * b[i];
    }
    return resultado;
}

// Precisa inverter na inversa(ta no nome ne) pra poder achar o menor autovalor, trouxe pra ca pra nao dar bloat no codigo
Matrix inversa(const Matrix& A) {
    int n = A.size();
    Matrix I = identidade(n);
    Matrix B = A;

    for (int i = 0; i < n; ++i) {
        double diag = B[i][i];
        if (fabs(diag) < 1e-10) {
            cerr << "Erro: pivo zero na inversao" << endl;
            exit(1);
        }

        // Normaliza a linha
        for (int j = 0; j < n; ++j) {
            B[i][j] /= diag;
            I[i][j] /= diag;
        }

        // Limpa outros elementos da coluna
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
//Matriz - Matriz, com um ifzinho so pra barrar de fazer uma operacao q nao pode fazer
Matrix subtrairmatrizes(const Matrix& A, const Matrix& B) {
    int n = A.size();
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        cerr << "Erro: Matrizes de tamanhos diferentes nao podem ser subtraidas" << endl;
        exit(1);
    }
    
    Matrix resultado(n, Vector(A[0].size(), 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < A[0].size(); ++j) {
            resultado[i][j] = A[i][j] - B[i][j];
        }
    }
    return resultado;
}

// Matriz x escalar, precisa pra deslocamento (achar quem ta perto de 4, 10, 20 etc etce tc)
Matrix multiplicarMatrizPorEscalar(const Matrix& A, double escalar) {
    int n = A.size();
    int m = A[0].size();
    Matrix resultado(n, Vector(m, 0.0));
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            resultado[i][j] = A[i][j] * escalar;
        }
    }
    return resultado;
}

// tive que colocar ela aqui pra poder usar ela na inversa e no deslocamento, ia ser so regular mas ai deu BO com a regular.cpp, entao virou regular pra inversa
pair<double, Vector> Regularparainversa(const Matrix& A, const Vector& v0, double epsilon) {
    Vector vk_new = normalize(v0);
    Vector vk_old = v0;
    double lambda_new = 0.0;
    double lambda_old = 0.0;
    Vector x1;
    int iteracoes = 0;
    
    do {
        vk_old = vk_new;
        x1 = multiplicarMatrizVetor(A, vk_old);
        lambda_old = lambda_new;
        lambda_new = produtoescalar(x1, vk_old);
        vk_new = normalize(x1);
        iteracoes++;
        
        if (iteracoes > 10000) {
            cerr << "Aviso: Maximo de iteracoes atingido na potencia!" << endl;
            break;
        }
    } while (fabs(lambda_new - lambda_old) > epsilon);
    
    return {lambda_new, vk_new};
}

#endif // UTILITARIOS_H
