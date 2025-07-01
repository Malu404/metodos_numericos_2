#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "../../utilitarios.h"//as funcoes comuns pertinentes a unidade 3 estao aqui

using namespace std;

// calcula os autovalor e autovetor pelo metodo da potencia regular
pair<double, Vector> Regular(const Matrix& A, const Vector& v0, double epsilon){
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
        
        // contador de iteracoes
        if (iteracoes % 100 == 0) {
            cout << "Iteracao " << iteracoes << ": λ = " << lambda_new << ", diferença = " << fabs(lambda_new - lambda_old) << endl;
        }
        
        if (iteracoes > 10000) {
            cout << "Aviso: Maximo de iteracoes atingido!" << endl;
            break;
        }
    } while (fabs(lambda_new - lambda_old) > epsilon);
    
    cout << fixed << setprecision(6);
    cout << "Autovalor dominante: " << lambda_new << " (apos " << iteracoes << " iteracoes)" << endl;
    cout << "Autovetor dominante: ";
    for (double val : vk_new) {
        cout << val << " ";
    }
    cout << endl;
    
    return {lambda_new, vk_new};
}
int main() {//valores que eu calculei pelo matrixcalc
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

    cout << "=== Testando matriz B ===" << endl;
    auto [lambdab, vetorb] = Regular(B, v0, epsilon);
    cout << "=== Testando matriz C ===" << endl;
    auto [lambdac, vetorc] = Regular(C, v1, epsilon);
    // Verificação: A*v deve ser aproximadamente λ*v
    Vector Av = multiplicarMatrizVetor(C, vetorc);
    cout << "\nVerificação (A*v): ";
    for (double val : Av) {
        cout << val << " ";
    }
    cout << "\nλ*v: ";
    for (double val : vetorc) {
        cout << lambdac * val << " ";
    }
    cout << endl;
    
    return 0;
}

    