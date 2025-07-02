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
            cout << "Iteracao " << iteracoes << ": lambda = " << lambda_new << ", diferenca = " << fabs(lambda_new - lambda_old) << endl;
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

    Vector v1 = {1, 1, 1, 1, 1};
    double epsilon = 1e-6;
    Matrix Binv = inversa(B);
    Matrix C = multiplicarMatrizes(Binv, A);
    //cout << "~~ Testando matriz B ~~" << endl;
    //auto [lambdab, vetorb] = Regular(Binv, v1, epsilon);
    cout << "~~ Testando matriz C ~~" << endl;
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
    