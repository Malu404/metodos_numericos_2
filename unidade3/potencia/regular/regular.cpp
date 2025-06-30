#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "../../utilitarios.h"//as funcoes comuns pertinentes a unidade 3 estao aqui

using namespace std;

// Função para calcular o método de potência regular
void Regular(const Matrix& A, Vector& v0, double epsilon){
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
    
    cout << fixed << setprecision(6);
    cout << "Autovalor dominante: " << lambda_new << endl;
    cout << "Autovetor dominante: ";
    for (double val : vk_new) {
        cout << val << " ";
    }
    cout << endl;
}
int main() {
    Matrix A = {
        {5.0, 2.0, 1.0},
        {2.0, 3.0, 1.0},
        {1.0, 1.0, 1.0}
    };

    Vector v0 = {1, 1, 1};
    double epsilon = 1e-6;

    Regular(A, v0, epsilon);
    return 0;
}

    