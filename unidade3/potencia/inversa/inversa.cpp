#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include "../../utilitarios.h"

using namespace std;

void potenciaInversa(const Matrix& A, const Vector& v0, double epsilon) {
    // Step 1: Inversão da matriz
    Matrix Ainv = inversa(A);

    // Step 2: Aplica método da potência regular sobre A⁻¹
    auto [lambda_dom_inversa, v_dom_inversa] = Regularparainversa(Ainv, v0, epsilon);
    
    // Step 3 e 4: λn = 1 / λdominante, xn = xdominante
    double lambda_n = 1.0 / lambda_dom_inversa;
    Vector xn = v_dom_inversa;

    // Step 5: Imprimir
    cout << fixed << setprecision(4);
    cout << "Autovalor mais proximo de zero (metodo da potencia inversa): " << lambda_n << endl;
    cout << "Autovetor correspondente: ";
    for (double val : xn) {
        cout << val << " ";
    }
    cout << endl;
}
int main() {
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

    potenciaInversa(C, v1, epsilon);
    return 0;
}