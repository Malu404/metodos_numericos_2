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
    cout << fixed << setprecision(6);
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
    };

    Vector v0 = {0.4, 3, 2};
    double epsilon = 1e-60;

    potenciaInversa(A, v0, epsilon);
    return 0;
}