#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include "../../utilitarios.h"

using namespace std;

// funcao do deslocamento
pair<double, Vector> deslocamento(const Matrix& A, const Vector& v0, double epsilon, double shift) {
    int n = A.size();
    Matrix I = identidade(n);
    
    // calculo minha matriz depois de contar com o shift do deslocamento
    Matrix Aa = subtrairmatrizes(A, multiplicarMatrizPorEscalar(I, shift));
    
    // depois aplico a potencia inversa em Aa
    Matrix Ainv = inversa(Aa);
    auto [lambda_dom_inversa, v_dom_inversa] = Regularparainversa(Ainv, v0, epsilon);
    
    // como eu inverti a matriz, desinverto o autovalor e adiciono o shift
    double lambda_deslocamento = 1.0 / lambda_dom_inversa + shift;
    Vector vec_deslocamento = v_dom_inversa;
    
    // print
    cout << fixed << setprecision(6);
    cout << "Autovalor mais proximo de " << shift << " (metodo do deslocamento): " << lambda_deslocamento << endl;
    cout << "Autovetor correspondente: ";
    for (double val : vec_deslocamento) {
        cout << val << " ";
    }
    cout << endl;
    
    return {lambda_deslocamento, vec_deslocamento};
}

int main() {//valores que eu calculei pelo matrixcalc
   Matrix B = {
        {2, 0, 0, 0, 0},
        {0, 3, 0, 0, 0},
        {0, 0, 4, 0, 0},
        {0, 0, 0, 5, 0},
        {0, 0, 0, 0, 6}
    };

    // Matriz A do problema
    Matrix prova = {
        {8, 4, 2, 6, 10},
        {6, 27, 9, 21, 24},
        {4, 12, 44, 24, 16},
        {15, 35, 30, 70, 35},
        {30, 48, 24, 42, 60}
    };

    Vector v1 = {1, 1, 1, 1, 1};
    double epsilon = 1e-6;
    double shift = 20.0; // valor alvo aproximado do autovalor

    // Transformacao do problema generalizado: C = B^(-1) * A
    Matrix Binv = inversa(B);
    Matrix C = multiplicarMatrizes(Binv, prova);

    cout << "~~ Metodo da Potencia com Deslocamento ~~" << endl;
    auto [lambda_dominante, vetor_dominante] = deslocamento(C, v1, epsilon, shift);

    // ~~~~~~~~~~~~~~~~~~~~=
    // Verificacao: A * x ~ lambda * B * x
    // ~~~~~~~~~~~~~~~~~~~~=
    cout << "\n~~ Verificacao: A * x ~ lambda * B * x ~~" << endl;

    Vector Ax = multiplicarMatrizVetor(prova, vetor_dominante);
    Vector Bx = multiplicarMatrizVetor(B, vetor_dominante);

    cout << "A*x:          ";
    for (double val : Ax) {
        cout << val << " ";
    }
    cout << "\nlambda*B*x:   ";
    for (double val : Bx) {
        cout << lambda_dominante * val << " ";
    }
    cout << endl;

    return 0;
}