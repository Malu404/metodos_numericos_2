#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include "../../utilitarios.h"
using namespace std;
// Função para calcular o método de deslocamento
void deslocamento(const Matrix& A, Vector& v0, double epsilon, double deslocamento){
    int n = A.size();
    matriz I = identidade(n);
    A = A - I*deslocamento;
}
