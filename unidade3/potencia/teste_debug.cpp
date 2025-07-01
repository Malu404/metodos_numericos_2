#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "../utilitarios.h"

using namespace std;

// Função melhorada para debug
pair<double, Vector> potenciaDebug(const Matrix& A, const Vector& v0, double epsilon, string nome) {
    cout << "\n=== " << nome << " ===" << endl;
    
    Vector vk = normalize(v0);
    double lambda_old = 0.0;
    double lambda_new = 0.0;
    
    for (int iter = 0; iter < 1000; iter++) {
        Vector x = multiplicarMatrizVetor(A, vk);
        lambda_old = lambda_new;
        lambda_new = produtoescalar(x, vk);
        vk = normalize(x);
        
        if (iter % 100 == 0 || iter < 10) {
            cout << "Iter " << iter << ": λ = " << fixed << setprecision(6) << lambda_new 
                 << ", |Δλ| = " << fabs(lambda_new - lambda_old) << endl;
        }
        
        if (fabs(lambda_new - lambda_old) < epsilon && iter > 5) {
            cout << "Convergiu em " << iter << " iterações!" << endl;
            break;
        }
    }
    
    cout << "Autovalor final: " << lambda_new << endl;
    cout << "Autovetor: ";
    for (double val : vk) {
        cout << val << " ";
    }
    cout << endl;
    
    // Verificação: A*v ≈ λ*v
    Vector Av = multiplicarMatrizVetor(A, vk);
    cout << "Verificação ||A*v - λ*v||: ";
    double erro = 0;
    for (int i = 0; i < vk.size(); i++) {
        double diff = Av[i] - lambda_new * vk[i];
        erro += diff * diff;
    }
    cout << sqrt(erro) << endl;
    
    return {lambda_new, vk};
}

int main() {
    Matrix B = {
        {-2.7083, -2.6824, 0.4543},
        {0.1913, 0.7269, 0.1007},
        {-0.3235, -0.4052, 5.0453}
    };
   
    Matrix C = {
    {40, 8, 4, 2, 1},
    {8, 30, 12, 6, 2},
    {4, 12, 20, 1, 2,},
    {2, 6, 1, 25, 4},
    {1, 2, 2, 4, 5};
    }    
    Vector v0 = {1, 1, 1};
    double epsilon = 1e-6;
    
    cout << "Autovalores esperados: 5.0187, -2.5313, 6.125" << endl;
    
    // Teste método regular
    auto [lambda_dom, v_dom] = potenciaDebug(B, v0, epsilon, "Potência Regular");
    
    // Teste potência inversa  
    cout << "\n=== Calculando Inversa ===" << endl;
    Matrix B_inv = inversa(B);
    auto [lambda_inv, v_inv] = potenciaDebug(B_inv, v0, epsilon, "Potência na Inversa");
    
    double lambda_menor = 1.0 / lambda_inv;
    cout << "\nAutovalor menor (1/λ_inv): " << lambda_menor << endl;
    
    return 0;
}
