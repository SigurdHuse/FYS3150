#include <bits/stdc++.h>
using namespace std;
#include <chrono>
const double E = 2.71828182845904523536;
const double val = pow(E, -10);

double f(double x){
    return 1.0 - (1 - val)*x - pow(E, -10*x); 
}


vector<double> solve_matrix_specialized(int n, vector<double> g){
    vector<double> v(n+1);

    const double w = - 0.5;
    const double divisor = (2.0 + w);
    for(int i = 1; i < n; ++i){
        g[i] = g[i] - w * g[i-1];
    }

    v[n-1] = g[n-1] / divisor;

    for(int i = n-2; i > 0; --i){
        v[i] = (g[i] + v[i+1]) / divisor;
    }
    return v;
}

vector<double> solve_matrix(int n, vector<double> a, vector<double> b, vector<double> c, vector<double> g){
    vector<double> v(n+1);

    for(int i = 1; i < n; ++i){
        double w = a[i] / b[i-1];
        b[i] = b[i] - w * c[i-1];
        g[i] = g[i] - w * g[i-1];
    }
    v[n-1] = g[n-1] / b[n-1];

    for(int i = n-2; i > 0; --i){
        v[i] = (g[i] - c[i] * v[i+1]) / b[i];
    }

    return v;
}

vector<double> generate_diag(int n, double value){
    vector<double> v(n, value);
    return v;
}

vector<double> generate_g(int n){
    double h = 1.0/n;
    vector<double> g(n);

    for(double i = 0;i < n; ++i){
        g[i] = -f((i+1)*h + h) + 2*f((i+1)*h) - f((i+1)*h - h);
    }
    return g;
}


int main ()
{   
        ofstream outfile;
        outfile.open("Run_times.txt");
        int width = 12, precision = 6;
        double duration_seconds_1 = 0, duration_seconds_2 = 0;
        
        const int runs = 20;
        for(int n = 10; n <= 1e6; n *= 10){
            for(int j = 0; j < runs; ++j){
                auto g = generate_g(n);
                auto t1_1 = std::chrono::high_resolution_clock::now();
                auto a = generate_diag(n, -1);
                auto b = generate_diag(n, 2);
                auto c = generate_diag(n, -1);
                auto an1 = solve_matrix(n, a, b, c, g);
                auto t2_1 = std::chrono::high_resolution_clock::now();

                duration_seconds_1 += std::chrono::duration<double>(t2_1 - t1_1).count();
                
                auto t1_2 = std::chrono::high_resolution_clock::now();
                auto ans = solve_matrix_specialized(n, g);
                auto t2_2 = std::chrono::high_resolution_clock::now();

                duration_seconds_2 += std::chrono::duration<double>(t2_2 - t1_2).count();
            } 
            outfile << setw(width) << setprecision(precision) << scientific << duration_seconds_1/runs <<", "
                    << setw(width) << setprecision(precision) << scientific << duration_seconds_2/runs << endl;
        }
}