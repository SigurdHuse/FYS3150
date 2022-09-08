#include <bits/stdc++.h>
using namespace std;

const double E = 2.71828182845904523536;
const double val = pow(E, -10);

double f(double x){
    return 1.0 - (1 - val)*x - pow(E, -10*x); 
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

vector<double> generate_x(int n){
    double h = 1.0/n;
    vector<double> x(n+1);
    for(double i = 0;i <= n; ++i) x[i] = i*h;
    return x;
}


int main()
{  
	for(int n = 10; n <= 1e4; n *= 10){
        vector<double> a, b, c;
        a = generate_diag(n, -1);
        b = generate_diag(n, 2);
        c = generate_diag(n, -1);
        auto g = generate_g(n);
        auto ans = solve_matrix(n, a,b,c, g);
        auto x = generate_x(n);
        ans[0] = 0;
        ans[n] = 0;
        ofstream outfile;
        outfile.open("numeric_solution_n_" + to_string(n) + ".txt");
        int width = 12, precision = 6;

        for(int i = 0; i <= n; ++i){
            outfile << setw(width) << setprecision(precision) << scientific << x[i] <<", "
                    << setw(width) << setprecision(precision) << scientific << ans[i] << endl;
        }
    }
}
