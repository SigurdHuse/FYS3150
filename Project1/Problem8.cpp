#include <bits/stdc++.h>
using namespace std; 

const long double E = 2.71828182845904523536;
const long double val = pow(E, -10);

long double u(long double x){
    return 1.0 - (1.0 - val)*x - pow(E, -10*x); 
}

int main()
{
    for(double n = 10; n <= 1e4; n *= 10){
        vector<long double> xvalues(n+1), yvalues(n+1);
        long double step = 1.0/n;

        for(double i = 0; i <= n; ++i){
            xvalues[i] = step*i;
            yvalues[i] = u(xvalues[i]);
        }
        ofstream outfile;
        outfile.open("exact_solution_n_" + to_string((int)n) + ".txt");

        int width = 20, precision = 12;

        for(int i = 0; i <= n; ++i){
            outfile << setw(width) << setprecision(precision) << scientific << xvalues[i] << ", "
                    << setw(width) << setprecision(precision) << scientific << yvalues[i] << endl;
        }
    }
}
