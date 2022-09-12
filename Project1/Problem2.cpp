#include <bits/stdc++.h>
using namespace std; 


const int len = 100000;
const double E = 2.71828182845904523536;
const double val = pow(E, -10);

double f(double x){
    return 1.0 - (1.0 - val)*x - pow(E, -10*x); 
}


int main()
{
	vector<double> xvalues(len + 1), yvalues(len+1);
    double step = 1.0/len;

    for(double i = 0; i < len; ++i){
        xvalues[i] = step*i;
        yvalues[i] = f(xvalues[i]);
    }
    xvalues[len] = 1.0;
    yvalues[len] = f(1.0);
    
    ofstream outfile;
    outfile.open("exact_solution.txt");

    int width = 12, precision = 4;

    for(int i = 0; i <= len; ++i){
        outfile << setw(width) << setprecision(precision) << scientific << xvalues[i]
                << setw(width) << setprecision(precision) << scientific << yvalues[i] << endl;
    }

}
