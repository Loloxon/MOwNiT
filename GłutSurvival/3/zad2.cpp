#include <iostream>
#include <cmath>
#include <armadillo>



using namespace std;
using namespace arma;


vec f(vec x){
    vec result(3);
    result(0) = x(0)*x(0) + x(1)*x(1) + x(2) - 1;
    result(1) = 2*x(0)*x(0) + x(1)*x(1) + x(2)*x(2)*x(2) -2;
    result(2) = 3*x(0) - 2*x(1)*x(1)*x(1) - 2*x(2)*x(2) - 3;


    return result;
}

Mat<double> df(vec x){
    Mat<double> result(3,3);
    result(0,0) = 2*x(0);
    result(0,1) = 2*x(1);
    result(0,2) = 1;
    result(1,0) = 4*x(0);
    result(1,1) = 2*x(1);
    result(1,2) = 3*x(2)*x(2);
    result(2,0) = 3;
    result(2,1) = -6*x(1)*x(1);
    result(2,2) = -4*x(2);
    return result;
}


vec newton1(vec x, double eps, vec (*f)(vec), Mat<double> (*df)(vec), int *iter){
    (*iter) = 1; 
    while(1){
        vec x_temp = solve(df(x),(-1)*f(x));
        
        if(max(abs(x_temp)) < eps){
            x=x_temp+x;
            return x;
        }
        x= x_temp + x;
        (*iter)++;  
    
    }
}

vec newton2(vec x, double eps, vec (*f)(vec), Mat<double> (*df)(vec), int *iter){
    (*iter) = 1; 
    while(1){
        vec x_temp = solve(df(x),(-1)*f(x));
        x= x_temp + x;
        if(max(abs(f(x))) < eps){
            return x;
        }
        
        (*iter)++;  
    
    }
}

int main(int argc, char* argv[]){
    double eps = stod(argv[1]);
    vec result1;
    vec result2;
    int iter1;
    int iter2;
    // int mode = atoi(argv[2]);
    arma_rng::set_seed_random();
    vec rand = randu<vec>(3);
    rand = (rand-0.5)*200;
    result1 = newton1(rand, eps, f, df, &iter1);
    result2 = newton2(rand, eps, f, df, &iter2);

    cout << eps << " ["<< rand(0) << "," << rand(1) << "," << rand(2) << "] " << iter1 <<   " ["<< result1(0) << "," << result1(1) << "," << result1(2) << "] ";
    cout << eps << " ["<< rand(0) << "," << rand(1) << "," << rand(2) << "] " << iter2 <<   " ["<< result2(0) << "," << result2(1) << "," << result2(2) << "] " << endl;







}