#include <iostream>
#include <cmath>



using namespace std;


double newton1(double (*f)(double) , double (*df)(double), double x0, double eps, int *iter){
    (*iter) = 0;
    double x1;
    while(1){
        x1 = x0 - f(x0)/df(x0);
        if(abs(x1-x0) < eps) break;
        (*iter)++;
        x0=x1;
    }
    return x1;
}

double newton2(double (*f)(double) , double (*df)(double), double x0, double eps, int *iter){
    (*iter) = 1;
    double x1;
    while(1){
        x1 = x0 - f(x0)/df(x0);
        if(abs(f(x0)) < eps) break;
        (*iter)++;
        x0=x1;
    }
    return x1;
}

double f(double x){
    return pow(x,10) - pow(1-x,15);
}

double df(double x){
    return 10*pow(x,9) + 15*pow(1-x,14);
}


int main(int argc, char* argv[]){

    int iter;
    double x0 = 0.1;
    double eps = stod(argv[1]);
    double result;
    int mode = atoi(argv[2]);
    while (x0 <= 2.15){
        if(mode == 1) result = newton1(f,df,x0,eps,&iter);
        else result = newton2(f,df,x0,eps,&iter);
        cout << x0 << " " << eps << " " << iter << " " << result << endl; 
        x0+=0.1;
    }
    /*
        cout << "Parameters: x0 = "<< x0 << "\t eps = "<< eps << endl;
        cout << "The first result is: " << result << " in " << iter << " iterations." << endl;
        result = newton2(f,df,x0,eps,&iter);
        cout << "The second result is: " << result << " in " << iter << " iterations." << endl;
        cout << endl<< endl;
    */



    return 0;
}