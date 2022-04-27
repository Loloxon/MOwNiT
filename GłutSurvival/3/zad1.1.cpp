#include <iostream>
#include <cmath>



using namespace std;


double secant1(double (*f)(double), double x0, double x1, double eps, int *iter){
    (*iter) = 1;
    double tmp;
    while(1){
        if(abs(x1-x0) < eps) break;
        tmp = x1 - (f(x1)*(x1-x0)/(f(x1)-f(x0)));
        x0 = x1;
        x1 = tmp;
        (*iter)++;
    }
    return x1;

}

double secant2(double (*f)(double), double x0, double x1, double eps, int *iter){
    (*iter) = 1;
    double tmp;
    while(1){
        if(abs(f(x0)) < eps) break;
        tmp = x1 - (f(x1)*(x1-x0)/(f(x1)-f(x0)));
        x0 = x1;
        x1 = tmp;
        (*iter)++;
    }
    return x1;

}

double f(double x){
    return pow(x,10) - pow(1-x,15);
}



int main(int argc, char* argv[]){

    int iter;
    double x0 = 0.1;
    double x1 = 2.1;
    double eps = stod(argv[1]);
    double result;
    int mode = atoi(argv[2]);
    int dir = atoi(argv[3]);
    if(dir == 1){

        while (x0 < 2.1){
            if(mode == 1) result = secant1(f,x0,x1,eps,&iter);
            else result = secant2(f,x0,x1,eps,&iter);
            cout << x0 << " " << x1 << " " << eps << " " << iter << " " << result << endl;
            x0+=0.1;
        }
    } else {
        while (x1 > 0.1){
            if(mode == 1) result = secant1(f,x0,x1,eps,&iter);
            else result = secant2(f,x0,x1,eps,&iter);
            cout << x0 << " " << x1 << " " << eps << " " << iter << " " << result << endl;
            x1-=0.1;
        }
    }

/*
 result = secant1(f,x0,x1,eps,&iter);
        cout << "Parameters: x0 = "<< x0 << "\t x1 = " << x1 << "\t eps = "<< eps << endl;
        cout << "The first result is: " << result << " in " << iter << " iterations." << endl;
        result = secant2(f,x0,x1,eps,&iter);
        cout << "The second result is: " << result << " in " << iter << " iterations." << endl;
        cout << endl<< endl;
        }




        result = secant1(f,x1,x0,eps,&iter);
        cout << "Parameters: x0 = "<< x1 << "\t x1 = " << x0 << "\t eps = "<< eps << endl;
        cout << "The first result is: " << result << " in " << iter << " iterations." << endl;
        result = secant2(f,x1,x0,eps,&iter);
        cout << "The second result is: " << result << " in " << iter << " iterations." << endl;
        cout << endl<< endl;
*/


    return 0;
}