#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <armadillo>

using namespace std;
using namespace arma;

vec SOR2(double omega, Mat<double> A, vec B, vec start, double ro, int *iter, int size){
        
    
        vec X = start;
        vec tmp(size);
        int it=1;
        while(1){
            for(int i=0; i<size; i++){
                tmp(i) =  (1- omega) * X(i);
                double b = 0.0;
                for(int j=0; j<i; j++){
                    b += A(i,j)*tmp(j);
                }
                for(int j=i+1; j<size; j++){
                    b += A(i,j)*X(j);
                }
                b = B(i) -b;
                b = b*omega/A(i,i);
                tmp(i) = tmp(i) + b;
            }
            if(max(abs(A*tmp - B)) < ro ){
                for(int i=0; i<size; i++){
                    X(i) = tmp(i);
                }
                (*iter) = it;
                break;
            }
            for(int i=0; i<size; i++){
                X(i) = tmp(i);
            }
            it++;
        }
        return X;
}   

vec SOR1(double omega, Mat<double> A, vec B, vec start, double ro, int *iter, int size){
        
        
        vec X = start;
        vec tmp(size);
        int it=1;
        while(1){
            for(int i=0; i<size; i++){
                tmp(i) =  (1- omega) * X(i);
                double b = 0.0;
                for(int j=0; j<i; j++){
                    b += A(i,j)*tmp(j);
                }
                for(int j=i+1; j<size; j++){
                    b += A(i,j)*X(j);
                }
                b = B(i) -b;
                b = b*omega/A(i,i);
                tmp(i) = tmp(i) + b;
            }
            if(max(abs(tmp-X)) < ro ){
                for(int i=0; i<size; i++){
                    X(i) = tmp(i);
                }
                (*iter) = it;
                break;
            }
            for(int i=0; i<size; i++){
                X(i) = tmp(i);
            }
            it++;
        }
        return X;
}   

double spec_radius(Mat<double> A, int size){
    vec N(size);
    for(int i=0; i<size; i++){
        N(i) = 1/A(i,i);
    }
    Mat<double> M(size,size);
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            if(j==i) M(i,j) = 0.0;
            else M(i,j) = - (A(i,j) * N(i));
        }
        }
    cx_vec eigens = eig_gen(M);
    double result = max(abs(eigens));

    return result;
}

Mat<double> Matrix_A(int size){
    Mat<double> A(size,size);
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            if(i==j) A(i,j) = 7.0;
            else A(i,j) = 1.0 /(size - i - 1 - j - 1 +0.5);
        }
    }

    return A;   
}

vec Matrix_X(int size){
    vec X(size);
    for(int i=0; i<size; i++){
        if(i%2==0) X(i) = 1.0;
        else X(i) = -1.0;
    }

    return X;
}


int main(int argc, char* argv[]){
    srand(time(0));
    int size = atoi(argv[1]);
    double ro = stod(argv[2]);
    int start_vec = atoi(argv[3]);
    double omega = stod(argv[4]);
    vec st_vec;
    vec X = Matrix_X(size);
    switch(start_vec){
        case 1: st_vec = zeros<vec>(size);
            break;
        case -1: st_vec = X;
            break;
        case 0: st_vec = X*-1;
            break;
        case 100: st_vec = X*100;
            break;
        case -100: st_vec = X*-100;
            break;

    }

    Mat<double> A = Matrix_A(size);
    vec B(size);
    B = A*X;
    vec X_calc_1;
    vec X_calc_2;
    int iterations1;
    int iterations2;
    cout.precision(10);
    clock_t start = clock();
    X_calc_1 = SOR1(omega, A, B , st_vec, ro,  &iterations1, size);
    clock_t time1 = clock() - start;
    start = clock();
    X_calc_2 = SOR2(omega, A, B, st_vec, ro, &iterations2, size);
    clock_t time2 = clock() - start;
    /*cout << "Size: " << size << "\tRo: " << ro << "\tomega: " << omega << endl;
    cout << "Starting vector is: ";
    switch(start_vec){
        case 1: cout << "same as the result" << endl;
            break;
        case -1: cout << "result vector multiplied by -1" << endl;
            break;
        case 0: cout << "only zeros" << endl;
            break;
        case 100: cout << "result vector mulitplied by 100" << endl;
            break;
        case -100: cout << "result vector mulitplied by -100" << endl;
            break;
    }
    
    cout << "First type of ending condition: " << endl;
    cout <<endl <<"The error is: " << max(abs(X_calc_1-X)) << endl;
    cout << "Number of iterations was " << iterations1 << endl;

    cout << endl << "Second type of ending condtition: "<< endl;
    cout <<endl <<"The error is: " << max(abs(X_calc_2-X)) << endl;
    cout << "Number of iterations was " << iterations2 << endl;

    cout << "Spectral Radius of iteration matrix is: "<< spec_radius(A,size) << endl;
    cout << "---------------------------------------------------------------" << endl;
    */
    cout << size << " " << ro << " " << start_vec << " ";
    cout << omega << " " << spec_radius(A,size) << " ";
    cout << max(abs(X_calc_1-X)) << " " << iterations1 << " " << time1/((double)CLOCKS_PER_SEC/1000) << " ";
    cout << max(abs(X_calc_2-X)) << " " << iterations2 << " " << time2/((double)CLOCKS_PER_SEC/1000) << " ";
    cout << endl;
    return 0;
}