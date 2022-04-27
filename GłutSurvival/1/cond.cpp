#include <iostream>
#include <cmath>
#include <vector>
#include <armadillo>

using namespace std;

void print(vector< vector<float> > A, int size) {
    int n = size;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            cout << A[i][j] << "\t";
            if (j == n-1) {
                cout << "| ";
            }
        }
        cout << "\n";
    }
    cout << endl;
}

vector<vector<float>> matrix_one_f(int size) {
    vector<float> row(size);
    vector<vector <float>> A(size,row); 
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            if(i==0) A[i][j] = 1.0;
            else A[i][j] = 1.0/(i+1+j+1-1);
            
        } 
    }
    return A;

}


vector<vector<float>> matrix_two_f(int size){
    vector<float> row(size);
    vector<vector <float>> A(size,row);
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            if(i<=j) A[i][j] = 2.0*(i+1)/(j+1);
            else A[i][j] = A[j][i];
        }
    }

    return A;
}


int main(int argc,  char** argv) {
    int size = atoi(argv[1]);
    cout.precision(5);  
    cout << "Calculating conditions for matrices of size:" << size << endl;
    vector<vector<float>> A = matrix_one_f(size);
    vector<vector<float>> B = matrix_two_f(size);
    arma::Mat<float> A_mat(size,size);
    arma::Mat<float> B_mat(size,size);
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            A_mat(i,j) = A[i][j];
            B_mat(i,j) = B[i][j];
        }
    }
    float cond_A = arma::cond(A_mat);
    float cond_B = arma::cond(B_mat);
    cout << "Condition number of Matrix A is: "<< cond_A << endl;
    cout << "Condition number of Matrix B is:"<< cond_B << endl;
    
    return 0;
}