#include <iostream>
#include <cmath>
#include <vector>

using namespace std;



void print(vector< vector<double> > A, int size) {
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

vector<float> thomas_f(vector<float> a, vector<float> b,
                     vector<float> c, vector<float> d){
    int n = a.size();
    vector<float> p;
    vector<float> q;
    p.push_back(c[0]/b[0]);
    q.push_back(d[0]/b[0]);
    for(int j=1; j<n; j++){
        if(j<n-1){
            float pj = c[j] / (b[j] - a[j]*p[j-1]);
            p.push_back(pj);
        }
        float qj = (d[j] - a[j]*q[j-1]) / (b[j] - a[j] * p[j-1]);
        q.push_back(qj);
    }
    vector<float> X(n,0.0);
    X[n-1] = q[n-1];   
    for(int j=n-2; j>=0; j--){
        float xj = q[j] - p[j]*X[j+1];
        X[j] = xj;
    }
    return X;
}

vector<double> thomas(vector<double> a, vector<double> b,
                     vector<double> c, vector<double> d){
    int n = a.size();
    vector<double> p;
    vector<double> q;
    p.push_back(c[0]/b[0]);
    q.push_back(d[0]/b[0]);
    for(int j=1; j<n; j++){
        if(j<n-1){
            double pj = c[j] / (b[j] - a[j]*p[j-1]);
            p.push_back(pj);
        }
        double qj = (d[j] - a[j]*q[j-1]) / (b[j] - a[j] * p[j-1]);
        q.push_back(qj);
    }
    vector<double> X(n,0.0);
    X[n-1] = q[n-1];   
    for(int j=n-2; j>=0; j--){
        double xj = q[j] - p[j]*X[j+1];
        X[j] = xj;
    }
    return X;
}

vector<double> gauss(vector< vector<double> > A) {
    int n = A.size();
    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -A[k][i]/A[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    vector<double> x(n);
    for (int i=n-1; i>=0; i--) {
        x[i] = A[i][n]/A[i][i];
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
}

vector<double> matrix_X(int size){
    vector<double> X(size);
    for(int i=0; i<size; i++){
        if(i%2==0) X[i] = -1.0;
        else X[i] = 1.0;
    }

    return X;
}

double euclid_norm(vector<double> X, vector<double> X_calc){
    double sum = 0;
    vector<double> matrix(X.size());
    for(int i=0; i<X.size(); i++){
        matrix[i] = X[i] - X_calc[i];
    }
    for(int i=0; i<matrix.size(); i++){
        sum += matrix[i]*matrix[i];
    }
    return sqrt(sum);
}

double max_norm(vector<double> X, vector<double> X_calc){
    vector<double> matrix(X.size());
    for(int i=0; i<X.size(); i++){
        matrix[i] = X[i] - X_calc[i];
    }
    double max = abs(matrix[0]);
    for(int i=1; i< matrix.size(); i++){
        if(abs(matrix[i])>max) max = abs(matrix[i]);
    }
    return max;
}

vector<float> gauss_f(vector< vector<float> > A) {
    int n = A.size();
    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        float maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            float tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            float c = -A[k][i]/A[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    vector<float> x(n);
    for (int i=n-1; i>=0; i--) {
        x[i] = A[i][n]/A[i][i];
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
}

vector<vector<double>> matrix_three(int size) {
    vector<double> row(size+1);
    vector<vector <double>> A(size,row); 
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            if(i==j) A[i][j] = 3.0;
            else if(i-1==j) A[i][j] = 3.0/(i+1+5+1);
                else if(i+1==j) A[i][j] = 1.0/(i+1+5);
                    else A[i][j] = 0.0;
            
        }
        A[i][size] = 0.0;
    }
    return A;

}



vector<vector<float>> matrix_three_f(int size) {
    vector<float> row(size+1);
    vector<vector <float>> A(size,row); 
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            if(i==j) A[i][j] = 3.0;
            else if(i-1==j) A[i][j] = 3.0/(i+1+5+1);
                else if(i+1==j) A[i][j] = 1.0/(i+1+5);
                    else A[i][j] = 0.0;
            
        }
        A[i][size] = 0.0;
    }
    return A;

}

vector<float> matrix_X_f(int size){
    vector<float> X(size);
    for(int i=0; i<size; i++){
        if(i%2==0) X[i] = -1.0;
        else X[i] = 1.0;
    }

    return X;
}

float euclid_norm_f(vector<float> X, vector<float> X_calc){
    float sum = 0;
    vector<float> matrix(X.size());
    for(int i=0; i<X.size(); i++){
        matrix[i] = X[i] - X_calc[i];
    }
    for(int i=0; i<matrix.size(); i++){
        sum += matrix[i]*matrix[i];
    }
    return sqrt(sum);
}

float max_norm_f(vector<float> X, vector<float> X_calc){
    vector<float> matrix(X.size());
    for(int i=0; i<X.size(); i++){
        matrix[i] = X[i] - X_calc[i];
    }
    float max = abs(matrix[0]);
    for(int i=1; i< matrix.size(); i++){
        if(abs(matrix[i])>max) max = abs(matrix[i]);
    }
    return max;
}


int main(int argc,  char** argv) {
    int size = atoi(argv[1]);
    cout.precision(5);  
    vector<vector<double>> A = matrix_three(size);
    vector<double> X = matrix_X(size);
    vector<double> B_T(size,0.0);
    cout << "Calculating using Gauss with size: " << size << "\n";
    //Calculate vector B
    for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
            B_T[i] += A[i][j] *X[j];
            A[i][size] += A[i][j] * X[j];
        }
    }   

    vector<double> X_calc(size);

    // Calculate solution
    clock_t start = clock();
    X_calc = gauss(A);
    clock_t time = clock() - start;
    cout << endl;
    //Calclating the norm
    double euclid = euclid_norm(X,X_calc);
    double max = max_norm(X,X_calc);
    cout << "Error using euclid norm: " << euclid << "\n";
    cout << "Error using max norm: " << max << "\n";
    cout.precision(15);
    cout << "Time to execute: " << time/((double)CLOCKS_PER_SEC/1000)<<endl<<endl;

    cout << "Calculating using Thomas with size: " << size << "\n";
    //Calculate vector B
    vector<double> subdiagonal(size);
    vector<double> diagonal(size);
    vector<double> superdiagonal(size);
    subdiagonal[0] = 0.0;
    superdiagonal[size-1] = 0.0;
    for(int i=0; i< size; i++){
        diagonal[i] = 3.0;
        if(i>0) subdiagonal[i] = 3.0/(i+1+5+1);
        if(i<size-1) superdiagonal[i] = 1.0/(i+1+5);
        //B_T[i] = A[i][size];
    }
    

    vector<double> X_calc_T(size);

    // Calculate solution
    start = clock();
    X_calc_T = thomas(subdiagonal,diagonal,superdiagonal,B_T);
    time = clock() - start;
    cout << endl;
    
    //Calclating the norm
    double euclid_T = euclid_norm(X,X_calc_T);
    double max_T = max_norm(X,X_calc_T);
    cout.precision(5);
    cout << "Error using euclid norm: " << euclid_T << "\n";
    cout << "Error using max norm: " << max_T << "\n";
    cout.precision(15);
    cout << "Time to execute: " << time/((double)CLOCKS_PER_SEC/1000) << endl;
    cout << endl << "--------------------------------------------------------------" << endl;
    return 0;
}