#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

void inverse(vector<vector<float>> A, vector<vector<float>> inverse);
int determinant(vector<vector<float>> A, int n);
void adjoint(vector<vector<float>> A,vector<vector<float>> adj);
void getCofactor(vector<vector<float>> A, vector<vector<float>> temp, int p, int q, int n);


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



vector<vector<double>> matrix_one(int size) {
    vector<double> row(size+1);
    vector<vector <double>> A(size,row); 
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            if(i==0) A[i][j] = 1.0;
            else A[i][j] = 1.0/(i+1+j+1-1);
            
        } 
        A[i][size] = 0.0;
    }
    return A;

}


vector<vector<double>> matrix_two(int size){
    vector<double> row(size+1);
    vector<vector <double>> A(size,row);
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            if(i<=j) A[i][j] = 2.0*(i+1)/(j+1);
            else A[i][j] = A[j][i];
        }
        A[i][size] = 0.0;
    }

    return A;
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



vector<vector<float>> matrix_one_f(int size) {
    vector<float> row(size+1);
    vector<vector <float>> A(size,row); 
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            if(i==0) A[i][j] = 1.0;
            else A[i][j] = 1.0/(i+1+j+1-1);
            
        } 
        A[i][size] = 0.0;
    }
    return A;

}


vector<vector<float>> matrix_two_f(int size){
    vector<float> row(size+1);
    vector<vector <float>> A(size,row);
    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            if(i<=j) A[i][j] = 2.0*(i+1)/(j+1);
            else A[i][j] = A[j][i];
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
    vector<vector<double>> A = matrix_one(size);
    vector<double> X = matrix_X(size);

    cout << "Calculating on doubles with size: " << size << "\n";
    //Calculate vector B
    for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
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

    cout << "Calculating on floats with size: " << size << "\n";
    //Calculate vector B
    vector<vector<float>> A_f = matrix_one_f(size);

    vector<float> X_f = matrix_X_f(size);
    for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
            A_f[i][size] += A_f[i][j] * X_f[j];
        }
    }   

    vector<float> X_calc_f(size);

    // Calculate solution
    start = clock();
    X_calc_f = gauss_f(A_f);
    time = clock() - start;
    cout << endl;
    //Calclating the norm
    float euclid_f = euclid_norm_f(X_f,X_calc_f);
    float max_f = max_norm_f(X_f,X_calc_f);
    cout.precision(5);
    cout << "Error using euclid norm: " << euclid_f << "\n";
    cout << "Error using max norm: " << max_f << "\n";
    cout.precision(15);
    cout << "Time to execute: " << time/((double)CLOCKS_PER_SEC/1000) << endl;

    return 0;
}



// Function to get cofactor of A[p][q] in temp[][]. n is current 
// dimension of A[][] 
void getCofactor(vector<vector<float>> A, vector<vector<float>> temp, int p, int q, int n) 
{ 
    int i = 0, j = 0; 
  
    // Looping for each element of the matrix 
    for (int row = 0; row < n; row++) 
    { 
        for (int col = 0; col < n; col++) 
        { 
            //  Copying into temporary matrix only those element 
            //  which are not in given row and column 
            if (row != p && col != q) 
            { 
                temp[i][j++] = A[row][col]; 
  
                // Row is filled, so increase row index and 
                // reset col index 
                if (j == n - 1) 
                { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    } 
} 
  
/* Recursive function for finding determinant of matrix. 
   n is current dimension of A[][]. */
int determinant(vector<vector<float>> A, int n) 
{   
    int N = A.size();
    int D = 0; // Initialize result 
  
    //  Base case : if matrix contains single element 
    if (n == 1) 
        return A[0][0]; 
    
    vector<float> row(0);
    vector<vector<float>> temp(N,row); // To store cofactors 
  
    int sign = 1;  // To store sign multiplier 
  
     // Iterate for each element of first row 
    for (int f = 0; f < n; f++) 
    { 
        // Getting Cofactor of A[0][f] 
        getCofactor(A, temp, 0, f, n); 
        D += sign * A[0][f] * determinant(temp, n - 1); 
  
        // terms are to be added with alternate sign 
        sign = -sign; 
    } 
  
    return D; 
} 
  
// Function to get adjoint of A[N][N] in adj[N][N]. 
void adjoint(vector<vector<float>> A,vector<vector<float>> adj) 
{ 
    int N = A.size();
    if (N == 1) 
    { 
        adj[0][0] = 1; 
        return; 
    } 
  
    // temp is used to store cofactors of A[][] 
    int sign = 1; 
    vector<float> row(0);
    vector<vector<float>> temp(N,row);
    for (int i=0; i<N; i++) 
    { 
        for (int j=0; j<N; j++) 
        { 
            // Get cofactor of A[i][j] 
            getCofactor(A, temp, i, j, N); 
  
            // sign of adj[j][i] positive if sum of row 
            // and column indexes is even. 
            sign = ((i+j)%2==0)? 1: -1; 
  
            // Interchanging rows and columns to get the 
            // transpose of the cofactor matrix 
            adj[j][i] = (sign)*(determinant(temp, N-1)); 
        } 
    } 
} 
  
// Function to calculate and store inverse, returns false if 
// matrix is singular 
void inverse(vector<vector<float>> A, vector<vector<float>> inverse) 
{ 
    int N = A.size();
    // Find determinant of A[][] 
    int det = determinant(A, N); 
    
  
    // Find adjoint 
    vector<float> row(0);
    vector<vector<float>> adj(N,row); 
    adjoint(A, adj); 
  
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
    for (int i=0; i<N; i++) 
        for (int j=0; j<N; j++) 
            inverse[i][j] = adj[i][j]/float(det); 
  
    
} 