#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f1(double x){
	return pow(x,8)-8*pow(x,7)+28*pow(x,6)-56*pow(x,5)+70*pow(x,4)-56*pow(x,3)+28*pow(x,2)-8*x+1;
}
double f2(double x){
	return (((((((x-8)*x+28)*x-56)*x+70)*x-56)*x+28)*x-8)*x+1;
}
double f3(double x){
	return pow((x-1),8);
}
double f4(double x){
	if(x!=1)
		return exp(8*log(abs(x-1)));
	return 0;
}

int main(){
    long long int x=9900;
    cout << setprecision(16);
	for(double i=0.99;i<=1.01;i+=0.0002){

		cout<<i<<" | ";
		cout<<f1(i)<<" | ";
		cout<<f2(i)<<" | ";
		cout<<f3(i)<<" | ";
//		cout<<f3t(x)<<" | ";
		cout<<f4(i)<<endl;
        x+=2;
	}
	return 0;
}
