// my first program in C++
#include <iostream>
#include <math.h>
using namespace std;

double u(double x, double y){
	return sin(x*M_PI)*cos(y*M_PI);
}

double f(double x, double y){
	return -2*M_PI*M_PI*sin(M_PI*x)*cos(M_PI*y);
}

int main()
{

	int N=320;
	double h=M_PI/N;
	double w=150./100;
	double x[N+1];
	for (int i = 0; i < N+1; ++i){
		x[i]=M_PI/2+h*i;
	}
	double num[N+1][N+1];
	double num1[N+1][N+1];
	for (int i = 0; i < N+1; ++i){
		for (int j = 0; j < N+1; ++j){
			num[i][j]=0;
			num1[i][j]=0;
		}
	}

	double eps=1;

	while(eps>1e-14){
		for (int i = 0; i < N+1; ++i){
			num1[0][i]=u(M_PI/2,x[i]);
			num1[N][i]=u(3*M_PI/2,x[i]);
			num1[i][0]=u(x[i],M_PI/2);
			num1[i][N]=u(x[i],3*M_PI/2);
		};
		for (int i = 1; i < N; ++i){
			for (int j = 1; j < N; ++j){
				num1[i][j]=-f(x[i],x[j])*h*h/4+(num1[i][j-1]+num1[i-1][j]+num[i][j+1]+num[i+1][j])/4;
			}
		};
		eps=0;
		for (int i = 1; i < N; ++i){
			for (int j = 1; j < N; ++j){
				num1[i][j]=w*num1[i][j]+(1-w)*num[i][j];
				eps=max(eps,fabs(num1[i][j]-num[i][j]));
			}
		}
		for (int i = 0; i < N+1; ++i){
			for (int j = 0; j < N+1; ++j){
				num[i][j]=num1[i][j];
			}
		}
		cout << eps << endl;
	}

	double eps1=0;

	for (int i = 1; i < N+1; ++i){
		for (int j = 1; j < N+1; ++j){
			eps1=max(eps1,fabs(num1[i][j]-u(x[i],x[j]))); 
		}
	}
	cout << eps1;

}




