#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#define precision 20

using namespace std;

double factorial (int x){
	
	double res = 1;
	
	for (int i = 1; i<=x; i++){
		res *= i;
	}	
	
	return res;
}

double potencia(double num, int exp){
	
	double res = 1;

	for (int i = 0; i < exp; i++){
		res *= num;
	}
	
	return res;
}

int main(int argc, char* argv[]){
	
	cout << setprecision(precision);

	int terms = atoi(argv[2]);
	double x =  atof (argv[1]);
	
	bool signo = true;
	
	double acum = 0;
				
	for (int i = 0; i < terms; i++){
		if (signo){
			acum +=  (potencia (x, (2*i)) / (double) (factorial (2*i)));
		}else{
			acum -=  (potencia (x, (2*i)) / (double) (factorial (2*i)));
		}
		signo = !signo;
	}
	cout << x << "," << terms << "," << acum << "," << cos(x) << endl;
	
	
};



