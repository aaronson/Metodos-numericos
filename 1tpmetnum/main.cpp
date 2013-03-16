#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>


#define TIPO long double
#define precision 20

using namespace std;

double factorial (int x){
	
	double res = 1;
	
	for (int i = 1; i<=x; i++){
		res *= i;
	}	
	
	return res;
}

TIPO potencia(TIPO num, int exp){
	
	TIPO res = 1;

	for (int i = 0; i < exp; i++){
		res *= num;
	}
	
	return res;
}

int main(int argc, char* argv[]){
	
	cout << setprecision(precision);
	
	int terms;
	TIPO x;
	
	bool signo = true;
	
	TIPO acum = 0;
	
	cin >> terms >> x;
			
	for (int i = 0; i < terms; i++){
		if (signo){
			acum +=  (potencia (x, (2*i)) / (TIPO) (factorial (2*i)));
		}else{
			acum -=  (potencia (x, (2*i)) / (TIPO) (factorial (2*i)));
		}
		signo = !signo;
		cout << acum << endl;
	}
	
	
};


