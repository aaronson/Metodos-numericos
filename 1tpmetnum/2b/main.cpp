#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#define precision 50
#define MAX_MANTISA 52


using namespace std;

int upper;
int lower;

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

double domask (double x){
	
	int* b = (int*) &x;
		
	*b = *b & lower;
	*(b+1) = *(b+1) & upper;
	
	return *((double *) b);
}

int main(int argc, char* argv[]){
	
	cout << setprecision(precision);

	int terms = atoi(argv[2]);
	double x =  atof (argv[1]);
	
	int mantisa_size = atoi(argv[3]);
	unsigned int cantshifts = MAX_MANTISA - mantisa_size;
	
	upper = 0xFFFFFFFF;
	lower = 0xFFFFFFFF;
		
	if (cantshifts > 32){
		for (unsigned int i = 0; i < max (cantshifts - 32, (unsigned int)0); i++){
			upper = upper << 1;
		}
	}
	
	for(unsigned int i = 0; i < min ((unsigned int)32, cantshifts); i++){
		lower = lower << 1;
	}
	
	
	bool signo = true;
	
	double acum = domask(0);
	double temppot, tempfact;
				
	for (int i = 0; i < terms; i++){
		if (signo){
			temppot = domask (potencia (x, (2*i)));
			tempfact = domask ((double) (factorial (2*i)));
			
			temppot = domask (temppot/tempfact);
			
			acum += temppot;
			
			acum = domask (acum);
			
		}else{
			temppot = domask (potencia (x, (2*i)));
			tempfact = domask ((double) (factorial (2*i)));
			
			temppot = domask (temppot/tempfact);
			
			acum -= temppot;
			
			acum = domask (acum);
		}
		signo = !signo;
	}
	

	cout << x << "," << terms << "," << mantisa_size  << "," << acum << "," << cos(x) << endl;
	
	return 0;
};




