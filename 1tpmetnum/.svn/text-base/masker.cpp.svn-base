#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <ieee754.h>

#define TIPO long double
#define BUF_SIZE 65
#define MAX_MANTISA 52


using namespace std;



char *int2bin(int lw, int up, char *buffer, int buf_size) {
    buffer += (buf_size - 1);

    for (int i = 63; i >= 32; i--) {
        *buffer-- = (up & 1) + '0';

        up >>= 1;
    }
    
    for (int i = 31; i >= 0; i--) {
        *buffer-- = (lw & 1) + '0';

        lw >>= 1;
    }

    return buffer;
}


void mostrar(float *x)
{
    unsigned char c0 = *((char*)x + 3);
    unsigned char c1 = *((char*)x + 2);
    unsigned char c2 = *((char*)x + 1);
    unsigned char c3 = *((char*)x);

    char todo[33];

    unsigned char pot = 128;
    for(int i=0; i<8; ++i)
    {
        todo[i] = '0' + (c0 / pot) % 2;
        pot /= 2;
    }

    pot = 128;
    for(int i=8; i<16; ++i)
    {
        todo[i] = '0' + (c1 / pot) % 2;
        pot /= 2;
    }

    pot = 128;
    for(int i=16; i<24; ++i)
    {
        todo[i] = '0' + (c2 / pot) % 2;
        pot /= 2;
    }

    pot = 128;
    for(int i=24; i<32; ++i)
    {
        todo[i] = '0' + (c3 / pot) % 2;
        pot /= 2;
    }

    todo[32] = 0;

    char signo = todo[0];
    char exponente[9];
    char mantisa[24];

    int _exponente = 0;
    float _mantisa = 1.0;

    for(int i=0, pot=128; i<8; ++i)
    {
        exponente[i] = todo[i+1];
        _exponente += (exponente[i] - '0') * pot;
        pot /= 2;
    }

    float _pot = 0.5;
    for(int i=0; i<23; ++i)
    {
        mantisa[i] = todo[i+9];
        _mantisa += (mantisa[i] - '0') * _pot;
        _pot /= 2;
    }

    exponente[8] = 0;
    mantisa[23] = 0;

    printf( "signo = %c, ", signo );
    printf( "exp = %s (%d), ", exponente, _exponente );
    printf( "mant = 1.%s \n", mantisa );
//  printf( "%3.5f = %s \n", *x, todo );
}



//ieee754_double generate_double(int mantissa_size){
	
	////mantissa_size;
	////int upper = 0xFFFFFFFF;
	////int lower = 0xFFFFFFFF;
	
	
	
	////unsigned int cantshifts = MAX_MANTISA - mantisa_size;
	
	////if (cantshifts > 32){
		////for (unsigned int i = 0; i < max (cantshifts - 32, (unsigned int)0); i++){
			////upper = upper << 1;
		////}
	////}
	
	////for(unsigned int i = 0; i < min ((unsigned int)32, cantshifts); i++){
		////lower = lower << 1;
	////}


	
	////char buffer[BUF_SIZE];
    ////buffer[BUF_SIZE - 1] = '\0';

    ////int2bin(upper, lower, buffer, BUF_SIZE - 1);
    //ieee754_double encoded;
    //encoded.ieee.mantissa0 = 0xFFFFF;
    //encoded.ieee.mantissa1 = 0xFFFFFFFF;


    //int cantshifts = MAX_MANTISA - mantissa_size;
    
    //printf("Son %u shifts\n", cantshifts);
        
    //int upper_mascara = 0xFFFFFFFF;
    //int lower_mascara = 0xFFFFFFFF;
    
    //printf("La primer mantisa es 0x%x, la segunda es 0x%x\n", encoded.ieee.mantissa0, encoded.ieee.mantissa1);
    
    //printf("La primer mascara es 0x%x, la segunda es 0x%x\n", upper_mascara, lower_mascara);
    
    
    
    //if (cantshifts > 32){
	//upper_mascara = upper_mascara << min(20, cantshifts-32);
    //}
    //if (cantshifts >= 32){
	//lower_mascara = 0;
    //}else{
	//lower_mascara = lower_mascara << cantshifts;
   //}
   
   
   //printf("La primer mascara es 0x%x, la segunda es 0x%x\n", upper_mascara, lower_mascara);
     
    ////encoded.d = a;
    //encoded.ieee.mantissa0 = encoded.ieee.mantissa0 & upper_mascara;
    //encoded.ieee.mantissa1 = encoded.ieee.mantissa1 & lower_mascara; 
       
       
    //// mantissa0: 20bits
    //// mantissa1: 32bits
     
     
    
     
     
     
    ////bla = bla < (32- (12+mantisa_size));
 
    ////printf("%d = %s\n", mantisa_size, buffer);
    //printf("La primer mantisa es 0x%x, la segunda es 0x%x\n", encoded.ieee.mantissa0, encoded.ieee.mantissa1);

////	mostrar(&bla);
	
	
	//return encoded;
//}

int main (int argv, char* argc[]){

	int mantisa_size = atoi(argc[1]);

	long long int mask = -1;
	
	unsigned int cantshifts = MAX_MANTISA - mantisa_size;
	
	for(unsigned int i = 0; i < cantshifts; i++){
		mask = mask << 1;
	}
	
	printf("0x%.llX\n",mask);
	return 0;	
	

	
}
