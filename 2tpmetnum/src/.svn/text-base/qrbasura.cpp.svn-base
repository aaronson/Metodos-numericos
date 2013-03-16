#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <iomanip>

#define precision 2

using namespace std;

typedef vector< vector<double> > vvd;
typedef vector<double> vd;



int cant_radios, cant_angulos, radio_int, radio_ext;
int cant_ecus;

double radio_dist, angulo_dist;


double obtener_valor_radio (int radio);
vvd matriz_identidad(int n);
vvd submatriz(int n, vvd mat);
vvd completar_identidad(int n, vvd mat);
void descomposicion_QR (vvd matriz);
vvd producto_mat (vvd m1, vvd m2);
void trasponer(vvd &mat);
vd producto_matvec(vvd m, vd v);
vd resolver_triangulado(vvd r, vd v);



vd temps_exterior;
vvd matrix;
vvd q;
vvd r;
vd vecsol;

void inicializar_matrix();
void imprimir_matrix(vvd mat);

int main(){

	cout << setprecision(precision);

	cin >> cant_radios >> cant_angulos >> radio_int >> radio_ext;

	//cant_radios = n, cant_angulos = m
	//double temps_exterior[cant_angulos];


	for(int i = 0; i<cant_angulos; i++){
		double temp;
		cin >> temp;
		temps_exterior.push_back(temp);
	}

	radio_dist = (radio_ext-radio_int) / (cant_radios-1);
	angulo_dist = 2*M_PI / cant_angulos;

	cant_ecus = cant_angulos * cant_radios;

	matrix.resize(cant_ecus);
	vecsol.resize(cant_ecus);
	
	for (int i = 0; i < cant_ecus; i++){
		vecsol[i] = 1;
	}

	inicializar_matrix();
	r = matrix;
	q = matriz_identidad(cant_ecus);


	cout << "ESTA es Q" << endl;
	imprimir_matrix(q);

	cout << "ESTA es R" << endl;
	imprimir_matrix(r);


	descomposicion_QR(matrix);

	//imprimir_matrix(matrix);
	
	//vvd res = producto_mat(q,r);

	trasponer(q);
	
	vd res = producto_matvec(q,vecsol);
	res[5] = 0.1f;
	for (unsigned int i = 0; i < res.size(); i++){
		cout << "Esta es res" << i << ": " << res[i] << endl;
	}
	
	vd soluciones = resolver_triangulado(r,res);
	
	cout << "Estas son las soluciones" << endl;

	for (unsigned int i = 0; i < soluciones.size(); i++){
		cout << "Esta es x" << i << ": " << soluciones[i] << endl;
	}
	
	
	
	return 0;
};

int convertir_coord (int rad, int ang){
	if (ang == -1){
		return rad * cant_angulos + (cant_angulos-1);
	}else if (ang == cant_angulos){
		return rad * cant_angulos;
	} else {
		return rad * cant_angulos + ang;
	}
}

void calcular_coefs(int neq){
		int nradio = neq / cant_radios;
		int nangulo = neq % cant_radios;

		int vradio = obtener_valor_radio(nradio);

		double adentro = (vradio - radio_dist)/(radio_dist*radio_dist*vradio);

		matrix[neq][convertir_coord(nradio - 1, nangulo)] =  adentro;

		double afuera = (1/(radio_dist*radio_dist));
		matrix[neq][convertir_coord(nradio + 1, nangulo)] =  afuera;

		double costado = (1/(vradio*vradio*angulo_dist*angulo_dist));
		matrix[neq][convertir_coord(nradio, nangulo+1)] =  costado;
		matrix[neq][convertir_coord(nradio, nangulo-1)] =  costado;

		double centro = ((-2*vradio*vradio*angulo_dist*angulo_dist) +
		(vradio * radio_dist * angulo_dist*angulo_dist) - (2*radio_dist*radio_dist)) / (vradio*vradio*angulo_dist*angulo_dist*radio_dist*radio_dist);
		matrix[neq][convertir_coord(nradio, nangulo)] =  centro;

}

double obtener_valor_radio (int radio){
	return radio*radio_dist+radio_int;
}

void inicializar_matrix(){

	for(int i = 0; i < cant_ecus; i++){
		matrix[i].resize(cant_ecus);
	}

	for (int i = 0; i < cant_angulos; i++){
		matrix[i][i] = 1;
		vecsol[i] = temps_exterior[i];
	}

	for (int i = 0; i < cant_angulos; i++){
		int indice = cant_angulos*(cant_radios-1) + i;
		matrix[indice][indice] = 1;
		vecsol[indice] = 1500;
	}


	for (int j = cant_angulos; j < cant_ecus - cant_angulos; j++ ){
		//~ for (int k = 0; k < cant_angulos; k++){
			//~
		//~ }
		calcular_coefs(j);
	}
}

vd extraer_col(vvd mat, int ncol){
	vd res;
	int tam = mat.size();
	for(int i = 0; i < tam; i++){
			res.push_back(mat[i][ncol]);
	}
	return res;
}


double calc_norma (vd vec){
	double suma_cuadrado = 0;
	for (unsigned int i = 0; i<vec.size(); i++){
        suma_cuadrado+= vec[i]*vec[i];
	}
	suma_cuadrado = sqrt(suma_cuadrado);

return suma_cuadrado;}

void trasponer(vvd &mat){
	for (unsigned int i = 0; i < mat.size(); i++){
		for (unsigned int j = 0; j < i; j++){
			double temp;
			temp = mat[i][j];
			mat[i][j] = mat[j][i];
			mat[j][i] = temp;
		}
	}
}

void imprimir_matrix(vvd mat){
    int tam = mat.size();
	for (int i = 0; i < tam; i++){
		for(int j = 0; j < tam; j++){
			if (mat[i][j] > 99){
				cout << setprecision(4);
				cout << mat[i][j] << "\t";
				cout << setprecision(precision);
				}else if ((mat[i][j] < 0.0001) && (mat[i][j] > -0.0001)){
					cout << 0 << "\t";
			} else
			cout << mat[i][j] << "\t";
	}
	if (vecsol[i] > 99){
		cout << setprecision(4);
		cout << "= " <<  vecsol[i]<< endl;
		cout << setprecision(precision);
	}else
		cout << "= " <<  vecsol[i]<< endl;
	}
}

vvd matriz_identidad(int n){
	vvd matriz;
	matriz.resize(n);

	for(int i = 0; i < n; i++){
		matriz[i].resize(n);
	}

	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			if (j==i){
          		matriz[i][j] = 1;
 			}else{
	    		matriz[i][j] = 0;
			}
		}
	}
	return matriz;
}

vd resolver_triangulado(vvd r, vd v){
	int tam = v.size();
	vd res(tam,0);
	res[tam-1] = v[tam-1];
	for (int i = tam-2; i >= 0; i--){
		if ((i>=cant_angulos) && (i<cant_ecus-cant_angulos) && (v[i] != r[i][i])){
			double acum = 0;
			for (int j = tam-1; j > i; j--){
				acum += r[i][j]*res[j];
			}
			res[i] = acum / (v[i]-r[i][i]);
		}else{
			res[i] = v[i];
		}
	}
	return res;
}

vd producto_matvec(vvd m, vd v){
	int tam = m.size();
	vd res(tam,0);
	for(int i = 0; i<tam; i++){
		double acum=0;
		for (int j=0; j<tam; j++){
			acum += m[i][j]*v[j];
		
		}
		res[i]=acum;
	}
	return res;
}

vvd producto_mat (vvd m1, vvd m2){
	vvd res;

	int tam = m1.size();

	res.resize(tam);

	for(int i = 0; i < tam; i++){
		res[i].resize(m2.size());
	}

	for(int i = 0; i < tam; i++){
		for(int j = 0; j < tam; j++){
			double acum = 0;
			for (int k = 0; k < tam; k++){
				acum += m1[i][k]*m2[k][j];
			}
			res[i][j] = acum;
		}
	}
  return res;
}


vvd resta (vvd m1, vvd m2){
	int tam = m1.size();
	vvd res;
	res.resize(tam);
	for(int i = 0; i < tam; i++){
		res[i].resize(tam);
	}
	for(int i = 0; i < tam; i++){
		for(int j = 0; j < tam; j++){
			res[i][j] = m1[i][j] - m2[i][j];
		}
	}

	return res;
}

vd resta_vec(vd v1, vd v2){
	int tam = v1.size();
	vd resta(tam,0);
	for(int i = 0; i < tam; i++){
			resta [i]=v1[i]-v2[i];
	}
	return resta;
}

vvd matricear(vd v){
	vvd res;
	int tam = v.size();
    res.resize(tam);
	for(int i = 0; i < tam; i++){
		res[i].resize(tam);
	}
	for(int i = 0; i < tam; i++){
		for(int j = 0; j < tam; j++){
			res[i][j] = v[i]*v[j];
		}
	}
	return res;

}

vvd mult_escalar (vvd m1, double n){
	vvd mult;
	int tam = m1.size();
	mult.resize(tam);
	for(int i = 0; i < tam; i++){
		mult[i].resize(tam);
	}
	for(int i = 0; i < tam; i++){
		for(int j = 0; j < tam; j++){
			mult[i][j]= m1[i][j] * n;
		}
	}
	return mult;
}

bool ya_triangulado(vd x){
    for (unsigned int i = 1; i<x.size(); i++){
    	if (x[i] != 0) return false;
    }
    return true;
}

vvd householdear (vvd aux){
	vd x = extraer_col(aux,0);
	if (ya_triangulado(x)) return aux;
	double norma = calc_norma(x);
	vd y (aux.size(),0);
	y[0]=norma;
	vd u(aux.size(),0);
	u = resta_vec(x,y); // u = x-y
	norma = calc_norma(u);
	vvd mult = matricear(u); // mult = u * ut quedando una matriz de nxn
	mult = mult_escalar(mult,2); // mult = 2 x u * ut           ut= u traspuesto
	double temp = norma*norma;
	mult = mult_escalar(mult,(1/temp)); // mult = 2 x u*ut / norma^2
	vvd h; // householder
	h = resta(matriz_identidad(aux.size()),mult); // householder de la primer columa
    cout << "ESTA ES H"<<endl;
    imprimir_matrix(h);
	return h;
}


void descomposicion_QR (vvd matriz){
	//vvd aux = matriz;
	vvd aux;
	int tam = matriz.size();
	//vvd aux3 = matriz_identidad(tam);
	for (int i=0; i<tam-1; i++){
		aux = completar_identidad(i,(householdear(submatriz(i,r))));
		r = producto_mat(aux, r);
		for (int j = i+1; j<tam; j++){
			r[j][i] = 0;
		}
		trasponer(aux);
		q = producto_mat(q, aux);


		cout << "ESTA es Q" << endl;
		imprimir_matrix(q);

		cout << "ESTA es R" << endl;
		imprimir_matrix(r);
		/*
		 Una veez que se tiene el householder lo que se hace es multiplicarlo por la matriz original. Aca obtenemos una matriz
		 * nueva con la columna 1 en 0 salvo el A [1][1], entonces se toma una submatriz de n-1 x n-1 y se aplica lo mismo que ahora.
		 * Vamos a obtener h2' q es de tamanio mas chico al que necesitamos, entonces se completa h2 quedando la matriz n x n con la forma
		 * de la identidad en las filas sacadas y h2' en el medio
		 *
		 *                         1 0 0
								   0 H2'
		                           0
		 ahi tenemos lo que seria H2 y asi hasta que queda la matriz triangulada..
		 Eso es lo que me esta costando hacer.
		 Es importante acumular las H que vamos obteniendo xq
		 * R = prod Hi * A
		 * Q = prod Hi trans
		 */



		/*
		r = producto_mat(r, resta(matriz_identidad(cant_ecus), mult));
		q = producto_mat(q, trasponer(resta(matriz_identidad(cant_ecus), mult)));
		*/
     }

}

vvd submatriz(int n, vvd mat){
	int tam = mat.size() - n;
	vvd res;
	res.resize(tam);

	for (int i = 0; i < tam; i++){
		res[i].resize(tam);
	}

	for (int i = 0; i < tam; i++){
		for (int j = 0; j < tam; j++){
			res[i][j] = mat[n+i][n+j];
		}
	}
	return res;
}

vvd completar_identidad(int n, vvd mat){
    if (n == 0) return mat;

    vvd res;
	int tam_orig = mat.size();

	res.resize(n+tam_orig);


	for (int i = 0; i < n+tam_orig; i++){
		res[i].resize(n+tam_orig);
	}


	for (int i = 0; i < n; i++){
			res[i][i] = 1;
	}

	for (int i = 0; i < tam_orig; i++){
		for (int j = 0; j < tam_orig; j++){
			res[i+n][j+n] = mat[i][j];
		}
	}
	return res;

}
