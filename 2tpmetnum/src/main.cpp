#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <iomanip>

#define precision 6

using namespace std;

typedef vector< vector<double> > vvd;
typedef vector<double> vd;



int cant_radios, cant_angulos;
double radio_int, radio_ext;
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
void gauss(vvd mat);
void calcular_iso(vd params);



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
		vecsol[i] = 0;
	}

	inicializar_matrix();

	//~ cout << "BEFORE GAUSS" << endl;
	//~ imprimir_matrix(matrix);
//~ 
	 gauss(matrix);
	//~ cout << "AFTER GAUSS" << endl;
//~ 
	//~ imprimir_matrix(matrix);
//~ 
	vd res= resolver_triangulado(matrix, vecsol);

	//~ for (unsigned int i =0; i< res.size(); i++ ){
        //~ cout << "El valor x" << i << " es: " << res[i] << endl;
    //~ }

    calcular_iso(res);

	return 0;
};

void pivotear(int i){
	double temp;
	int j = i;
	while (matrix[i][j] == 0) j++;

	for (int k = i; k < cant_ecus; k++){
		temp = matrix[i][k];
		matrix[i][k] = matrix[j][k];
		matrix[j][k] =  temp;
	}

	temp = vecsol[i];
	vecsol[j] = temp;
	vecsol[i] = vecsol[j];
}

void gauss(vvd mat){\
	double factor;
	int tam = matrix.size();
	for (int i = 0; i < tam-1; i++){
		if (matrix[i][i] == 0) pivotear(i);
		for (int k = i+1; k < tam; k++){
			if (matrix[k][i] != 0){
				factor = matrix[k][i]/matrix[i][i];
				for (int j= 0; j < tam; j++){
					matrix[k][j] -= matrix[i][j]*factor;
				}
				vecsol[k] -= vecsol[i]*factor;
			}
		}
	}
}

void calcular_iso(vd params){
	vd res;
	vvd mat;
	mat.resize(cant_radios);
	for(int i =0; i<cant_radios; i++){
		mat[i].resize(cant_angulos);
	}
	for(int i =0; i<cant_radios; i++){
		for(int j =0; j<cant_angulos; j++){
			mat[i][j] = params[i*cant_angulos+j];
		}
	}

	for(int j =0; j<cant_angulos; j++){
		int i=0;
		while (mat[i][j]<400) i++;
		if (mat[i][j]==400){
			res.push_back(obtener_valor_radio(i));
		}else{
			double vmenor = mat[i-1][j];
			double vmayor = mat[i][j];
			double rmenor = obtener_valor_radio(i-1);
			double rmayor = obtener_valor_radio(i);

			double radio400 = rmenor + (400-vmenor)*(rmayor-rmenor)/(vmayor-vmenor);

			res.push_back(radio400);

		}
	}

	for (unsigned int i =0; i< res.size(); i++ ){
        cout << res[i] << endl;
    }


}

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
		int nradio = neq / cant_angulos;
		int nangulo = neq % cant_angulos;

		double vradio = obtener_valor_radio(nradio);

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
	return (cant_radios-1-radio)*radio_dist+radio_int;
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
		calcular_coefs(j);
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
			acum-=vecsol[i];
			res[i] = acum / (-r[i][i]);
		}else{
			res[i] = v[i];
		}
	}
	return res;
}
