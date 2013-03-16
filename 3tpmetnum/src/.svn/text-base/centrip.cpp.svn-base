#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <iomanip>

#define precision 6
#define forn(i,n) for(int i = 0;i<(int)(n);i++)
#define CANT_SIMPSON 6
#define ERROR 0.1

using namespace std;

double delta;
int cant_puntos;

class Punto {
public:
	double x;
	double y;
	Punto(double _x, double _y) :
			x(_x), y(_y) {
	}
	Punto operator =(const Punto &p) {
		x = p.x;
		y = p.y;
		return *this;
	}
};

typedef vector<double> vd;
typedef vector<vd> vvd;
typedef vector<Punto> vp;
vd x, y, t;

void do_centripeta();
double calc_distancia(int i);
void gauss(vvd &mat, vd &v);
vd resolver_triangulado(vvd mat, vd v);
void imprimir_matrix(vvd mat);
void imprimir_vector(vd vec);

class Spline {
	vd a, b, c, d;
	int tam;
public:
	Spline(const Spline &s) :
			a(s.a), b(s.b), c(s.c), d(s.d), tam(s.tam) {
	}
	Spline(vd param);
	double evaluar(double x);
	double derivar(double x);
};

Spline::Spline(vd param) {
	tam = cant_puntos;

	vvd matrix(tam);
	for (int i = 0; i < tam; i++) {
		matrix[i].resize(tam);
	}

	vd h(tam - 1), alfa(tam), u(tam), l(tam), z(tam);
	alfa.resize(tam);
	a = param;
	b.resize(tam - 1);
	d.resize(tam - 1);

	for (int i = 0; i < tam - 1; i++) {
		h[i] = t[i + 1] - t[i];
	}

	for (int i = 1; i < tam - 1; i++) {
		alfa[i] = (3 * (param[i + 1] - param[i]) / h[i])
				- (3 * (param[i] - param[i - 1]) / h[i - 1]);
	}

	matrix[0][0] = 1;
	for (int i = 1; i < tam - 1; i++) {
		matrix[i][i - 1] = h[i - 1];
		matrix[i][i] = 2 * (h[i - 1] + h[i]);
		matrix[i][i + 1] = h[i];
	}
	matrix[tam - 1][tam - 1] = 1;

	gauss(matrix, alfa);

	c = resolver_triangulado(matrix, alfa);

	for (int i = tam - 2; i >= 0; i--) {
		b[i] = ((param[i + 1] - param[i]) / h[i])
				- (h[i] * (c[i + 1] + 2 * c[i]) / 3);
		d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
	}

	this->a = a;
	this->b = b;
	c.pop_back();
	this->c = c;
	this->d = d;
	this->tam = tam;
}

double Spline::evaluar(double x) {
	int i = 1;
	while (x > t[i]) {
		i++;
	}
	double x1 = x - t[i - 1];
	double acum = 0;
	acum += x1 * x1 * x1 * d[i - 1];
	acum += x1 * x1 * c[i - 1];
	acum += x1 * b[i - 1];
	acum += a[i - 1];
	return acum;

}

double Spline::derivar(double x) {
	int i = 1;
	while (x > t[i]) {
		i++;
	}
	double x1 = x - t[i - 1];

	double acum = 0;
	acum += x1 * x1 * 3 * d[i - 1];
	acum += x1 * 2 * c[i - 1];
	acum += b[i - 1];
	return acum;
}

class Curva {
public:

	Spline sx;
	Spline sy;

	Curva(Spline x, Spline y) :
			sx(x), sy(y) {
	}
	Punto evaluar(double x);
	double funcion_a_integrar(double x);
	double simpson_compuesto(double a, double b);
	vd resolver_spline();
	double aproximar_con_bb(double ant);
	vp resolver_param_original();
	void imprimir_vector(vd vec);
	void imprimir_puntos(vp ps);

}
;

vp Curva::resolver_param_original() {
	int i = 0;
	vp res;
	while (delta * i < t[cant_puntos - 1]) {
		res.push_back(evaluar(delta * i));
		i++;
	}
	return res;
}

Punto Curva::evaluar(double x) {
	return Punto(sx.evaluar(x), sy.evaluar(x));
}
;

double Curva::funcion_a_integrar(double x) {
	double dx = sx.derivar(x);
	double dy = sy.derivar(x);
	double bla = sqrt(dx * dx + dy * dy);
	return bla;
}

double Curva::simpson_compuesto(double a, double b) {
	int n = CANT_SIMPSON;
	double h = (b - a) / n; // toma el punto intermedio
	double x0 = funcion_a_integrar(a) + funcion_a_integrar(b);
	double x1 = 0;
	double x2 = 0;

	for (int i = 1; i < (int) (n); i++) {
		double x = a + i * h;
		if ((i % 2) == 0) {
			x2 += funcion_a_integrar(x);
		} else {
			x1 += funcion_a_integrar(x);
		}
	}

	x1 = h * (x0 + (2 * x2) + (4 * x1)) / 3;
	return x1;
}


double Curva::aproximar_con_bb(double ant) {
	double b;
	if (t[cant_puntos - 1] / 2 < ant) {
		b = t[cant_puntos - 1];
	} else {
		b = t[cant_puntos - 1] / 2 + ant;
	}
	double half = (b - ant);
	double res = simpson_compuesto(ant, b);
	if (b == t[cant_puntos - 1] && (res < delta))
		return 0;
	while (fabs(res - delta) > ERROR) {
		half /= 2;
		if (res > delta) {
			b -= half;
		} else {
			b += half;
		}
		res = simpson_compuesto(ant, b);
	}
	return b;
}

vd Curva::resolver_spline() {
	vd vsol;
	double sig = 0;
	vsol.push_back(sig);
	while (sig < t[cant_puntos - 1]) {
		sig = aproximar_con_bb(sig);
		if (sig == 0)
			break;
		vsol.push_back(sig);
	}
	return vsol;
}

vd resolver_triangulado(vvd mat, vd v) {
	int tam = v.size();
	vd res(tam, 0);
	res[tam - 1] = v[tam - 1];
	for (int i = tam - 2; i >= 0; i--) {
		double acum = 0;
		for (int j = tam - 1; j > i; j--) {
			acum += mat[i][j] * res[j];
		}
		acum -= v[i];
		res[i] = acum / (-mat[i][i]);
	}
	return res;
}

void gauss(vvd &mat, vd &v) {
	double factor;
	int tam = mat.size();
	for (int i = 0; i < tam - 1; i++) {
		for (int k = i + 1; k < tam; k++) {
			if (fabs(mat[k][i]) > 0.00001f) {
				factor = mat[k][i] / mat[i][i];
				for (int j = 0; j < tam; j++) {
					mat[k][j] -= mat[i][j] * factor;
				}
				v[k] -= v[i] * factor;
			}
		}
	}
}

//void do_hardcodeado() {
//	t.resize(21);
//	t[0] = 0.9;
//	t[1] = 1.3;
//	t[2] = 1.9;
//	t[3] = 2.1;
//	t[4] = 2.6;
//	t[5] = 3.0;
//	t[6] = 3.9;
//	t[7] = 4.4;
//	t[8] = 4.7;
//	t[9] = 5.0;
//	t[10] = 6.0;
//	t[11] = 7.0;
//	t[12] = 8.0;
//	t[13] = 9.2;
//	t[14] = 10.5;
//	t[15] = 11.3;
//	t[16] = 11.6;
//	t[17] = 12.0;
//	t[18] = 12.6;
//	t[19] = 13.0;
//	t[20] = 13.3;
//
//	Spline sx(x);
//	Spline sy(y);
//	TwoSpline s(sx, sy);
//
//
//	vd sol = s.resolver_spline();
//	imprimir_vector(sol);
//
//	vd resx,resy;
//
//	for(int i = 0; i < sol.size(); i++){
//		Punto p = s.evaluar(sol[i]);
//		resx.push_back(p.x);
//		resy.push_back(p.y);
//	}
//	imprimir_vector(resx);
//	//imprimir_vector(resy);
//}

void t_distancias() {
	for (int i = 1; i < t.size(); ++i) {
		cout << t[i] - t[i - 1] << " ";
	}
	cout << endl << endl;
}

void do_centripeta() {
	t.push_back(0);
	double acum = 0;
	for (int i = 1; i < cant_puntos; i++) {
		acum += sqrt(calc_distancia(i));
		t.push_back(acum);
	}
	Spline sx(x);
	Spline sy(y);
	Curva s(sx, sy);

	t_distancias();

	vp par_orig = s.resolver_param_original();
	s.imprimir_puntos(par_orig);
	vd sol = s.resolver_spline();

	s.imprimir_vector(sol);

}

double calc_distancia(int i) {
	double distx = (x[i] - x[i - 1]);
	distx *= distx;
	double disty = (y[i] - y[i - 1]);
	disty *= disty;
	return sqrt(distx + disty);
}

void imprimir_matrix(vvd mat) {
	int tam = mat.size();
	for (int i = 0; i < tam; i++) {
		for (int j = 0; j < tam; j++) {
			if (mat[i][j] > 99) {
				cout << setprecision(4);
				cout << mat[i][j] << "\t";
				cout << setprecision(precision);
			} else if ((mat[i][j] < 0.0001) && (mat[i][j] > -0.0001)) {
				cout << 0 << "\t";
			} else
				cout << mat[i][j] << "\t";
		}
		cout << endl;
	}
}

void Curva::imprimir_vector(vd vec) {
	int tam = vec.size();
	cout << tam << endl;
	for (int i = 0; i < tam; i++) {
		cout << sx.evaluar(vec[i]) << " " << sy.evaluar(vec[i]) << " "
				<< sx.derivar(vec[i]) << " " << sy.derivar(vec[i]) << " "
				<< endl;
	}
}

void imprimir_vector(vd vec) {
	int tam = vec.size();
	cout << "[";
	for (int i = 0; i < tam; i++) {
		cout << vec[i] << ",";
	}
	cout << "]" << endl;
}

void Curva::imprimir_puntos(vp ps) {
	int tam = ps.size();
	cout << tam << endl;
	;
	for (int i = 0; i < tam; i++) {
		cout << ps[i].x << " " << ps[i].y << " " << sx.derivar(delta * i) << " "
				<< sy.derivar(delta * i) << " " << endl;
	}
	cout << endl;
}

int main() {

	cout << setprecision(precision);

	cin >> cant_puntos;
	cin >> delta;

	double temp;

	forn(i, cant_puntos) {
		cin >> temp;
		x.push_back(temp);
		cin >> temp;
		y.push_back(temp);
	}

	do_centripeta();

	return 0;
}

