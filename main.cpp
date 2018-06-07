#include <iostream>
#include "gnuplot_i.hpp"
#include <cmath>

using namespace std;

// variables used for different variations of the problem
double w0,Beta,A,omega,phi; 

double rk4(double(*f)(double, double, double), double h, double x, double y1, double y2);
void plot(vector<double> x_axis, vector<double> y_axis, double n);
double pattern1(double x, double y1, double y2);
double pattern2(double x, double y1, double y2);
void k(double k1[], double k2[],double x1, double y1, double y2, double h);
double RungeKutta(double k[], double h);
void yi(double x, double y1, double y2, double h, double sol1[], double sol2[], int j);

int main(){

    double h = 0.01;
    double y1,y2;
    double x1 = 0;
    int n = 100;
    y1 = 1;
    y2 = 0;
    // y1 = x(t), y2 = v(t)

    cout << "1. Oscylator prosty (bez tlumienia i wymuszenia\n" << "2. Oscylator slabo tlumiony (bez wymuszenia)\n" 
         << "3. Oscylator silnie tlumiony (bez wymuszenia)\n" << "4. Oscylator z wymuszeniem\n" << "5. Oscylator z wymuszeniem (rezonans)\n" 
         << "9. Exit\n" << "Twoj wybor: ";

    int a;
    cin >> a;
    cin.get();

    switch(a){
        // x(0) = 1, v(0) = 0
        case 1:{
            w0 = 20*M_PI;
            Beta = 0;
            A = 0;
            break;
        }
        // x(0) = 1, v(0) = 0
        case 2:{
            w0 = 20*M_PI;
            Beta = 0.5;
            A = 0;
            break;
        }
        // x(0) = 1, v(0) = 0
        case 3:{
            w0 = 20*M_PI;
            Beta = 10;
            A = 0;
            break;
        }
        // x(0) = 0, v(0) = 0
        case 4:{
            w0 = 20*M_PI;
            Beta = 0.5;
            A = 50;
            omega = 2*M_PI;
            phi = 0;
            break;
        }
        // x(0) = 0, v(0) = 0
        case 5:{
            w0 = 20*M_PI;
            Beta = 0.5;
            A = 50;
            omega = 20*M_PI;
            phi = 0;
            break;
        }
    }





    double * rk1 = new double[n + 1];
	double * rk2 = new double[n + 1];
	double * rk3 = new double[n + 1];

	rk1[0] = y1;
	rk2[0] = y2;
	rk3[0] = x1;

	for (int j = 1; j <= n; j++)
	{
		yi(x1, rk1[j-1], rk2[j-1], h, rk1, rk2, j);
		x1 = x1 + h;
		rk3[j] = x1;
	}

    for (int i = 0; i < n; i++)
	{
		cout << endl << "Wartosc dla x = " << rk3[i] << endl <<  "Wyniki funkcji pierwszej: "  << rk1[i] << "   "  << "drugiej: " << rk2[i] << endl ;
	}

    vector<double> x_axis(n);
    vector<double> y_axis(n);
    double tmp = x1;
    for(int i=0;i<n;i++){
        x_axis[i] = rk3[i];
    }
    for(int i=0;i<n;i++){
        y_axis[i] = rk2[i];
    }

    plot(x_axis,y_axis, n);

    return 0;
}

// Funkcja f(x)
double F(double t){
    return A*sin(omega*t + phi);
}



/* double rk4(double(*f)(double, double, double), double h, double x, double y1, double y2){
    double yi1;
    double k1,k2,k3,k4;
    k1 = h*f(x, y1, y2);
    k2 = h*f(x+h/2.0, y1+k1/2.0, y2+k1/2.0);
    k3 = h*f(x+h/2.0, y1+k2/2.0, y2+k2/2.0);
    k4 = h*f(x+h, y1+k3, y2+k3);
    yi1 = y1 + (k1+2.0*k2+2.0*k3+k4)/6.0;
    return yi1;
} */



double pattern1(double x, double y1, double y2)
{
	return A*sin(omega*x + phi);
}

double pattern2(double x, double y1, double y2)
{
	double tmp;
    tmp = F(x) - w0*w0*y1 - 2*Beta*y2;
    return tmp;
}

void k(double k1[], double k2[],double x1, double y1, double y2, double h)
{
	cout << endl;	

	k1[0] = h*(pattern1(x1, y1, y2));
	k2[0] = h*(pattern2(x1, y1, y2));

	k1[1] = h*(pattern1((x1 + (h / 2.0)), (y1 + (k1[0] / 2.0)), (y2 + (k2[0] / 2.0))));
	k2[1] = h*(pattern2((x1 + (h / 2.0)), (y1 + (k1[0] / 2.0)), (y2 + (k2[0] / 2.0))));

	k1[2] = h*(pattern1((x1 + h), (y1 + (2.0 * k1[1]) - k1[0]), (y2 + (2.0 * k2[1]) - k2[0])));
	k2[2] = h*(pattern2((x1 + h), (y1 + (2.0 * k1[1]) - k1[0]), (y2 + (2.0 * k2[1]) - k2[0])));

}

double RungeKutta(double k[], double h)
{
	double yi1;
	double a = k[0], b = (4.0 * k[1]), c = k[2];
	double d = (a + b + c);

	yi1 = d/6.0;
	return yi1;
}


void yi(double x, double y1, double y2, double h, double sol1[], double sol2[], int j)
{
	double y3, y4;
	double k2[3], k3[3];

	k(k2, k3, x, y1, y2, h);

	y3 = RungeKutta(k2, h);
	y4 = RungeKutta(k3, h);
	sol1[j] = y1 + y3;
	sol2[j] = y2 + y4;
}

void plot(vector<double> x_axis, vector<double> y_axis, double n){
    Gnuplot graph;
    graph.set_grid();
    graph.set_xrange(-10, 10);
    graph.set_yrange(y_axis[0], y_axis[n-1]);
    graph.set_xlabel("Os x");
    graph.set_ylabel("Os y");
    graph.set_title("Metoda Rungego-Kutty4");
    graph.set_style("lines");

    graph.plot_xy(x_axis, y_axis, "Metoda Kungego-Kutty4");
    getchar();
}