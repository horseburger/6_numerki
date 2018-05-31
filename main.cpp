#include <iostream>
#include "gnuplot_i.hpp"
#include <cmath>

using namespace std;

// variables used for different variations of the problem
#define X_POINT 0
double w0,Beta,A,omega,phi; 

double F(double t);
double dY2(double x, double y1, double y2);
double rk4(double(*f)(double, double, double), double h, double x, double y1, double y2);
void plot(vector<double> x_axis, vector<double> y_axis);

int main(){

    double h = 0.01;
    double y1,y2;
    // y1 = x(t), y2 = v(t)

   /*  cout << "1. Oscylator prosty (bez tlumienia i wymuszenia\n" << "2. Oscylator slabo tlumiony (bez wymuszenia)\n" 
         << "3. Oscylator silnie tlumiony (bez wymuszenia)\n" << "4. Oscylator z wymuszeniem\n" << "5. Oscylator z wymuszeniem (rezonans)\n" 
         << "9. Exit\n" << "Twoj wybor: ";

    int a;
    cin >> a;

    do{
        cin >> a;
        switch(a){
            // x(0) = 1, v(0) = 0
            case 1:{
                w0 = 20*M_PI;
                Beta = 0;
                A = 0;
            }
            // x(0) = 1, v(0) = 0
            case 2:{
                w0 = 20*M_PI;
                Beta = 0.5;
                A = 0;
            }
            // x(0) = 1, v(0) = 0
            case 3:{
                w0 = 20*M_PI;
                Beta = 10;
                A = 0;
            }
            // x(0) = 0, v(0) = 0
            case 4:{
                w0 = 20*M_PI;
                Beta = 0.5;
                A = 50;
                omega = 2*M_PI;
                phi = 0;
            }
            // x(0) = 0, v(0) = 0
            case 5:{
                w0 = 20*M_PI;
                Beta = 0.5;
                A = 50;
                omega = 20*M_PI;
                phi = 0;
            }
        }
    }while(a!=9);
 */

    double tmp;
    y1 = 1;
    y2 = 0;
    w0 = 20*M_PI;
    Beta = 0;
    omega = 1;
    A = 0;

    tmp = rk4(dY2, h, X_POINT, y1, y2);
    cout << tmp;
    return 0;
}

// Funkcja f(x)
double F(double t){
    return A*sin(omega*t + phi);
}

// Druga pochodna f''(x)
double dY2(double x, double y1, double y2){
    double tmp;
    tmp = F(x) - w0*w0*y1 - 2*Beta*y2;
    cout << tmp;
    getchar();
    return tmp;
}


double rk4(double(*f)(double, double, double), double h, double x, double y1, double y2){
    double yi1;
    double k1,k2,k3,k4;
    k1 = h*f(x, y1, y2);
    k2 = h*f(x+h/2.0, y1+k1/2.0, y2+k1/2.0);
    k3 = h*f(x+h/2.0, y1+k2/2.0, y2+k2/2.0);
    k4 = h*f(x+h, y1+k3, y2+k3);
    yi1 = y1 + (k1+2.0*k2+2.0*k3+k4)/6.0;
    getchar();
    return yi1;
}

void plot(vector<double> x_axis, vector<double> y_axis){
    Gnuplot graph;
    graph.set_grid();
    graph.set_xrange(-10, 10);
    graph.set_yrange(-5, 5);
    graph.set_xlabel("Os x");
    graph.set_ylabel("Os y");
    graph.set_title("Metoda Rungego-Kutty4");
    graph.set_style("lines");

    graph.plot_xy(x_axis, y_axis, "Metoda Kungego-Kutty4");
}
