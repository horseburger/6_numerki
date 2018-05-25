#include <iostream>
#include "gnuplot_i.hpp"
#include <cmath>

using namespace std;

double Y1(double t, double A, double omega, double phi);
double Y2(double t, double A, double omega, double phi);
void plot();
double rk(double(*f)(double, double), double h, double x, double y);

int main(){

    double w0,beta,A,omega,phi;

    double rk1[4], rk2[4];
    double h = 0.5;


    cout << "1. Oscylator prosty (bez tlumienia i wymuszenia\n" << "2. Oscylator slabo tlumiony (bez wymuszenia)\n" 
         << "3. Oscylator silnie tlumiony (bez wymuszenia)\n" << "4. Oscylator z wymuszeniem\n" << "5. Oscylator z wymuszeniem (rezonans)\n" 
         << "9. Exit\n" << "Twoj wybor: ";

    int a;
    cin >> a;

    do{
        cin >> a;
        switch(a){
            case 1:{
                w0 = 20*M_PI;
                beta = 0;
                A = 0;
                rk1[0] = 1;
            }
            case 2:{
                w0 = 20*M_PI;
                beta = 0.5;
                A = 0;
                rk1[0] = 1;
            }
            case 3:{
                w0 = 20*M_PI;
                beta = 10;
                A = 0;
                rk1[0] = 1;
            }
            case 4:{
                w0 = 20*M_PI;
                beta = 0.5;
                A = 50;
                omega = 2*M_PI;
                phi = 0;
                rk1[0] = 0;
            }
            case 5:{
                w0 = 20*M_PI;
                beta = 0.5;
                A = 50;
                omega = 20*M_PI;
                phi = 0;
                rk1[0] = 0;
            }
        }
    }while(a!=9);

    return 0;
}

double Y1(double t, double A, double omega, double phi){
    return A*sin(omega*t + phi);
}

double Y2(double t, double A, double omega, double phi){
    return A*cos(omega*t + phi)*omega;
}


void plot();

double k(double h, double rk1[], double rk2[]){
    rk1[1] = h*rk2[0];
    rk2[1] = h*
} 