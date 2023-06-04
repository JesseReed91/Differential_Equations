/*
Copyright (C) Jesse Reed 2023, see end of file for extended copyright information

This program was written to find numerical approximations differential equations without analytical solutions,
the methods available are Euler's Method, Heun's Method (Improved Euler's Method), and the Runge-Kutta Method (RK4).
This was originally written as a side project for my differential equations class to compare the accuracy of the different methods for numerical approximations.
There are four required inputs from the user, the differential equation (Fxy), the initial conditions (xInitial, yInitial), the step size (h), 
and the y value we are approximating to (TARGETVALUE). These must be manipulated by the user in the .cpp file itself
*/

//TODO: implement a menu that allows a user to input the values instead of having to manipulate the .cpp file and add documentation for functions

#include <iostream>
#include <cmath>

using namespace std;

double Fxy(double x, double y); //the differential equation
double Xn(double x, double h); //calculate Xn+1
double eulerYn(double x, double y, double h); //calculate Yn+1 for Euler's method
double heunYn(double x, double y, double h); //calculate Yn+1 for Heun's method
double rkYn(double x, double y, double h); //calculate Yn+1 for Runge-Kutta method
double k1(double x, double y, double h); //calculate k1 for Runge-Kutta method
double k2(double x, double y, double h); //calculate k2 for Runge-Kutta method
double k3(double x, double y, double h); //calculate k3 for Runge-Kutta method
double k4(double x, double y, double h); //calculate k4 for Runge-Kutta method
double sumK(double x, double y, double h); //calculate (k1 + 2k2 + 2k3 + k4) * (1 / 6) for Runge-Kutta method
void eulersApprox(double x, double y, double h, int n); //recursive method to find Euler's Method for the numerical solution to an ODE
void heunApprox(double x, double y, double h, int n); //recusrive method to find Heun's Method for the numerical solution to an ODE
void rungeKuttaApprox(double x, double y, double h, int n); //recursive method to find the Runge Kutta approximation for the numerical solution to the ODE


double const EPSILON = 0.0000001; //this is required to evaluate equalities and inequalities for doubles because computers can only approximate real numbers
double const TARGETVALUE = 0.5; //the y value that you'd like to approximate to

int main()
{
    double xInitial = 0; //intial x value
    double yInitial = 0.5; //initial y value
    double const h = 0.1; //step size
    int n = 0; //counter for each step

    std::cout << "Euler's Method" << endl;
    eulersApprox(xInitial, yInitial, h, n);
    std::cout << endl;
    std::cout << "Heun's Method" << endl;
    heunApprox(xInitial, yInitial, h, n);
    std::cout << endl;
    std::cout << "Runge-Kutta Method" << endl;
    rungeKuttaApprox(xInitial, yInitial, h, n);
}

double Fxy(double x, double y)
{
    return (x - y) * (x - y);
}

double Xn(double x, double h)
{
    return x + h;
}

double eulerYn(double x, double y, double h)
{
    return y + (h * Fxy(x, y));
}

double heunYn(double x, double y, double h)
{
    return y + h * ((Fxy(x, y) + Fxy(Xn(x, h), eulerYn(x, y, h))) / 2);
}

double rkYn(double x, double y, double h)
{
    return y + sumK(x, y, h);
}

double k1(double x, double y, double h)
{
    return h * Fxy(x, y);
}

double k2(double x, double y, double h)
{
    return (h * Fxy((x + ((0.5) * h)), (y + ((0.5) * k1(x, y, h)))));
}

double k3(double x, double y, double h)
{
    return (h * Fxy((x + ((0.5) * h)), (y + ((0.5) * k2(x, y, h)))));
}

double k4(double x, double y, double h)
{
    return (h * Fxy((x + h), (y + k3(x, y, h))));
}

double sumK(double x, double y, double h)
{
    //this is done this way to improve readability
    double oneSixth = static_cast<double>(1) / 6;
    double sum = k1(x, y, h);
    sum += 2 * k2(x, y, h);
    sum += 2 * k3(x, y, h);
    sum += k4(x, y, h);
    sum = oneSixth * sum;
    return sum;
}

void eulersApprox(double x, double y, double h, int n)
{
    if (x - TARGETVALUE < EPSILON)
    {
        std::cout << n << ": " << "x: " << x << ", y: " << y << endl;
        y = eulerYn(x, y, h);
        x = Xn(x, h);
        n++;
        eulersApprox(x, y, h, n);
    }

}

void heunApprox(double x, double y, double h, int n)
{
    if (x - TARGETVALUE < EPSILON)
    {
        std::cout << n << ": " << "x: " << x << ", y: " << y << endl;
        y = heunYn(x, y, h);
        x = Xn(x, h);
        n++;
        heunApprox(x, y, h, n);
    }

}

void rungeKuttaApprox(double x, double y, double h, int n)
{
    if (x - TARGETVALUE < EPSILON)
    {
        std::cout << n << ": " << "x: " << x << ", y: " << y << endl;
        y = rkYn(x, y, h);
        x = Xn(x, h);
        n++;
        rungeKuttaApprox(x, y, h, n);
    }
}


/*
Copyright (C) Jesse Reed

This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License or (at your option) any later version.

This is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see: https://www.gnu.org/licenses
*/