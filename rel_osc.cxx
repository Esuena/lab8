// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx, double* k1, double* k2, double* k3, double* k4);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
  const int dim = 2;
	double dx = 0.001,x=0;
	const double L = 100;
  for(double p0 = 0.001; p0 < 5; p0 += 0.1){  //Amplitude des Oszi wird verändert
	
	double y0[dim] = {p0, 0};//p0 muss bei Maximum starten, d.h. {1.0, 0}, damit der faellt: Aenderung der Anfangsbedingungen
	double yn[dim];
        double k1[dim], k2[dim], k3[dim], k4[dim];
//   out << x << "\t" << y0[0] << "\t" << y0[1] << "\t"  << endl;
	x = 0;
	while(x<=L)
	{
		x += dx;
		RKstep(yn, y0, x, dx, k1, k2 , k3, k4);
		
		  if(yn[1]<0 && y0[1]>0)
		    break;
		  
                for(int i=0; i<dim; i++) y0[i] = yn[i];
// 		out << x << "\t" << y0[0] << "\t" << y0[1] << "\t" << endl;
	}
	cout << x - dx << endl;
	double theta;
	double theta_l=0;
	double theta_r=1;
	
	while(abs(yn[1])>1e-8){
	  theta = (theta_l+theta_r)/2.;
	  
	  double b1, b2, b3, b4;
	
	  b1 = theta - 3*theta*theta/2 + 2*theta*theta*theta/3;
	  b2 = b3 = theta*theta - 2*theta*theta*theta/3;
	  b4 = -theta*theta/2 + 2*theta*theta*theta/3;

	  yn[1] = y0[1] + dx*(b1*k1[1]+b2*k2[1]+b3*k3[1]+b4*k4[1]);
	  
	  if(yn[1]>0)
	    theta_l = theta;
	  else
	    theta_r = theta;
	  
// 	  cout <<theta_r <<"\t"<<theta_l <<"\t" << yn[1] << endl;
	  
	}
	
	out <<p0 << "\t" << x -dx + theta*dx << endl; 
	
  }
	out.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx,
	    double* k1, double* k2, double* k3, double* k4)
{
	const int dim = 2;//dimension muss geaendert werden, da dimension schon in main fct. festgelegt wurde fuer k1, k2, k3 und k4
	//danger of overwriting !!!
	//Alternative waere die Reihenfolge (s.oben) zu aendern, d.h. k4[dim], k3[dim], usw.

  for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
  f(k2, x+0.5*dx);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3, x+0.5*dx);

  for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4,  x+dx);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}
//-------------------
//Rel. Oszi
void f(double* const y0, const double x)
{
// 	const double a = 10;
// 	const double b = 28;
// 	const double c = 8.0/3.0;
	double y[2] = { y0[0], y0[1] };

        y0[0] = y[1];
	y0[1] = -y[0]/sqrt(1.+y[0]*y[0]);
// 	y0[2] = y[0]*y[1] - c*y[2];
}
