/*
#ifndef NUMERICALSOLVERS.H
#define NUMERICALSOLVERS.H
*/

#include <math.h>
#include <stdio.h>


/*Function definition to perform integration by Trapezoidal Rule */
double trapezoidal(double f(double x), double a, double b, int n){
  double x,h,sum=0,integral;
  int i;
  h=fabs(b-a)/n;
  for(i=1;i<n;i++){
    x=a+i*h;
    sum=sum+f(x);
  }
  integral=(h/2)*(f(a)+f(b)+2*sum);
  return integral;
}
/*Function definition to perform integration by Simpson's 1/3rd Rule */
double simpsons(double f(double x), double a,double b,double n){
  double h,integral,x,sum=0;
  int i;
  h=fabs(b-a)/n;
  for(i=1;i<n;i++){
    x=a+i*h;
    if(i%2==0){
      sum=sum+2*f(x);
    }
    else{
      sum=sum+4*f(x);
    }
  }
  integral=(h/3)*(f(a)+f(b)+sum);
  return integral;
}

  double euler(double f(double x, double y), double x0, double y0, double x, double h){
      double y;
      while(fabs(x-x0)>0.0000000001){
          y=y0+h*f(x0,y0);
          y0=y;
          x0=x0+h;
      }
      return y;
  }
  double RK1(double f(double x, double y), double x0, double y0, double x, double h){
      double y,k1,k2;
      while(fabs(x-x0)>0.0000000001){
          k1=h*f(x0,y0);
          k2=h*f(x0+h/2.0,y0+k1/2.0);
          y=y0+k2;
          y0=y;
          x0=x0+h;
      }
      return y;
  }

  double RK2(double f(double x, double y), double x0, double y0, double x, double h){
      double y,k1,k2,k3,k4;
      while(fabs(x-x0)>0.0000000001){
          k1=h*f(x0,y0);
          k2=h*f(x0+h/2.0,y0+k1/2.0);
          k3=h*f(x0+h/2.0,y0+k2/2.0);
          k4=h*f(x0+h,y0+k3);
          y=y0+1/6.0*(k1+2*k2+2*k3+k4);
          y0=y;
          x0=x0+h;
      }
      return y;
  }

//#endif /*NUMERICALSOLVERS.H */
