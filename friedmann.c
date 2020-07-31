#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numericalsolvers.h"
#include "cosmological_const.h"

double friedmann_mass(double a,double t){
    return H*(sqrt(Om/a ));
}

double friedmann_rad(double a, double t){
  return H*(sqrt(Or/(a*a)));

}

double friedmann_de(double a, double t ){
  return H*(sqrt(Ol*(a*a)));
}

double friedmann_curvature(double a, double t){
  return H*(sqrt(Ok));
}

double friedmann_matter_de(double a,double t){
  return H*(sqrt( Ol*(a*a) + Om/a));

}

double friedmann(double a, double t){
  return H*(sqrt(Or/(a*a) + Om/a + Ol*(a*a)));
}

int main()
{
  FILE *fp, *fp1, *fp2, *fp3, *fp4, *fp5;
  double sol_friedmman, t;
  // Make the files to write the solutions
  fp = fopen("friedmann_matter_single_k=0.dat","w");
  fp1 = fopen("friedmann_rad_single_k=0.dat","w");
  fp2 = fopen("friedmann_de_single_k=0.dat", "w");
  fp3 = fopen("friedmann_curv_single.dat","w");
  fp4 = fopen("friedmann_two_matter_de.dat","w");
  fp5 = fopen("friedmann_complete.dat","w");
  /* The initial condition is:
      a(13) = 1
  */
  double a0 = 1 , t0 = 13;

  for(t= 1; t <= 13; t = t + 0.01){
    sol_friedmman = RK1(friedmann_mass,a0,t0, t , 0.005);
    fprintf(fp,"%lf\t%lf\n",t,sol_friedmman);

    sol_friedmman = RK1(friedmann_rad,a0, t0 , t, 0.005);
    fprintf(fp1,"%lf\t%lf\n",t,sol_friedmman);

    sol_friedmman = RK1(friedmann_de,a0, t0 , t, 0.005);
    fprintf(fp2,"%lf\t%lf\n",t,sol_friedmman);

    sol_friedmman = RK1(friedmann_curvature,a0, t0 , t, 0.005);
    fprintf(fp3,"%lf\t%lf\n",t,sol_friedmman);

    sol_friedmman = RK1(friedmann_matter_de,a0, t0 , t, 0.005);
    fprintf(fp4,"%lf\t%lf\n",t,sol_friedmman);

    sol_friedmman = RK1(friedmann,a0, t0 , t, 0.005);
    fprintf(fp5,"%lf\t%lf\n",t,sol_friedmman);

    }
  }
