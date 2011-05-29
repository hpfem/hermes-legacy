#include "definitions.h"

/* Custom function f */

double CustomFunction::value(double x, double y) const 
{
  return (-sin(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),4)
           + 2*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),2)*pow((pow(x,2) + pow(y,2)),(1.0/2.0)))
           + pow(x,2)*sin(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),4)*(pow(x,2) + pow(y,2)))
           + pow(y,2)*sin(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),4)*(pow(x,2) + pow(y,2)))
           - pow(x,2)*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),2)*pow((pow(x,2) + pow(y,2)),(3.0/2.0)))
           - pow(y,2)*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),2)*pow((pow(x,2) + pow(y,2)),(3.0/2.0)))
           - 2*pow(x,2)*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),3)*(pow(x,2) + pow(y,2)))
           - 2*pow(y,2)*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),3)*(pow(x,2) + pow(y,2))));
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord(10);
}

/* Exact solution */

double CustomExactSolution::value(double x, double y) const 
{
  double r = sqrt(x*x + y*y);
  return sin(1/(alpha + r));
}

void CustomExactSolution::derivatives(double x, double y, scalar& dx, scalar& dy) const 
{
  double r = sqrt(x*x + y*y);
  double h = 1/(alpha + r);
  dx = -cos(h) * h * h * x / r;
  dy = -cos(h) * h * h * y / r;
}

Ord CustomExactSolution::ord(Ord x, Ord y) const 
{
  return Ord(10);
}

/* Weak forms */

CustomWeakForm::CustomWeakForm(CustomFunction* f) : WeakForm(1) 
{
  // Jacobian.
  add_matrix_form(new CustomMatrixFormVol(0, 0, f->alpha));

  // Residual.
  add_vector_form(new CustomVectorFormVol(0, f));
}
