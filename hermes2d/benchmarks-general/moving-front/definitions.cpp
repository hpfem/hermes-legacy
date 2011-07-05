#include "definitions.h"

double CustomExactSolution::value (double x, double y) const 
{
  double t = *t_ptr;
  return 0;
}

void CustomExactSolution::derivatives(double x, double y, scalar& dx, scalar& dy) const 
{
  double t = *t_ptr;
  
  dx = 0;
  dy = 0;
}

Ord CustomExactSolution::ord(Ord x, Ord y) const 
{
  return Ord(20);
}

double CustomFunction::value(double x, double y) const 
{
  double S = s;
  double C = c;
  double t = *t_ptr;
  double f = -2*pow(S, 3)*pow(x, 2)*(-t + sqrt(pow(x, 2) + pow(y, 2)))*(x - x0)*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*(pow(x, 2) + pow(y, 2))*pow(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1, 2)) - 2*pow(S, 3)*pow(y, 2)*(-t + sqrt(pow(x, 2) + pow(y, 2)))*(x - x0)*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*(pow(x, 2) + pow(y, 2))*pow(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1, 2)) - S*pow(x, 2)*(x - x0)*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*pow(pow(x, 2) + pow(y, 2), 3.0/2.0)*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) + 2*S*x*(x - x0)*(y - y0)*(y - y1)*atan(t)/(C*sqrt(pow(x, 2) + pow(y, 2))*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) + 2*S*x*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*sqrt(pow(x, 2) + pow(y, 2))*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) - S*pow(y, 2)*(x - x0)*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*pow(pow(x, 2) + pow(y, 2), 3.0/2.0)*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) + 2*S*y*(x - x0)*(x - x1)*(y - y0)*atan(t)/(C*sqrt(pow(x, 2) + pow(y, 2))*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) + 2*S*y*(x - x0)*(x - x1)*(y - y1)*atan(t)/(C*sqrt(pow(x, 2) + pow(y, 2))*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) + S*(x - x0)*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) + 2*S*(x - x0)*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*sqrt(pow(x, 2) + pow(y, 2))*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) - 2*(x - x0)*(x - x1)*(-atan(S*(-t + sqrt(pow(x, 2) + pow(y, 2)))) + M_PI/2)*atan(t)/C - 2*(y - y0)*(y - y1)*(-atan(S*(-t + sqrt(pow(x, 2) + pow(y, 2)))) + M_PI/2)*atan(t)/C + (x - x0)*(x - x1)*(y - y0)*(y - y1)*(-atan(S*(-t + sqrt(pow(x, 2) + pow(y, 2)))) + M_PI/2)/(C*(pow(t, 2) + 1));

  return f;
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord(20);
}
