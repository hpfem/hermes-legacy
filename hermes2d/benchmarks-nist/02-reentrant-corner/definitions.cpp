#include "definitions.h"

double CustomExactSolution::value(double x, double y) const
{
  return (pow(sqrt(x*x + y*y), alpha) * sin(alpha * get_angle(y, x)));
};

void CustomExactSolution::derivatives (double x, double y, scalar& dx, scalar& dy) const
{
  double a = sqrt(x*x + y*y);
  double b = pow(a, (alpha - 1.0));
  double c = pow(a, alpha);
  double d = ((y*y)/(x*x) + 1.0 );

  dx = (((alpha* x* sin(alpha * get_angle(y,x)) *b)/a)
       - ((alpha *y *cos(alpha * get_angle(y, x)) * c)/(pow(x, 2.0) *d)));
  dy = (((alpha* cos(alpha* get_angle(y, x)) *c)/(x * d))
       + ((alpha* y* sin(alpha* get_angle(y, x)) *b)/a));
};

Ord CustomExactSolution::ord(Ord x, Ord y) const
{
  return Ord(10);
}

