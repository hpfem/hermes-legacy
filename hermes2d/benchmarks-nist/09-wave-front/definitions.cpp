#include "definitions.h"

double CustomFunction::value(double x, double y) const
{
  double a = pow(x - x_loc, 2);
  double b = pow(y - y_loc, 2);
  double c = sqrt(a + b);
  double d = ((alpha * x - (alpha * x_loc)) * (2*x - (2 * x_loc)));
  double e = ((alpha * y - (alpha * y_loc)) * (2*y - (2 * y_loc)));
  double f = (pow(alpha * c - (alpha * r_zero), 2) + 1.0);
  double g = (alpha * c - (alpha * r_zero));

  return (((alpha/(c * f)) - (d/(2 * pow(a + b, 1.5) * f))
         - ((alpha * d * g)/((a + b) * pow(f, 2))) +
         (alpha/(c * f)) - (e/(2 * pow(a + b, 1.5) * f))
	    - ((alpha * e * g)/((a + b) * pow(f, 2)))));
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord(8);
}


double CustomExactSolution::value(double x, double y) const
{
  return atan(alpha * (sqrt(pow(x - x_loc, 2) + pow(y - y_loc, 2)) - r_zero));
}

void CustomExactSolution::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  double a = pow(x - x_loc, 2);
  double b = pow(y - y_loc, 2);
  double c = sqrt(a + b);
  double d = (alpha * x - (alpha * x_loc));
  double e = (alpha * y - (alpha * y_loc));
  double f = (pow(alpha * c - (alpha * r_zero), 2) + 1.0);

  dx = (d/(c * f));
  dy = (e/(c * f));
}

Ord CustomExactSolution::ord(Ord x, Ord y) const 
{
  return Ord(8);
}
