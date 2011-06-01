#include "definitions.h"

double CustomFunction::value(double x, double y) const
{
  double a_P = (-alpha * pow((x - x_loc), 2) - alpha * pow((y - y_loc), 2));

  return (4 * exp(a_P) * alpha * (alpha * (x - x_loc) * (x - x_loc)
         + alpha * (y - y_loc) * (y - y_loc) - 1));
}
 
Ord CustomFunction::value(Ord x, Ord y) const
{
  return Ord(8);
}

double CustomExactSolution::value(double x, double y) const
{
  return exp(-alpha * (pow((x - x_loc), 2) + pow((y - y_loc), 2)));
}

void CustomExactSolution::derivatives (double x, double y, scalar& dx, scalar& dy) const
{
  double a = -alpha * ( (x - x_loc) * (x - x_loc) + (y - y_loc) * (y - y_loc));
  dx = -exp(a) * (2 * alpha * (x - x_loc));
  dy = -exp(a) * (2 * alpha * (y - y_loc));
}

Ord CustomExactSolution::ord(Ord x, Ord y) const
{
  return Ord(8);
}
