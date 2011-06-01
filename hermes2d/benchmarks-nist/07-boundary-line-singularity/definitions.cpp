#include "definitions.h"

double CustomFunction::value(double x, double y) const
{
  return alpha * (alpha - 1.) * pow(x, alpha - 2.);
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord((int)(alpha + 3.1));
}


double CustomExactSolution::value(double x, double y) const
{
  return pow(x, alpha);
};

void CustomExactSolution::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  dx = alpha * pow(x, alpha - 1.);
  dy = 0;
};

Ord CustomExactSolution::ord(Ord x, Ord y) const
{
  return Ord((int)(alpha + 1));
};
