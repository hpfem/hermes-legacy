#include "definitions.h"

void CustomExactSolution::derivatives(double x, double y, scalar& dx, scalar& dy) const 
{
  dx = cos(x)*sin(y);
  dy = sin(x)*cos(y);
}

double CustomExactSolution::value(double x, double y) const 
{
  return sin(x)*sin(y);
}

Ord CustomExactSolution::ord(Ord x, Ord y) const 
{
  return Ord(7);
}


double CustomFunction::value(double x, double y) const 
{
  return -2*sin(x)*sin(y);
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord(7);
}
