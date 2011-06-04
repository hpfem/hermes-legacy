#include "definitions.h"

scalar ExactSolutionCustom::value (double x, double y) const 
{
  return pow(x*x + y*y, 0.25);
}

void ExactSolutionCustom::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  dx = 0.25 * pow(x*x + y*y, -0.75) * 2 * x;
  dy = 0.25 * pow(x*x + y*y, -0.75) * 2 * y;
}

Ord ExactSolutionCustom::ord(Ord x, Ord y) const 
{
  return pow(x*x + y*y, 0.25);
}
