#include "definitions.h"

/* Right-hand side */

double CustomExactFunction::fn(double x, double y) 
{
  if (x <= 0) return cos(k * y);
  else return cos(k * y) + pow(x, alpha);
}

/* Custom function f */

double CustomFunction::value(double x, double y) const 
{
  if (x < 0) return -cef->fn(x, y) * k * k;
  else return -(cef->fn(x, y) * k * k
              - alpha *(alpha - 1) * pow(x, alpha - 2.)
              - k * k * pow(x, alpha));
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord(20);
}

/* Exact solution */

double CustomExactSolution::value(double x, double y) const 
{
  return cef->fn(x, y);
}

void CustomExactSolution::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  if (x <= 0) dx = 0;
  else dx = alpha * pow(x, alpha - 1);
  dy = -sin(k * y) * k;
}

Ord CustomExactSolution::ord(Ord x, Ord y) const 
{
  return Ord(20);
}
