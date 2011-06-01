#include "definitions.h"

double CustomExactFunction::fn(double x, double y) 
{
  if (x <= 0) return cos(k * y);
  else return cos(k * y) + pow(x, alpha);
}


CustomFunction::CustomFunction(double k, double alpha) : HermesFunction(), k(k), alpha(alpha)
{
  cef = new CustomExactFunction(k, alpha);
}

CustomFunction::~CustomFunction()
{
  delete cef;
}

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


CustomExactSolution::CustomExactSolution(Mesh* mesh, double k, double alpha) : ExactSolutionScalar(mesh), k(k), alpha(alpha)
{
  cef = new CustomExactFunction(k, alpha);
}

CustomExactSolution::~CustomExactSolution() 
{ 
  delete cef; 
}

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
