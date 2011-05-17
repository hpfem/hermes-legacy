#include "definitions.h"


double CustomExactSolution::value(double x, double y) const
{
  
  return (pow(sqrt(x*x + y*y), ALPHA) * sin(ALPHA * get_angle(y, x)));
  };

  void CustomExactSolution::derivatives (double x, double y, scalar& dx, scalar& dy) const
  {
    double a = sqrt(x*x + y*y);
    double b = pow(a, (ALPHA - 1.0));
    double c = pow(a, ALPHA);
    double d = ((y*y)/(x*x) + 1.0 );

    dx = (((ALPHA* x* sin(ALPHA * get_angle(y,x)) *b)/a)
         - ((ALPHA *y *cos(ALPHA * get_angle(y, x)) * c)/(pow(x, 2.0) *d)));
    dy = (((ALPHA* cos(ALPHA* get_angle(y, x)) *c)/(x * d))
         + ((ALPHA* y* sin(ALPHA* get_angle(y, x)) *b)/a));
  };

  Ord CustomExactSolution::ord(Ord x, Ord y) const
  {
    return Ord(10);
  }

