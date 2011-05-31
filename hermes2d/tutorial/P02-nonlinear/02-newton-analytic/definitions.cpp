#include "definitions.h"

/* Nonlinearity lambda(u) = pow(u, alpha) */

CustomNonlinearity::CustomNonlinearity(double alpha): HermesFunction()
{
  this->is_const = false;
  this->alpha = alpha;
}

double CustomNonlinearity::value(double u) const
{
  return 1 + pow(u, alpha);
}

Ord CustomNonlinearity::value(Ord u) const
{
  // If alpha is not an integer, then the function
  // is non-polynomial. 
  // NOTE: Setting Ord to 10 is safe but costly,
  // one could save here by looking at special cases 
  // of alpha. 
  return Ord(10);
}

double CustomNonlinearity::derivative(double u) const
{
  return alpha * pow(u, alpha - 1.0);
}

Ord CustomNonlinearity::derivative(Ord u) const
{
  // Same comment as above applies.
  return Ord(10);
}

/* Initial condition for the Newton's method */

double CustomInitialCondition::value(double x, double y) const 
{
  return (x+10) * (y+10) / 100. + 2;
}

void CustomInitialCondition::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{   
  dx = (y+10) / 100.;
  dy = (x+10) / 100.;
}

Ord CustomInitialCondition::ord(Ord x, Ord y) const 
{
  return x*y;
}

/* Essential boundary conditions */

double CustomEssentialBCNonConst::value(double x, double y, double n_x, double n_y, 
                                        double t_x, double t_y) const
{
  return (x+10) * (y+10) / 100.;
}
