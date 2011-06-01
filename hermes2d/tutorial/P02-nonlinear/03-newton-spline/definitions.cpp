#include "definitions.h"

void CustomInitialCondition::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  dx = (y+10) / 100.;
  dy = (x+10) / 100.;
}

double CustomInitialCondition::value (double x, double y) const 
{
  return (x+10) * (y+10) / 100. + 2;
}

Ord CustomInitialCondition::ord(Ord x, Ord y) const 
{
  return x * y;
}

EssentialBoundaryCondition::EssentialBCValueType CustomEssentialBCNonConst::get_value_type() const 
{ 
  return EssentialBoundaryCondition::BC_FUNCTION; 
}

double CustomEssentialBCNonConst::value(double x, double y, double n_x, double n_y, 
                                        double t_x, double t_y) const
{
  return (x+10) * (y+10) / 100.;
}
