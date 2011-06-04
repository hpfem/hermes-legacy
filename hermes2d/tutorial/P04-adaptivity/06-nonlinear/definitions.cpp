#include "definitions.h"

/* Initial condition */

scalar InitialSolutionHeatTransfer::value (double x, double y) const 
{
  return (x + 10) * (y + 10) / 100. + 2.;
}

void InitialSolutionHeatTransfer::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  dx = (y + 10) / 10.;
  dy = (x + 10) / 10.;
}

Ord InitialSolutionHeatTransfer::ord(Ord x, Ord y) const 
{
  return (x + 10) * (y + 10) / 100. + 2.;
}

/* Essential BC */
EssentialBoundaryCondition::EssentialBCValueType CustomEssentialBCNonConst::get_value_type() const 
{
  return EssentialBoundaryCondition::BC_FUNCTION;
}

scalar CustomEssentialBCNonConst::value(double x, double y, double n_x, double n_y, double t_x, double t_y) const 
{
  return (x + 10) * (y + 10) / 100.;
}
