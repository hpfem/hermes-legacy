#include "definitions.h"

void CustomExactSolution::derivatives(double x, double y, scalar& dx, scalar& dy) const 
{
  dx = cos(x);
  dy = 0;
}

double CustomExactSolution::value(double x, double y) const 
{
  return sin(x);
}

Ord CustomExactSolution::ord(Ord x, Ord y) const 
{
  return Ord(7);
}


double CustomFunction::value(double x, double y) const 
{
  return -sin(x);
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord(7);
}


CustomWeakFormPoisson::CustomWeakFormPoisson(std::string bdy_marker_right) : WeakForm(1) 
{
  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0));
   // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0));
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol(0, HERMES_ANY, new CustomFunction));
  add_vector_form_surf(new WeakFormsH1::DefaultVectorFormSurf(0, bdy_marker_right, new HermesFunction(1.0)));
}
