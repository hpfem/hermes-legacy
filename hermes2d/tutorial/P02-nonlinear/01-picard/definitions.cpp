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

/* Weak forms */

CustomWeakFormPicard::CustomWeakFormPicard(Solution* prev_iter_sln, HermesFunction* lambda, HermesFunction* f) 
                    : WeakForm(1)
{
  // Jacobian.
  CustomJacobian* matrix_form = new CustomJacobian(0, 0, lambda);
  matrix_form->ext.push_back(prev_iter_sln);
  add_matrix_form(matrix_form);

  // Residual.
  CustomResidual* vector_form = new CustomResidual(0, lambda, f);
  vector_form->ext.push_back(prev_iter_sln);
  add_vector_form(vector_form);
}

/* Essential boundary conditions */

double CustomEssentialBCNonConst::value(double x, double y, double n_x, double n_y, 
                       double t_x, double t_y) const
{
  return (x+10) * (y+10) / 100.;
}
