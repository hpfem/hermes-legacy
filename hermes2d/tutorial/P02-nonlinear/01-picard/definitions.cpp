#include "definitions.h"

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
  return Ord(10);
}

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

scalar CustomWeakFormPicard::CustomJacobian::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
{
  scalar result = 0;
  for (int i = 0; i < n; i++) 
  {
    result += wt[i] * lambda->value(ext->fn[0]->val[i]) 
                    * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  }
  return result;
}

Ord CustomWeakFormPicard::CustomJacobian::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  Ord result = 0;
  for (int i = 0; i < n; i++) 
  {
    result += wt[i] * lambda->value(ext->fn[0]->val[i]) 
                    * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  }
  return result;
}

double CustomEssentialBCNonConst::value(double x, double y, double n_x, double n_y, 
                                        double t_x, double t_y) const
{
  return (x+10) * (y+10) / 100.;
}
