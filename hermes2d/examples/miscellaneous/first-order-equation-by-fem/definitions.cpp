#include "definitions.h"

/* Weak forms */

CustomWeakForm::CustomWeakForm() : WeakForm(1)
{
  // Jacobian form.
  add_matrix_form(new CustomJacobian());

  // Residual form.
  add_vector_form(new CustomResidual());
};

scalar CustomWeakForm::CustomJacobian::value(int n, double *wt, Func<scalar> *u_ext[], 
                                             Func<double> *u, Func<double> *v,
                                             Geom<double> *e, ExtData<scalar> *ext) const
{
  scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * (u->dx[i] - u->dy[i]) * v->val[i];
    //result += wt[i] * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]);
  }
  return result;
}
  
Ord CustomWeakForm::CustomJacobian::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                        Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  Ord result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * (u->dx[i] - u->dy[i]) * v->val[i];
    //result += wt[i] * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]);
  }
  return result;
}

scalar CustomWeakForm::CustomResidual::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                             Geom<double> *e, ExtData<scalar> *ext) const
{
  scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * (u_ext[0]->dx[i] - u_ext[0]->dy[i]) * v->val[i];
    //result += wt[i] * (u_ext[0]->dx[i]*v->dx[i] + u_ext[0]->dy[i]*v->dy[i]);
  }
  return result;
}

Ord CustomWeakForm::CustomResidual::ord(int n, double *wt, Func<Ord> *u_ext[],
                                        Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  Ord result = 0;
  for (int i = 0; i < n; i++) {
    result += wt[i] * (u_ext[0]->dx[i] - u_ext[0]->dy[i]) * v->val[i];
    //result += wt[i] * (u_ext[0]->dx[i]*v->dx[i] + u_ext[0]->dy[i]*v->dy[i]);
  }
  return result;
}

/* Custom non-constant Dirichlet condition */

CustomDirichletCondition::CustomDirichletCondition(Hermes::vector<std::string> markers, 
                                                   double A, double B, double C)
  : EssentialBoundaryCondition(markers), A(A), B(B), C(C) 
{ 
}

EssentialBoundaryCondition::EssentialBCValueType CustomDirichletCondition::get_value_type() const
{ 
  return EssentialBoundaryCondition::BC_FUNCTION; 
}

scalar CustomDirichletCondition::value(double x, double y, double n_x, double n_y, 
                                       double t_x, double t_y) const 
{
  return A*x + B*y + C;
}
