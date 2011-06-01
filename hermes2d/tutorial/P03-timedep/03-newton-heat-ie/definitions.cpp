#include "definitions.h"

CustomNonlinearity::CustomNonlinearity(double alpha): HermesFunction()
{
  this->is_const = false;
  this->alpha = alpha;
}

scalar CustomNonlinearity::value(double u) const
{
  return 1 + pow(u, alpha);
}

Ord CustomNonlinearity::value(Ord u) const
{
  return Ord(10);
}

scalar CustomNonlinearity::derivative(double u) const
{
  return alpha * pow(u, alpha - 1.0);
}

Ord CustomNonlinearity::derivative(Ord u) const
{
  return Ord(10);
}


CustomWeakForm::CustomWeakForm(HermesFunction* lambda, HermesFunction* f, double tau, Solution* sln_prev_time) 
      : WeakForm(1) 
{
  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 0, HERMES_ANY, new HermesFunction(1. / tau)));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, HERMES_ANY, lambda));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualVol(0, HERMES_ANY, new HermesFunction(1. / tau)));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, HERMES_ANY, lambda));
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol(0, HERMES_ANY, f));
  CustomVectorFormVol* vector_form = new CustomVectorFormVol(0, tau);
  vector_form->ext.push_back(sln_prev_time);
  add_vector_form(vector_form);
}

scalar CustomWeakForm::CustomVectorFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                                  Geom<double> *e, ExtData<scalar> *ext) const 
{
  scalar result = 0;
  Func<scalar>* u_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result -= wt[i] * u_prev_time->val[i] * v->val[i] / tau;
  return result;
}

Ord CustomWeakForm::CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                             Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  Ord result = 0;
  Func<Ord>* u_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result -= wt[i] * u_prev_time->val[i] * v->val[i] / tau;
  return result;
}


EssentialBCNonConst::EssentialBCNonConst(std::string marker) : EssentialBoundaryCondition(Hermes::vector<std::string>())
{
  markers.push_back(marker);
}

EssentialBoundaryCondition::EssentialBCValueType EssentialBCNonConst::get_value_type() const
{ 
  return EssentialBoundaryCondition::BC_FUNCTION; 
}

scalar EssentialBCNonConst::value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
{
  return (x+10)*(y+10)/100.;
}


void CustomInitialCondition::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  dx = (y+10)/100.;
  dy = (x+10)/100.;
}

scalar CustomInitialCondition::value (double x, double y) const 
{
  return (x+10)*(y+10)/100.;
}

Ord CustomInitialCondition::ord(Ord x, Ord y) const 
{
  return (x+10)*(y+10)/100.;
}
