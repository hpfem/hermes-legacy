#include "definitions.h"

// Problem parameters.
double k_s = 20.464;
double alpha = 0.001;
double theta_r = 0;
double theta_s = 0.45;

// The pressure head is raised by H_OFFSET 
// so that the initial condition can be taken
// as the zero vector. Note: the resulting 
// pressure head will also be greater than the 
// true one by this offset.
double H_OFFSET = 1000;

scalar K(double h)
{
  if (h < 0) return k_s * exp(alpha * h);
  else return k_s;    
}

scalar dKdh(double h)
{
  if (h < 0) return k_s * alpha * exp(alpha * h);
  else return 0;
}

scalar ddKdhh(double h)
{
  if (h < 0) return k_s * alpha * alpha * exp(alpha * h);
  else return 0;
}

scalar C(double h)
{
  if (h < 0) return alpha * (theta_s - theta_r) * exp(alpha * h);
  else return alpha * (theta_s - theta_r);    
}

scalar dCdh(double h)
{
  if (h < 0) return alpha * (theta_s - theta_r) * alpha * exp(alpha * h);
  else return 0;    
}

scalar ddCdhh(double h)
{
  if (h < 0) return alpha * alpha * (theta_s - theta_r) * alpha * exp(alpha * h);
  else return 0;    
}

/* Custom non-constant Dirichlet condition */

EssentialBoundaryCondition::EssentialBCValueType CustomEssentialBCNonConst::get_value_type() const
{ 
  return EssentialBoundaryCondition::BC_FUNCTION; 
}

double CustomEssentialBCNonConst::value(double x, double y, double n_x, double n_y, 
                                        double t_x, double t_y) const 
{
  return x*(100. - x)/2.5 * y/100 - 1000. + H_OFFSET;
}

/* Custom weak forms */

CustomWeakFormRichardsIE::CustomWeakFormRichardsIE(double time_step, Solution* h_time_prev) : WeakForm(1)
{
  // Jacobian volumetric part.
  CustomJacobianFormVol* jac_form_vol = new CustomJacobianFormVol(0, 0, time_step);
  jac_form_vol->ext.push_back(h_time_prev);
  add_matrix_form(jac_form_vol);

  // Residual - volumetric.
  CustomResidualFormVol* res_form_vol = new CustomResidualFormVol(0, time_step);
  res_form_vol->ext.push_back(h_time_prev);
  add_vector_form(res_form_vol);
}

double CustomWeakFormRichardsIE::CustomJacobianFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                              Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
  {
    double h_val_i = h_prev_newton->val[i] - H_OFFSET;
    result += wt[i] * (   dCdh(h_val_i) * u->val[i] * (h_prev_newton->val[i] - h_prev_time->val[i]) 
                          * v->val[i] + C(h_val_i) * u->val[i] * v->val[i] 
			  + dKdh(h_val_i) * u->val[i] * (h_prev_newton->dx[i] * v->dx[i] 
							 + h_prev_newton->dy[i] * v->dy[i]) * time_step
                        + K(h_val_i) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) * time_step
			- ddKdhh(h_val_i) * u->val[i] * h_prev_newton->dy[i] * v->val[i] * time_step
                        - dKdh(h_val_i) * u->dy[i] * v->val[i] * time_step
		      );
  }
  return result;
}

Ord CustomWeakFormRichardsIE::CustomJacobianFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                                                         Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(20);
}

WeakForm::MatrixFormVol* CustomWeakFormRichardsIE::CustomJacobianFormVol::clone() 
{
  return new CustomJacobianFormVol(*this);
}

scalar CustomWeakFormRichardsIE::CustomResidualFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                                                              ExtData<scalar> *ext) const 
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
  {
    double h_val_i = h_prev_newton->val[i] - H_OFFSET;
    result += wt[i] * (   C(h_val_i) * (h_val_i - (h_prev_time->val[i] - H_OFFSET)) * v->val[i]
                        + K(h_val_i) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i]) * time_step
                        - dKdh(h_val_i) * h_prev_newton->dy[i] * v->val[i] * time_step
                      );
  }
  return result;
}

Ord CustomWeakFormRichardsIE::CustomResidualFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(20);
}

WeakForm::VectorFormVol* CustomWeakFormRichardsIE::CustomResidualFormVol::clone() 
{
  return new CustomResidualFormVol(*this);
}


