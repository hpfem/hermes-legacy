#include "definitions.h"

// Problem parameters.
double k_s = 20.464;
double alpha = 0.001;
double theta_r = 0;
double theta_s = 0.45;

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
  return x*(100. - x)/2.5 * pow(y/100., y_power) - 1000.;
}

/* Custom Initial condition */

void CustomInitialCondition::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  dx = (100. - 2*x)/2.5 * pow(y/100., y_power);
  dy = x*(100. - x)/2.5 * pow(y/100., y_power - 1) * 1./100.;
}

double CustomInitialCondition::value (double x, double y) const 
{
  return x*(100. - x)/2.5 * pow(y/100., y_power) - 1000.;
}

Ord CustomInitialCondition::ord(Ord x, Ord y) const 
{
  return Ord(10);
}

CustomWeakFormRichardsIE::CustomWeakFormRichardsIE(double time_step, Solution* u_time_prev) : WeakForm(1)
{
  // Jacobian volumetric part.
  add_matrix_form(new CustomFormMatrixFormVol(0, 0, time_step));

  // Residual - volumetric.
  add_vector_form(new CustomResidualFormVol(0, time_step));
  CustomVectorFormVol* vec_form_vol = new CustomVectorFormVol(0, time_step);
  vec_form_vol->ext.push_back(u_time_prev);
  add_vector_form(vec_form_vol);
}


template<typename Real, typename Scalar>
Scalar CustomWeakFormRichardsIE::CustomVectorFormVol::vector_form_ie(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                                                     Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Func<double>* temp_prev_time = ext->fn[0];
  return -C(temp_prev_time->val[i]) * int_u_v<double, scalar>(n, wt, temp_prev_time, v) / time_step;
}

scalar CustomWeakFormRichardsIE::CustomVectorFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                                                            ExtData<scalar> *ext) const 
{
  return vector_form_ie<double, scalar>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormRichardsIE::CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);
}

WeakForm::VectorFormVol* CustomWeakFormRichardsIE::CustomVectorFormVol::clone() 
{
  return new CustomVectorFormVol(*this);
}


template<typename Real, typename Scalar>
Scalar CustomWeakFormRichardsIE::CustomFormMatrixFormVol::matrix_form_ie(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * ( C(u->val[i]) * u->val[i] * v->val[i] / time_step
                      + K(u->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                      - dKdh(u->val[i]) * u->dy[i] * v->val[i]);
  }
  return result;
}

scalar CustomWeakFormRichardsIE::CustomFormMatrixFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                                                                Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return matrix_form_ie<double, scalar>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakFormRichardsIE::CustomFormMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                                                           Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);
}

WeakForm::MatrixFormVol* CustomWeakFormRichardsIE::CustomFormMatrixFormVol::clone() 
{
  return new CustomFormMatrixFormVol(*this);
}


template<typename Real, typename Scalar>
Scalar CustomWeakFormRichardsIE::CustomResidualFormVol::vector_form_ie(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                                                       Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  Func<Scalar>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * ( C(h_prev_newton->val[i]) * h_prev_newton->val[i] * v->val[i] / time_step
                      + K(h_prev_newton->val[i]) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                      - dKdh(h_prev_newton->val[i]) * h_prev_newton->dy[i] * v->val[i]);
  }
  return result;
}

scalar CustomWeakFormRichardsIE::CustomResidualFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                                                              ExtData<scalar> *ext) const 
{
  return vector_form_ie<double, scalar>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormRichardsIE::CustomResidualFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);
}

WeakForm::VectorFormVol* CustomWeakFormRichardsIE::CustomResidualFormVol::clone() 
{
  return new CustomResidualFormVol(*this);
}


