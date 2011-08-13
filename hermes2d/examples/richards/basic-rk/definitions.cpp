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

EssentialBoundaryCondition::EssentialBCValueType CustomDirichletCondition::get_value_type() const
{ 
  return EssentialBoundaryCondition::BC_FUNCTION; 
}

scalar CustomDirichletCondition::value(double x, double y, double n_x, double n_y, 
                                       double t_x, double t_y) const 
{
  return x*(100 - x)/2.5 * pow(y/100, y_power) - 1000;
}

/* Custom Initial condition */

void CustomInitialCondition::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  dx = (100 - 2*x)/2.5 * pow(y/100, y_power);
  dy = x*(100 - x)/2.5 * pow(y/100, y_power - 1) * 1./100;
}

scalar CustomInitialCondition::value (double x, double y) const 
{
  return x*(100 - x)/2.5 * pow(y/100, y_power) - 1000;
}

Ord CustomInitialCondition::ord(Ord x, Ord y) const 
{
  return x*(100 - x)/2.5 * pow(y/100, y_power) - 1000;
}

CustomWeakFormRichardsRK::CustomWeakFormRichardsRK() : WeakForm(1)
{
  // Jacobian volumetric part.
  add_matrix_form(new CustomFormMatrixFormVol(0, 0));

  // Residual - volumetric.
  add_vector_form(new CustomFormVectorFormVol1(0));
  //add_vector_form(new CustomFormVectorFormVol2(0));
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormRichardsRK::CustomFormMatrixFormVol::matrix_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * 
              (- (K(u->val[i]) / C(u->val[i])) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
               + (K(u->val[i]) * dCdh(u->val[i])) / (C(u->val[i]) * C(u->val[i])) * (u->dx[i] * u->dx[i] + u->dy[i] * u->dy[i]) * v->val[i]
               + (dKdh(u->val[i]) / C(u->val[i])) * u->dy[i] * v->val[i]);

  }
  return result;
}

scalar CustomWeakFormRichardsRK::CustomFormMatrixFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                                                                Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return matrix_form_rk<double, scalar>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakFormRichardsRK::CustomFormMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                                                           Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);
}

WeakForm::MatrixFormVol* CustomWeakFormRichardsRK::CustomFormMatrixFormVol::clone() 
{
  return new CustomFormMatrixFormVol(*this);
}


template<typename Real, typename Scalar>
Scalar CustomWeakFormRichardsRK::CustomFormVectorFormVol1::vector_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                                                         Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  Func<Scalar>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * 
              (- (K(h_prev_newton->val[i]) / C(h_prev_newton->val[i])) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
               + (K(h_prev_newton->val[i]) * dCdh(h_prev_newton->val[i])) / (C(h_prev_newton->val[i]) * C(h_prev_newton->val[i])) * (h_prev_newton->dx[i] * h_prev_newton->dx[i] + h_prev_newton->dy[i] * h_prev_newton->dy[i]) * v->val[i]
               + (dKdh(h_prev_newton->val[i]) / C(h_prev_newton->val[i])) * h_prev_newton->dy[i] * v->val[i]);

  }
  return result;
}

scalar CustomWeakFormRichardsRK::CustomFormVectorFormVol1::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                                                                ExtData<scalar> *ext) const 
{
  return vector_form_rk<double, scalar>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormRichardsRK::CustomFormVectorFormVol1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);
}

WeakForm::VectorFormVol* CustomWeakFormRichardsRK::CustomFormVectorFormVol1::clone() 
{
  return new CustomFormVectorFormVol1(*this);
}


template<typename Real, typename Scalar>
Scalar CustomWeakFormRichardsRK::CustomFormVectorFormVol2::vector_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                                                         Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
  {
    Scalar C2 = C(h_prev_newton->val[i]) * C(h_prev_newton->val[i]);
    Scalar a1_1 = -(dKdh(h_prev_newton->val[i]) * C(h_prev_newton->val[i]) - K(h_prev_newton->val[i]) * dCdh(h_prev_newton->val[i])) / C2;
    Scalar a1_2 = -(K(h_prev_newton->val[i]) / C(h_prev_newton->val[i]));
    Scalar a2_1 = (dKdh(h_prev_newton->val[i]) * dCdh(h_prev_newton->val[i]) + K(h_prev_newton->val[i]) * ddCdhh(h_prev_newton->val[i])) / C2
                   - K(h_prev_newton->val[i]) * dCdh(h_prev_newton->val[i]) * 2 * C(h_prev_newton->val[i]) * dCdh(h_prev_newton->val[i]) / (C2 * C2);
    Scalar a2_2 = K(h_prev_newton->val[i]) * dCdh(h_prev_newton->val[i]) / C2;
    Scalar a3_1 = (ddKdhh(h_prev_newton->val[i]) * C(h_prev_newton->val[i]) - dKdh(h_prev_newton->val[i]) * dCdh(h_prev_newton->val[i])) / C2;
    Scalar a3_2 = dKdh(h_prev_newton->val[i]) / C(h_prev_newton->val[i]);

    result += wt[i] * (a1_1 * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i]) * v->val[i] 
                       + a1_2 * (v->dx[i] * v->dx[i] + v->dy[i] * v->dy[i]) 
                       + a2_1 * (h_prev_newton->dx[i] * h_prev_newton->dx[i] + h_prev_newton->dy[i] * h_prev_newton->dy[i]) * v->val[i] * v->val[i] 
                       + a2_2 * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i]) * v->val[i] * 2
                       + a3_1 * h_prev_newton->dy[i] * v->val[i] * v->val[i]
                       + a3_2 * v->dy[i] * v->val[i] * v->val[i]
                      );

  }
  return result;
}

scalar CustomWeakFormRichardsRK::CustomFormVectorFormVol2::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                                                                ExtData<scalar> *ext) const 
{
  return vector_form_rk<double, scalar>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormRichardsRK::CustomFormVectorFormVol2::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);
}

WeakForm::VectorFormVol* CustomWeakFormRichardsRK::CustomFormVectorFormVol2::clone() 
{
  return new CustomFormVectorFormVol2(*this);
}
