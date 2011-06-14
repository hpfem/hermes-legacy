#include "definitions.h"

double lambda(double x, double y) 
{
  return 1.0;
}

template<typename Real>
Real T_fire_x(Real x) 
{
  return -1./32 * x*x*x + 3./16 * x*x;
}

template<typename Real>
Real T_fire_t(Real t) 
{
  if (0 <= t  &&  t <= 100) return 0;
  if (100 <= t  &&  t <= 600) return 980. / 500 * (t - 100.);
  if (600 <= t  &&  t <= 1800) return 980;
  if (1800 <= t  &&  t <= 3000) return 980 - 980. / 1200 * (t - 1800.);
  return 0.;
}


CustomWeakFormHeatRK::CustomWeakFormHeatRK(std::string bdy_fire, std::string bdy_air,
                                           double alpha_fire, double alpha_air, double rho, double heatcap,
                                           double temp_ext_air, double temp_init, double* current_time_ptr) : WeakForm(1)
{
  // Jacobian - volumetric part.
  add_matrix_form(new CustomJacobianVol(0, 0, rho, heatcap));

  // Jacobian - surface part.
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf(0, 0, bdy_fire, new HermesFunction(-alpha_fire/(rho*heatcap))));
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf(0, 0, bdy_air, new HermesFunction(-alpha_air/(rho*heatcap))));

  // Residual - volumetric part.
  add_vector_form(new CustomFormResidualVol(0, rho, heatcap));

  // Surface residual - bottom boundary.
  CustomFormResidualSurfFire* vec_form_surf_1
    = new CustomFormResidualSurfFire(0, bdy_fire, alpha_fire, rho, heatcap, current_time_ptr);
  add_vector_form_surf(vec_form_surf_1);

  // Surface residual - top boundary.
  add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf(0, HERMES_ANY, new HermesFunction(-alpha_air / (rho*heatcap))));
  add_vector_form_surf(new WeakFormsH1::DefaultVectorFormSurf(0, HERMES_ANY, new HermesFunction(alpha_air* temp_ext_air / (rho*heatcap))));
}

scalar CustomWeakFormHeatRK::CustomJacobianVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, 
                                                      Geom<double> *e, ExtData<scalar> *ext) const 
{
  scalar result = 0;
  for (int i = 0; i < n; i++) 
  {
    result += wt[i] * lambda(e->x[i], e->y[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  }

  return -result / heatcap / rho;
}

Ord CustomWeakFormHeatRK::CustomJacobianVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                 Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  // Return the polynomial order of the gradient increased by five to account for lambda(x, y).
  return (u->dx[0] * v->dx[0] + u->dy[0] * v->dy[0]) * Ord(5);
}

WeakForm::MatrixFormVol* CustomWeakFormHeatRK::CustomJacobianVol::clone() 
{
  return new CustomJacobianVol(*this);
}


scalar CustomWeakFormHeatRK::CustomFormResidualVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                                                          Geom<double> *e, ExtData<scalar> *ext) const 
{
  Func<scalar>* u_prev_newton = u_ext[0];
  scalar result = 0;
  for (int i = 0; i < n; i++) 
  {
    result += wt[i] * lambda(e->x[i], e->y[i])
                    * (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i]);
  }

  return -result / heatcap / rho;
}

Ord CustomWeakFormHeatRK::CustomFormResidualVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                                     Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  Func<Ord>* u_prev_newton = u_ext[0];
  // Return the polynomial order of the gradient increased by five to account for lambda(x, y).
  return (u_prev_newton->dx[0] * v->dx[0] + u_prev_newton->dy[0] * v->dy[0]) * Ord(5);
}

WeakForm::VectorFormVol* CustomWeakFormHeatRK::CustomFormResidualVol::clone() 
{
  return new CustomFormResidualVol(*this);
}


template<typename Real, typename Scalar>
Scalar CustomWeakFormHeatRK::CustomFormResidualSurfFire::vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                                                          Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Func<Scalar>* sln_prev = u_ext[0];

  Scalar result = 0;
  for (int i = 0; i < n; i++) 
  {
    result += wt[i] * (T_fire(e->x[i], *current_time_ptr) - sln_prev->val[i]) * v->val[i];
  }

  return result / heatcap / rho * alpha_fire;
}

scalar CustomWeakFormHeatRK::CustomFormResidualSurfFire::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
Geom<double> *e, ExtData<scalar> *ext) const 
{
  return vector_form_surf<double, scalar>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormHeatRK::CustomFormResidualSurfFire::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                                          Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  // Return the polynomial order of the test function 'v' plus three for T_fire_x(x)
  return v->val[0] * Ord(3);
}

WeakForm::VectorFormSurf* CustomWeakFormHeatRK::CustomFormResidualSurfFire::clone() 
{
  return new CustomFormResidualSurfFire(*this);
}

template<typename Real>
Real CustomWeakFormHeatRK::CustomFormResidualSurfFire::T_fire(Real x, Real t) const 
{
  return T_fire_x(x) * T_fire_t(t) + 20;
}

