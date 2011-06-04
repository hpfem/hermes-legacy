#include "definitions.h"

CustomWeakFormHeatRK1::CustomWeakFormHeatRK1(std::string bdy_air, double alpha, double lambda, double heatcap, double rho,
                                             double time_step, double* current_time_ptr, double temp_init, double t_final,
                                             Solution* prev_time_sln) : WeakForm(1)
{
  /* Jacobian */
  // Contribution of the time derivative term.
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 0, HERMES_ANY, new HermesFunction(1.0 / time_step)));
  // Contribution of the diffusion term.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, HERMES_ANY, new HermesFunction(lambda / (rho * heatcap))));
  // Contribution of the Newton boundary condition.
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf(0, 0, bdy_air, new HermesFunction(alpha / (rho * heatcap))));

  /* Residual */
  // Contribution of the time derivative term.
  add_vector_form(new WeakFormsH1::DefaultResidualVol(0, HERMES_ANY, new HermesFunction(1.0 / time_step)));
  // Contribution of the diffusion term.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, HERMES_ANY, new HermesFunction(lambda / (rho * heatcap))));
  CustomVectorFormVol* vec_form_vol = new CustomVectorFormVol(0, time_step);
  vec_form_vol->ext.push_back(prev_time_sln);
  add_vector_form(vec_form_vol);
  // Contribution of the Newton boundary condition.
  add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf(0, bdy_air, new HermesFunction(alpha / (rho * heatcap))));
  // Contribution of the Newton boundary condition.
  add_vector_form_surf(new CustomVectorFormSurf(0, bdy_air, alpha, rho, heatcap,
                       time_step, current_time_ptr, temp_init, t_final));
}

scalar CustomWeakFormHeatRK1::CustomVectorFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  Func<double>* temp_prev_time = ext->fn[0];
  return -int_u_v<double, scalar>(n, wt, temp_prev_time, v) / time_step;
}

Ord CustomWeakFormHeatRK1::CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  Func<Ord>* temp_prev_time = ext->fn[0];
  return -int_u_v<Ord, Ord>(n, wt, temp_prev_time, v) / time_step;
}

scalar CustomWeakFormHeatRK1::CustomVectorFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return -alpha / (rho * heatcap) * temp_ext(*current_time_ptr + time_step) * int_v<double>(n, wt, v);
}

Ord CustomWeakFormHeatRK1::CustomVectorFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return -alpha / (rho * heatcap) * temp_ext(*current_time_ptr + time_step) * int_v<Ord>(n, wt, v);
}

// Time-dependent exterior temperature.
template<typename Real>
Real CustomWeakFormHeatRK1::CustomVectorFormSurf::temp_ext(Real t) const 
{
  return temp_init + 10. * sin(2*M_PI*t/t_final);
}
