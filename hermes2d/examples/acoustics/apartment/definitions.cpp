#include "definitions.h"

CustomWeakFormAcoustics::CustomWeakFormAcoustics(std::string bdy_newton, double rho, double sound_speed, double omega) : WeakForm(1) 
{
  scalar ii =  cplx(0.0, 1.0);

  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, HERMES_ANY, new HermesFunction(1.0/rho), HERMES_SYM));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 0, HERMES_ANY, new HermesFunction(-sqr(omega) / rho / sqr(sound_speed)), HERMES_SYM));
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf(0, 0, bdy_newton, new HermesFunction(-ii * omega / rho / sound_speed)));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, HERMES_ANY, new HermesFunction(1.0/rho)));
  add_vector_form(new WeakFormsH1::DefaultResidualVol(0, HERMES_ANY, new HermesFunction(-sqr(omega) / rho / sqr(sound_speed))));
  add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf(0, bdy_newton, new HermesFunction(-ii * omega / rho / sound_speed)));
}

