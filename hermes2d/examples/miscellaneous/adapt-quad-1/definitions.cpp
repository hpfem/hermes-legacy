#include "definitions.h"

/* Weak forms */

CustomWeakFormPoisson::CustomWeakFormPoisson(std::string mat_al, HermesFunction* lambda_al,
                                             std::string mat_cu, HermesFunction* lambda_cu,
                                             HermesFunction* src_term, bool adapt_eval, 
                                             int adapt_order_increase, double adapt_rel_error_tol) 
  : WeakForm(1)
{
  // Create weak forms.
  // Jacobian forms.
  MatrixFormVol* jacobian_diffusion_al = new WeakFormsH1::DefaultJacobianDiffusion(0, 0, mat_al, lambda_al);
  MatrixFormVol* jacobian_diffusion_cu = new WeakFormsH1::DefaultJacobianDiffusion(0, 0, mat_cu, lambda_cu);
  // Residual forms.
  VectorFormVol* residual_diffusion_al = new WeakFormsH1::DefaultResidualDiffusion(0, mat_al, lambda_al);
  VectorFormVol* residual_diffusion_cu = new WeakFormsH1::DefaultResidualDiffusion(0, mat_cu, lambda_cu);
  VectorFormVol* src_term_form = new WeakFormsH1::DefaultVectorFormVol(0, HERMES_ANY, src_term);

  // Set adaptive quadrature parameters.
  if(adapt_eval) 
  {
    // Order increase.
    jacobian_diffusion_al->adapt_order_increase = adapt_order_increase;
    jacobian_diffusion_cu->adapt_order_increase = adapt_order_increase;
    residual_diffusion_al->adapt_order_increase = adapt_order_increase;
    residual_diffusion_cu->adapt_order_increase = adapt_order_increase;
    src_term_form->adapt_order_increase = adapt_order_increase;
    // Relative error tolerance.
    jacobian_diffusion_al->adapt_rel_error_tol = adapt_rel_error_tol;
    jacobian_diffusion_cu->adapt_rel_error_tol = adapt_rel_error_tol;
    residual_diffusion_al->adapt_rel_error_tol = adapt_rel_error_tol;
    residual_diffusion_cu->adapt_rel_error_tol = adapt_rel_error_tol;
    src_term_form->adapt_rel_error_tol = adapt_rel_error_tol;
  }

  // Register weak forms.
  add_matrix_form(jacobian_diffusion_al);
  add_matrix_form(jacobian_diffusion_cu);
  add_vector_form(residual_diffusion_al);
  add_vector_form(residual_diffusion_cu);
  add_vector_form(src_term_form);
};
