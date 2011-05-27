#include "definitions.h"

CustomWeakForm::CustomWeakForm(double K_squared) : WeakForm(1) 
{
  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 0, HERMES_ANY, new HermesFunction(K_squared)));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0));
  add_vector_form(new WeakFormsH1::DefaultResidualVol(0, HERMES_ANY, new HermesFunction(K_squared)));
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol(0, HERMES_ANY, new HermesFunction(-K_squared)));
}


