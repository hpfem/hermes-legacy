#include "definitions.h"

CustomWeakFormLinearElasticity::CustomWeakFormLinearElasticity(double E, double nu, double rho_g) : WeakForm(2)
{
  double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
  double mu = E / (2*(1 + nu));
       
#ifdef USE_MULTICOMPONENT_FORMS
  // Jacobian matrix.
  // There is one multi-component and one single-component form since we want to exploit symmetry of the forms.
  add_multicomponent_matrix_form(new WeakFormsElasticity::DefaultJacobianElasticity_00_11(
        Hermes::vector<std::pair<unsigned int, unsigned int> >(make_pair(0, 0), make_pair(1, 1)), 
        HERMES_ANY, lambda, mu));
  add_matrix_form(new WeakFormsElasticity::DefaultJacobianElasticity_0_1(0, 1, lambda, mu));

  // Residual.
  add_multicomponent_vector_form(new WeakFormsElasticity::DefaultResidualElasticity_00_11(
				  Hermes::vector<unsigned int>(0, 1), lambda, mu));
  add_vector_form(new WeakFormsElasticity::DefaultResidualElasticity_0_1(0, lambda, mu));
  add_vector_form(new WeakFormsElasticity::DefaultResidualElasticity_1_0(1, lambda, mu));

  // Gravity loading.
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol(1, HERMES_ANY, new HermesFunction(-rho_g)));

#else 
  // SINGLE-COMPONENT FORMS. USEFUL FOR MULTIMESH, DO NOT REMOVE.
  // Jacobian.
  add_matrix_form(new WeakFormsElasticity::DefaultJacobianElasticity_0_0(0, 0, lambda, mu));
  add_matrix_form(new WeakFormsElasticity::DefaultJacobianElasticity_0_1(0, 1, lambda, mu));
  add_matrix_form(new WeakFormsElasticity::DefaultJacobianElasticity_1_1(1, 1, lambda, mu));

  // Residual - first equation.
  add_vector_form(new WeakFormsElasticity::DefaultResidualElasticity_0_0(0, HERMES_ANY, lambda, mu));
  add_vector_form(new WeakFormsElasticity::DefaultResidualElasticity_0_1(0, HERMES_ANY, lambda, mu));

  // Residual - second equation.
  add_vector_form(new WeakFormsElasticity::DefaultResidualElasticity_1_0(1, HERMES_ANY, lambda, mu));
  add_vector_form(new WeakFormsElasticity::DefaultResidualElasticity_1_1(1, HERMES_ANY, lambda, mu));
  // Gravity loading in the second vector component.
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol(1, HERMES_ANY, new HermesFunction(-rho_g)));

#endif
}
