#include "definitions.h"

const double FuelProperties::Tref = 0.0;              // Reference temperature for defining temperature dependence of parameters.
const double FuelProperties::neutron_velocity = 5e3;  // Averaged neutron speed.
const double FuelProperties::D = 1.268;               // Neutron diffusion coefficient.
const double FuelProperties::nu = 2.41;               // Number of neutrons emitted per fission event.
const double FuelProperties::Sigma_f = 0.00191244;    // Fission cross section.
const double FuelProperties::kappa = 1.0e-6;          // Energy per fission.
const double FuelProperties::rho = 1.0;               // Fuel density.
const double FuelProperties::cp = 1.0;                // Fuel heat capacity.
const double FuelProperties::k0 = 3.0e-3;             // Constant coefficient of heat conduction dependence on temperature.
const double FuelProperties::k1 = 2.0e-4;             // Linear coefficient of heat conduction dependence on temperature.
const double FuelProperties::Sigma_r_ref = 0.0349778; // Absorption cross section for reference temperature Tref.
const double FuelProperties::doppler_coeff = 1.0e-5;  // Coefficient for modelling the dependence of absorption cross section
                                                      // on temperature (Doppler effect).

// Heat source.
template<typename Real>
Real TemperatureField::SourceTerm::val(Real x, Real y, double t) const
{
  double dTdt = TemperatureField::transient_profile_derivative(t);
  double Tt = TemperatureField::transient_profile(t);
  double PHIt = NeutronField::transient_profile(t);
  
  double PI_sqr = sqr(M_PI);
  Real sx = sin((M_PI*x)/LX);
  Real sy = sin((M_PI*y)/LY);
  
  return CT*dTdt*sx*sy 
          - (NeutronField::CF*fuel.kappa*PHIt*x*fuel.Sigma_f*y*sx*sy)/(LX*LY*fuel.rho*fuel.cp) 
          - ( -(CT*PI_sqr*Tt*sx*sy)/sqr(LX) - (CT*PI_sqr*Tt*sx*sy)/sqr(LY) ) 
          * (fuel.k0 + fuel.k1*(-fuel.Tref + CT*Tt*sx*sy)) / (fuel.rho*fuel.cp);
}

// Neutron source.
template<typename Real>
Real NeutronField::SourceTerm::val(Real x, Real y, double t) const
{
  double PHIt = NeutronField::transient_profile(t);
  double dPHIdt = NeutronField::transient_profile_derivative(t);
  double Tt = TemperatureField::transient_profile(t);
  
  Real PI_sqr = sqr(M_PI);
  Real sx = sin((M_PI*x)/LX);
  Real sy = sin((M_PI*y)/LY);
  
  return (CF*dPHIdt*x*y*sx*sy)/(LX*LY) - fuel.neutron_velocity * (
            fuel.D * ( 
              (2*CF*PHIt*M_PI*x*cos((M_PI*y)/LY)*sx)/(LX*sqr(LY)) + (2*CF*PHIt*M_PI*y*cos((M_PI*x)/LX)*sy)/(sqr(LX)*LY) 
              -(CF*PHIt*PI_sqr*x*y*sx*sy)/(LX*pow(LY,3)) - (CF*PHIt*PI_sqr*x*y*sx*sy)/(pow(LX,3)*LY)
            ) + (
              CF*PHIt*x*y*sx*sy*(
                fuel.nu*fuel.Sigma_f-fuel.Sigma_r_ref*( 1 + fuel.doppler_coeff*(-sqrt(fuel.Tref) + sqrt(TemperatureField::CT*Tt*sx*sy)) )
              )
            )/(LX*LY)
          );
}

template<typename Real, typename Scalar>
Scalar TemperatureField::Residual::vector_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0.0;
  
  Func<Scalar>* T_prev_newton = u_ext[0];
  Func<Scalar>* phi_prev_newton = u_ext[1];
  
  double t = get_current_stage_time();
  
  for (int i = 0; i < n; i++)
    result += ( minus_lambda.value(T_prev_newton->val[i]) * (T_prev_newton->dx[i] * v->dx[i] + T_prev_newton->dy[i] * v->dy[i])
                + fission_source_coeff * phi_prev_newton->val[i] * v->val[i]
                + qT.value(e->x[i], e->y[i], t) * v->val[i] ) * wt[i];
  
  return result;
}

template<typename Real, typename Scalar>
Scalar NeutronField::NeutronFluxDerivative::matrix_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* u, Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0.0;
  Func<Scalar>* T_prev_newton = u_ext[0];
  
  for (int i = 0; i < n; i++)
    result += ( - diffusion_coeff * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                - Sigma_r.value(T_prev_newton->val[i]) * u->val[i] * v->val[i]
                + fission_yield_coeff * u->val[i] * v->val[i] ) * wt[i];
  
  return result;
}

template<typename Real, typename Scalar>
Scalar NeutronField::TemperatureDerivative::matrix_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* u, Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0.0;
  
  Func<Scalar>* T_prev_newton = u_ext[0];
  Func<Scalar>* phi_prev_newton = u_ext[1];
  
  for (int i = 0; i < n; i++)
    result += wt[i] * (-Sigma_r.derivative(T_prev_newton->val[i]) * u->val[i] * phi_prev_newton->val[i] * v->val[i]);
  
  return result;  
}

template<typename Real, typename Scalar>
Scalar NeutronField::Residual::vector_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0.0;
  
  Func<Scalar>* T_prev_newton = u_ext[0];
  Func<Scalar>* phi_prev_newton = u_ext[1];
 
  double t = get_current_stage_time();
  
  for (int i = 0; i < n; i++)
    result += ( - diffusion_coeff * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i])
                - Sigma_r.value(T_prev_newton->val[i]) * phi_prev_newton->val[i] * v->val[i] 
                + fission_yield_coeff * phi_prev_newton->val[i] * v->val[i]
                + qF.value(e->x[i], e->y[i], t) * v->val[i] ) * wt[i];
  
  return result;
}

CustomWeakForm::CustomWeakForm(double lx, double ly) : WeakForm(2)
{
  // Jacobian 
  
  // dF_T/dT
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, HERMES_ANY, 
                                                            new FuelProperties::NegativeConductivityTemperatureProfile, 
                                                            HERMES_SYM));

  // dF_T/dPhi
  double coeff = FuelProperties::kappa / (FuelProperties::rho * FuelProperties::cp) * FuelProperties::Sigma_f;
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 1, HERMES_ANY, new HermesFunction(coeff)));
  
  // dF_Phi/dT
  add_matrix_form(new NeutronField::TemperatureDerivative);
  
  // dF_Phi/dPhi
  add_matrix_form(new NeutronField::NeutronFluxDerivative);
  
  // Residual
  
  // F_T
  add_vector_form(new TemperatureField::Residual(lx, ly));
  
  // F_Phi
  add_vector_form(new NeutronField::Residual(lx, ly));
}

Views::Views(unsigned int width, unsigned int height) 
{
  // Initialize solution views (their titles will be updated in each time step).
  sview_T = new ScalarView("", new WinGeom(0, 0, width, height));
  sview_T->fix_scale_width(50);
  sview_phi = new ScalarView("", new WinGeom(0, height+50, width, height));
  sview_phi->fix_scale_width(50);
  sview_T_exact = new ScalarView("", new WinGeom(width, 0, width, height));
  sview_T_exact->fix_scale_width(50);
  sview_phi_exact = new ScalarView("", new WinGeom(width, height+50, width, height));
  sview_phi_exact->fix_scale_width(50);
}

Views::~Views()
{
  delete sview_T;
  delete sview_phi;
  delete sview_T_exact;
  delete sview_phi_exact;
}

void Views::show_solutions(double current_time, Hermes::vector<Solution*> solutions)
{
  // Show the new time level solution.
  sprintf(title, "Approx. solution for T, t = %g s", current_time);
  sview_T->set_title(title); 
  sview_T->show(solutions[0]);
  
  sprintf(title, "Approx. solution for phi, t = %g s", current_time);
  sview_phi->set_title(title);
  sview_phi->show(solutions[1]);
}

void Views::show_exact(double current_time, Hermes::vector< Solution* > exact)
{
  // Show exact solution.
  sprintf(title, "Exact solution for T, t = %g s", current_time);
  sview_T_exact->set_title(title);
  sview_T_exact->show(exact[0]);
  
  sprintf(title, "Exact solution for phi, t = %g s", current_time);
  sview_phi_exact->set_title(title);
  sview_phi_exact->show(exact[1]);
}


