#include "definitions.h"
#include "weakform_library/weakforms_h1.h"

//////////////////////////////////////////////////// CONTINUOUS APPROXIMATION ////////////////////////////////////////////////////

const double ConstFlowField::b1 = 0.;
const double ConstFlowField::b2 = 1.;

CustomWeakFormContinuousGalerkin::CustomWeakFormContinuousGalerkin(GalerkinMethod method, double epsilon)
  : WeakForm(1), fn_epsilon(new HermesFunction(epsilon)), 
    fn_b1(new HermesFunction(ConstFlowField::b1)), fn_b2(new HermesFunction(ConstFlowField::b2))
{ 
  add_matrix_form(new WeakFormsH1::DefaultJacobianAdvection(0, 0, HERMES_ANY, fn_b1, fn_b2));
  add_vector_form(new WeakFormsH1::DefaultResidualAdvection(0, HERMES_ANY, fn_b1, fn_b2));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, HERMES_ANY, fn_epsilon, HERMES_SYM));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, HERMES_ANY, fn_epsilon));
  
  if (method != CG) 
  {
    add_matrix_form(new StabilizationJacobian(method, epsilon));
    add_vector_form(new StabilizationResidual(method, epsilon));
  }
}

scalar CustomWeakFormContinuousGalerkin::StabilizationJacobian::value(int n, double *wt, 
                                                                      Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                                                                      Geom<double> *e, ExtData<scalar> *ext) const
{
  scalar result = 0;
  double h_e = e->diam;
  double Pe = element_Peclet_number(h_e);
  
  switch(method)
  {
    case CG_STAB_SUPG:
      // Streamline upwind Petrov-Galerkin stabilization, stabilization parameter according to Codina 
      // (derived from nodal exactness).
      for (int i=0; i < n; i++) 
      {
        double alpha = 1 + 2./(exp(2*Pe) - 1) - 1./Pe;  // coth(Pe)-1/Pe
        double tau = alpha*h_e / (2*b_norm);
        
        result += wt[i] * tau * ConstFlowField::dot_grad<double>(v,e,i) * 
                               (ConstFlowField::dot_grad<double>(u,e,i) - epsilon * u->laplace[i]);
      }
      
      break; 
    case CG_STAB_GLS_1:
    case CG_STAB_GLS_2:
      // Galerkin least-squares stabilization, stabilization parameter according to Codina.
      for (int i=0; i < n; i++) 
      {
        double C1 = (method == CG_STAB_GLS_1) ? 1./3. : 1./9.;
        double C2 = (method == CG_STAB_GLS_1) ? 1.    : 1./2.;
        double tau = std::min(C1*Pe, C2)*h_e / (2*b_norm);
        
        result += wt[i] * tau * (ConstFlowField::dot_grad<double>(v,e,i) - epsilon * v->laplace[i])
                              * (ConstFlowField::dot_grad<double>(u,e,i) - epsilon * u->laplace[i]);
      }
      break;
    case CG_STAB_SGS:
      // Subgrid scale stabilization, stabilization parameter from discrete max. principle, according to Codina.
      for (int i=0; i < n; i++) 
      {
        double tau = 1. / (4*epsilon/sqr(h_e) + 2*b_norm/h_e);
        result += wt[i] * tau * (ConstFlowField::dot_grad<double>(v,e,i) + epsilon * v->laplace[i])
                              * (ConstFlowField::dot_grad<double>(u,e,i) - epsilon * u->laplace[i]);
      }
    case CG_STAB_SGS_ALT:
      // Subgrid scale stabilization, stabilization parameter according to Shakib.
      for (int i=0; i < n; i++) 
      {
        double tau = 1. / sqrt(9*sqr(4*epsilon/sqr(h_e)) + sqr(2*b_norm/h_e));
        result += wt[i] * tau * (ConstFlowField::dot_grad<double>(v,e,i) + epsilon * v->laplace[i])
                              * (ConstFlowField::dot_grad<double>(u,e,i) - epsilon * u->laplace[i]);
      }      
  }
  
  return result;
}

Ord CustomWeakFormContinuousGalerkin::StabilizationJacobian::ord(int n, double* wt, 
                                                                 Func< Ord >* u_ext[], Func< Ord >* u, Func< Ord >* v, 
                                                                 Geom< Ord >* e, ExtData< Ord >* ext) const
{
  if (method == CG_STAB_SUPG)
    return ConstFlowField::dot_grad<Ord>(v,e,0) * (ConstFlowField::dot_grad<Ord>(u,e,0) - epsilon * u->laplace[0]);   
  else
    return (ConstFlowField::dot_grad<Ord>(v,e,0) + epsilon * v->laplace[0]) *
           (ConstFlowField::dot_grad<Ord>(u,e,0) + epsilon * u->laplace[0]);
}

scalar CustomWeakFormContinuousGalerkin::StabilizationResidual::value(int n, double *wt, 
                                                                      Func<scalar> *u_ext[], Func<double> *v,
                                                                      Geom<double> *e, ExtData<scalar> *ext) const
{
  scalar result = 0;
  double h_e = e->diam;
  double Pe = element_Peclet_number(h_e);
  
  switch(method)
  {
    case CG_STAB_SUPG:
      // Streamline upwind Petrov-Galerkin stabilization, stabilization parameter according to Codina 
      // (derived from nodal exactness).
      for (int i=0; i < n; i++) 
      {
        double alpha = 1 + 2./(exp(2*Pe) - 1) - 1./Pe;  // coth(Pe)-1/Pe
        double tau = alpha*h_e / (2*b_norm);
        
        result += wt[i] * tau * ConstFlowField::dot_grad<double>(v,e,i) * 
                               (ConstFlowField::dot_grad<double>(u_ext[0],e,i) - epsilon * u_ext[0]->laplace[i]);
      }
      
      break; 
    case CG_STAB_GLS_1:
    case CG_STAB_GLS_2:
      // Galerkin least-squares stabilization, stabilization parameter according to Codina.
      for (int i=0; i < n; i++) 
      {
        double C1 = (method == CG_STAB_GLS_1) ? 1./3. : 1./9.;
        double C2 = (method == CG_STAB_GLS_1) ? 1.    : 1./2.;
        double tau = std::min(C1*Pe, C2)*h_e / (2*b_norm);
        
        result += wt[i] * tau * (ConstFlowField::dot_grad<double>(v,e,i) - epsilon * v->laplace[i])
                              * (ConstFlowField::dot_grad<double>(u_ext[0],e,i) - epsilon * u_ext[0]->laplace[i]);
      }
      break;
    case CG_STAB_SGS:
      // Subgrid scale stabilization, stabilization parameter from discrete max. principle, according to Codina.
      for (int i=0; i < n; i++) 
      {
        double tau = 1. / (4*epsilon/sqr(h_e) + 2*b_norm/h_e);
        result += wt[i] * tau * (ConstFlowField::dot_grad<double>(v,e,i) + epsilon * v->laplace[i])
                              * (ConstFlowField::dot_grad<double>(u_ext[0],e,i) - epsilon * u_ext[0]->laplace[i]);
      }
    case CG_STAB_SGS_ALT:
      // Subgrid scale stabilization, stabilization parameter according to Shakib.
      for (int i=0; i < n; i++) 
      {
        double tau = 1. / sqrt(9*sqr(4*epsilon/sqr(h_e)) + sqr(2*b_norm/h_e));
        result += wt[i] * tau * (ConstFlowField::dot_grad<double>(v,e,i) + epsilon * v->laplace[i])
                              * (ConstFlowField::dot_grad<double>(u_ext[0],e,i) - epsilon * u_ext[0]->laplace[i]);
      }      
  }
  
  return result;
}

Ord CustomWeakFormContinuousGalerkin::StabilizationResidual::ord(int n, double* wt, 
                                                                 Func< Ord >* u_ext[], Func< Ord >* v, 
                                                                 Geom< Ord >* e, ExtData< Ord >* ext) const
{
  if (method == CG_STAB_SUPG)
    return ConstFlowField::dot_grad<Ord>(v,e,0) * (ConstFlowField::dot_grad<Ord>(u_ext[0],e,0) - epsilon * u_ext[0]->laplace[0]);   
  else
    return (ConstFlowField::dot_grad<Ord>(v,e,0) + epsilon * v->laplace[0]) *
           (ConstFlowField::dot_grad<Ord>(u_ext[0],e,0) + epsilon * u_ext[0]->laplace[0]);
}

//////////////////////////////////////////////////// DISCONTINUOUS APPROXIMATION ////////////////////////////////////////////////////

CustomWeakFormDiscontinuousGalerkin::CustomWeakFormDiscontinuousGalerkin(const EssentialBCs& boundary_values,
                                                                        double epsilon, int theta, int C_W)
  : WeakForm(1), theta(theta), C_W(C_W), fn_epsilon(new HermesFunction(epsilon))
{
  add_matrix_form(new Advection::VolumetricJacobian);
  add_vector_form(new Advection::VolumetricResidual);
  add_matrix_form_surf(new Advection::BoundaryJacobian);
  add_vector_form_surf(new Advection::BoundaryResidual(boundary_values));
  add_matrix_form_surf(new Advection::InterfaceJacobian);
  add_vector_form_surf(new Advection::InterfaceResidual);
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0,0,HERMES_ANY,fn_epsilon,HERMES_SYM));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0,HERMES_ANY,fn_epsilon));
  add_matrix_form_surf(new Diffusion::BoundaryJacobian(epsilon, theta, C_W));
  add_vector_form_surf(new Diffusion::BoundaryResidual(boundary_values, epsilon, theta, C_W));
  add_matrix_form_surf(new Diffusion::InterfaceJacobian(epsilon, theta, C_W));
  add_vector_form_surf(new Diffusion::InterfaceResidual(epsilon, theta, C_W));
}


// Weak forms.

#define JUMP(w)       ( w->get_val_central(i) - w->get_val_neighbor(i) )
#define AVG_GRAD(w)   ( 0.5 * ( (w->get_dx_central(i) + w->get_dx_neighbor(i))*e->nx[i] + \
                                (w->get_dy_central(i) + w->get_dy_neighbor(i))*e->ny[i] ) )

//--- ADVECTION

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::Advection::VolumetricJacobian::matrix_form(int n, double *wt, 
                                                                            Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                                                                            Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += -wt[i] * u->val[i] * ConstFlowField::dot_grad<Real>(v,e,i);
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::Advection::VolumetricResidual::vector_form(int n, double *wt, 
                                                                            Func<Scalar> *u_ext[], Func<Real> *v, 
                                                                            Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += -wt[i] * u_ext[0]->val[i] * ConstFlowField::dot_grad<Real>(v,e,i);
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::Advection::BoundaryJacobian::matrix_form(int n, double* wt,
                                                                                     Func< Scalar >* u_ext[], Func< Real >* u, Func< Real >* v, 
                                                                                     Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++)
    result += wt[i] * upwind_flux(u->val[i], 0, ConstFlowField::dot_n<Real>(e,i)) * v->val[i];
  
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::Advection::BoundaryResidual::vector_form(int n, double* wt, 
                                                                                     Func< Scalar >* u_ext[], Func< Real >* v, 
                                                                                     Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0;
  
  std::string marker = wf->get_boundary_markers_conversion()->get_user_marker(e->edge_marker);
  
  for (int i = 0; i < n; i++) 
  {
    Real beta_dot_n = ConstFlowField::dot_n<Real>(e,i);
    WeaklyImposableBC* bc = static_cast<WeaklyImposableBC*>(boundary_values.get_boundary_condition(marker));
    result += ( upwind_flux(u_ext[0]->val[i], 0, beta_dot_n) 
                + upwind_flux(0, bc->value(e->x[i], e->y[i]), beta_dot_n) ) * wt[i] * v->val[i];
  }
  
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::Advection::InterfaceJacobian::matrix_form(int n, double* wt, 
                                                                                      Func< Scalar >* u_ext[], Func< Real >* u, Func< Real >* v, 
                                                                                      Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++)
    result += wt[i] * upwind_flux( u->get_val_central(i), u->get_val_neighbor(i), ConstFlowField::dot_n<Real>(e,i) ) * JUMP(v);
  
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::Advection::InterfaceResidual::vector_form(int n, double* wt, 
                                                                                      Func< Scalar >* u_ext[], Func< Real >* v,
                                                                                      Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++)
    result += upwind_flux( u_ext[0]->get_val_central(i), u_ext[0]->get_val_neighbor(i), ConstFlowField::dot_n<Real>(e,i) )
              * JUMP(v) * wt[i];
  
  return result;
}

//--- DIFFUSION AND PENALTIES

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::Diffusion::BoundaryJacobian::matrix_form(int n, double* wt, 
                                                                                     Func< Scalar >* u_ext[], Func< Real >* u, Func< Real >* v, 
                                                                                     Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0;
  //Real sigma = C_W * EPSILON / e->diam;
  Real edge_len = 0.;
  for (int i = 0; i < n; i++)
    edge_len += wt[i];
  
  Real sigma = C_W * epsilon / (0.5*edge_len);
  
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * epsilon * ( -dot2<Real>(u->dx[i], u->dy[i], e->nx[i], e->ny[i]) * v->val[i]
                                  +dot2<Real>(v->dx[i], v->dy[i], e->nx[i], e->ny[i]) * theta * u->val[i] );   
    result += wt[i] * sigma * u->val[i] * v->val[i];                    
  }
  
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::Diffusion::BoundaryResidual::vector_form(int n, double* wt, 
                                                                                     Func< Scalar >* u_ext[], Func< Real >* v, 
                                                                                     Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0;
  //Real sigma = C_W * EPSILON / e->diam;
  Real edge_len = 0.;
  for (int i = 0; i < n; i++)
    edge_len += wt[i];
  
  Real sigma = C_W * epsilon / (0.5*edge_len);
  
  std::string marker = wf->get_boundary_markers_conversion()->get_user_marker(e->edge_marker);
  
  for (int i = 0; i < n; i++)
  {
    WeaklyImposableBC* bc = static_cast<WeaklyImposableBC*>(boundary_values.get_boundary_condition(marker));
    Scalar bnd_diff = u_ext[0]->val[i] - bc->value(e->x[i], e->y[i]);
    result += wt[i] * epsilon * ( -dot2<Real>(u_ext[0]->dx[i], u_ext[0]->dy[i], e->nx[i], e->ny[i]) * v->val[i]
                                  +dot2<Real>(v->dx[i], v->dy[i], e->nx[i], e->ny[i]) * theta * bnd_diff );   
    result += wt[i] * sigma * bnd_diff * v->val[i];                        
  }
  
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::Diffusion::InterfaceJacobian::matrix_form(int n, double* wt, 
                                                                                      Func< Scalar >* u_ext[], Func< Real >* u, Func< Real >* v, 
                                                                                      Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0;
  //Real sigma = 2 * C_W / (e->diam + e->get_neighbor_diam());
  Real edge_len = 0.;
  for (int i = 0; i < n; i++)
    edge_len += wt[i];
  
  Real sigma = C_W * epsilon / (0.5*edge_len);
  
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * epsilon * (-AVG_GRAD(u) * JUMP(v) + theta * AVG_GRAD(v) * JUMP(u)); // diffusion
    result += wt[i] * sigma * JUMP(u) * JUMP(v);                                          // interior discontinuity penalization
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::Diffusion::InterfaceResidual::vector_form(int n, double* wt,
                                                                                      Func< Scalar >* u_ext[], Func< Real >* v, 
                                                                                      Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0;
  //Real sigma = 2 * C_W / (e->diam + e->get_neighbor_diam());
  Real edge_len = 0.;
  for (int i = 0; i < n; i++)
    edge_len += wt[i];
  
  Real sigma = C_W * epsilon / (0.5*edge_len);
  
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * epsilon * (-AVG_GRAD(u_ext[0]) * JUMP(v) + theta * AVG_GRAD(v) * JUMP(u_ext[0])); 
    result += wt[i] * sigma * JUMP(u_ext[0]) * JUMP(v);                                                 
  }
  return result;
}