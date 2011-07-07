#include "definitions.h"
#include <cstdlib>

// Weighted integral over the specified boundary edge.
double BoundaryIntegral::value(MeshFunction* sln) const
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  RefMap rm;
  
  Mesh* mesh = sln->get_mesh();
  
  double integral = 0.0;
  Element* el;
  for_all_active_elements(el, mesh)
  {
    for (int e = 0; e < el->nvert; e++)
    {
      if (el->en[e]->bnd && 
          el->en[e]->marker == mesh->get_boundary_markers_conversion().get_internal_marker(marker))
      {
        sln->set_active_element(el);
        rm.set_active_element(el);
        
        // Set quadrature order.
        int eo = quad->get_edge_points(e, g_safe_max_order);
        sln->set_quad_order(eo, H2D_FN_VAL);
        
        // Obtain quadrature points, corresponding solution values and tangent vectors (for computing outward normal).
        double* x = rm.get_phys_x(eo);
        double* y = rm.get_phys_y(eo);
        scalar* z = sln->get_fn_values();
        double3* t = rm.get_tangent(e,eo);
        
        // Add contribution from current edge.
        int np = quad->get_num_points(eo);
        double3* pt = quad->get_points(eo);
        
        for (int i = 0; i < np; i++)
        {
          double beta_dot_n = beta.dot( x[i], y[i], t[i][1], -t[i][0] );
          // Weights sum up to two on every edge, therefore the division by two must be present.
          integral += 0.5 * pt[i][2] * t[i][2] * weight_fn(x[i],y[i]) * z[i] * beta_dot_n; 
        } 
      }
    }
  }
  
  return integral;
}

SemiAnalyticSolution::SemiAnalyticSolution(std::string file) : NA(false)
{
  std::string exact_solution_not_available;
  exact_solution_not_available.append("Exact solution could not be read from the supplied file:");
  exact_solution_not_available.append(file);
  
  std::ifstream ifs(file.c_str());
  
  if (ifs.fail())
  {
    warning(exact_solution_not_available.c_str());
    NA = true;
    return;
  }
  
  std::string ns, xs, ys, us, ws;
  
  std::getline(ifs, ns);
  n = strtoul(ns.c_str(), NULL, 0);
  x.reserve(n);
  y.reserve(n);
  u.reserve(n);
  w.reserve(n);
  
  while(!std::getline(ifs, xs, '\t').eof())
  {
    x.push_back(strtold(xs.c_str(), NULL));
    std::getline(ifs, ys, '\t');
    y.push_back(strtold(ys.c_str(), NULL));
    std::getline(ifs, us, '\t');
    u.push_back(strtold(us.c_str(), NULL));
    std::getline(ifs, ws);
    w.push_back(strtold(ws.c_str(), NULL));
  }
  
  if (x.size() != n || y.size() != n || u.size() != n || w.size() != n)
  {
    warning(exact_solution_not_available.c_str());
    NA = true;
  }
  
  ifs.close();
}

double SemiAnalyticSolution::get_l2_norm()
{
  if (NA) return -1;
  
  long double res = 0.0;
  
  for (unsigned int i = 0; i < n; i++)
    res += w[i] * u[i] * u[i];
  
  return sqrt(res);
}

double SemiAnalyticSolution::get_l2_rel_err(Solution *approx_sln)
{
//   if (NA) return -1;
  
  long double res = 0.0, nrm = 0.0;
  
  for (unsigned int i = 0; i < n; i++)
  {
    nrm += w[i] * u[i] * u[i];
    res += w[i] * pow(u[i] - approx_sln->get_pt_value(x[i],y[i]),2);
  }
  
  return sqrt(res/nrm);
}

//////////////////////////////////////////////////// CONTINUOUS APPROXIMATION ////////////////////////////////////////////////////

CustomWeakFormContinuousGalerkin::CustomWeakFormContinuousGalerkin(bool supg_stabilization) : WeakForm(1)
{
  DiscretizationMethod disc_method = supg_stabilization ? SUPG : CG;
  
  add_matrix_form(new VolumetricJacobian(disc_method));
  add_vector_form(new VolumetricResidual(disc_method));
  
  add_matrix_form_surf(new InflowBoundaryJacobian("nonzero_constant_inflow", new HermesFunction(1.0)));
  add_vector_form_surf(new InflowBoundaryResidual("nonzero_constant_inflow", new HermesFunction(1.0)));
  add_matrix_form_surf(new InflowBoundaryJacobian("zero_inflow", new HermesFunction(0.0)));
  add_vector_form_surf(new InflowBoundaryResidual("zero_inflow", new HermesFunction(0.0)));
  add_matrix_form_surf(new InflowBoundaryJacobian("varying_inflow", new InflowVariation));
  add_vector_form_surf(new InflowBoundaryResidual("varying_inflow", new InflowVariation));
  
  if (supg_stabilization)
  {
    add_matrix_form(new StabilizationJacobian);
    add_vector_form(new StabilizationResidual);
  }
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormContinuousGalerkin::InflowBoundaryJacobian::matrix_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* u, Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0;
  
  for (int i=0; i < n; i++)
    result += wt[i] * u->val[i] * v->val[i] * (-FlowField<Real>::dot_n(e,i)); // beta_dot_n < 0 for the inflow boundary
  
  return result;
}
template<typename Real, typename Scalar>
Scalar CustomWeakFormContinuousGalerkin::InflowBoundaryResidual::vector_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0;
  
  for (int i=0; i < n; i++)
    result += wt[i] * v->val[i] * (-FlowField<Real>::dot_n(e,i)) * (u_ext[0]->val[i] - inflow->value(e->x[i], e->y[i]));
  
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormContinuousGalerkin::VolumetricJacobian::matrix_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* u, Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0;
  for (int i=0; i < n; i++)
    result += wt[i] * v->val[i] * ( FlowField<Real>::dot_grad(u,e,i) + mu->value(e->x[i], e->y[i]) * u->val[i] );
  
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormContinuousGalerkin::VolumetricResidual::vector_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  Scalar result = 0;
  for (int i=0; i < n; i++)
    result += wt[i] * v->val[i] * ( FlowField<Real>::dot_grad(u_ext[0],e,i) + mu->value(e->x[i], e->y[i]) * u_ext[0]->val[i] );
  
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormContinuousGalerkin::StabilizationJacobian::matrix_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* u, Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  static FlowField<Real> beta;
  static ReactionTerm mu(SUPG);
 
  Scalar result = 0;
  Real norm_a_sq = 0.;
  Real norm_b_sq = 0.;
  for (int i=0; i < n; i++) 
  {
    result += wt[i] * beta.dot_grad(v,e,i) * ( beta.dot_grad(u,e,i) + mu.value(e->x[i], e->y[i]) * u->val[i] );
    norm_a_sq += 0.5 * wt[i] * sqr(beta.fn_a(e->x[i], e->y[i]));
    norm_b_sq += 0.5 * wt[i] * sqr(beta.fn_b(e->x[i], e->y[i]));
  }
  
  return result * sqr(e->diam)/(4*(norm_a_sq + norm_b_sq));
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormContinuousGalerkin::StabilizationResidual::vector_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  static FlowField<Real> beta;
  static ReactionTerm mu(SUPG);
  
  Scalar result = 0;
  Real norm_a_sq = 0.;
  Real norm_b_sq = 0.;
  for (int i=0; i < n; i++) 
  {
    result += wt[i] * beta.dot_grad(v,e,i) * ( beta.dot_grad(u_ext[0],e,i) + mu.value(e->x[i], e->y[i]) * u_ext[0]->val[i] );
    norm_a_sq += 0.5 * wt[i] * sqr(beta.fn_a(e->x[i], e->y[i]));
    norm_b_sq += 0.5 * wt[i] * sqr(beta.fn_b(e->x[i], e->y[i]));
  }
  
  return result * sqr(e->diam)/(4*(norm_a_sq + norm_b_sq));
}

//////////////////////////////////////////////////// DISCONTINUOUS APPROXIMATION ////////////////////////////////////////////////////

// Scalar average, vector jump.

#define AVG(w)      ( 0.5 * (w->get_val_central(i) + w->get_val_neighbor(i)) )

#define JUMP(w)     w->get_val_central(i)*e->nx[i] - w->get_val_neighbor(i)*e->nx[i],\
                    w->get_val_central(i)*e->ny[i] - w->get_val_neighbor(i)*e->ny[i] 

// Weak forms:

CustomWeakFormDiscontinuousGalerkin::CustomWeakFormDiscontinuousGalerkin(double theta) : WeakForm(1)
{
  add_matrix_form(new VolumetricJacobian);
  add_vector_form(new VolumetricResidual);
  add_matrix_form_surf(new BoundaryJacobian("nonzero_constant_inflow", new HermesFunction(1.0)));
  add_vector_form_surf(new BoundaryResidual("nonzero_constant_inflow", new HermesFunction(1.0)));
  add_matrix_form_surf(new BoundaryJacobian("zero_inflow", new HermesFunction(0.0)));
  add_vector_form_surf(new BoundaryResidual("zero_inflow", new HermesFunction(0.0)));
  add_matrix_form_surf(new BoundaryJacobian("varying_inflow", new InflowVariation));
  add_vector_form_surf(new BoundaryResidual("varying_inflow", new InflowVariation));
  add_matrix_form_surf(new BoundaryJacobian("outflow"));
  add_vector_form_surf(new BoundaryResidual("outflow"));
  add_matrix_form_surf(new InterfaceJacobian(theta));
  add_vector_form_surf(new InterfaceResidual(theta));
}

scalar CustomWeakFormDiscontinuousGalerkin::BoundaryJacobian::value(int n, double* wt, Func< scalar >* u_ext[], Func< double >* u, Func< double >* v, Geom< double >* e, ExtData< scalar >* ext) const
{
  scalar result = 0;
  
  for (int i = 0; i < n; i++) 
  {
    double beta_dot_n = FlowField<double>::dot_n(e,i);
    if (beta_dot_n >= 0)   // outflow
      result += wt[i] * u->val[i] * beta_dot_n * v->val[i];
  }
  
  return result;
}

Ord CustomWeakFormDiscontinuousGalerkin::BoundaryJacobian::ord(int n, double* wt, Func< Ord >* u_ext[], Func< Ord >* u, Func< Ord >* v, Geom< Ord >* e, ExtData< Ord >* ext) const
{
  return u->val[0] * FlowField<Ord>::dot_n(e,0) * v->val[0];
}

scalar CustomWeakFormDiscontinuousGalerkin::BoundaryResidual::value(int n, double* wt, Func< scalar >* u_ext[], Func< double >* v, Geom< double >* e, ExtData< scalar >* ext) const
{
  scalar result = 0;
  
  for (int i = 0; i < n; i++) 
  {
    double beta_dot_n = FlowField<double>::dot_n(e,i);
    if (beta_dot_n >= 0)   // outflow
      result += wt[i] * u_ext[0]->val[i] * beta_dot_n * v->val[i];
    else  // inflow
      result += wt[i] * (-beta_dot_n) * (-inflow->value(e->x[i], e->y[i])) * v->val[i];
  }
  
  return result;
}

Ord CustomWeakFormDiscontinuousGalerkin::BoundaryResidual::ord(int n, double* wt, Func< Ord >* u_ext[], Func< Ord >* v, Geom< Ord >* e, ExtData< Ord >* ext) const
{
  return FlowField<Ord>::dot_n(e,0) * (u_ext[0]->val[i] + inflow->value(e->x[0], e->y[0])) * v->val[0];
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::VolumetricJacobian::matrix_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* u, Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  static ReactionTerm mu(DG);
  
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * ( mu.value(e->x[i], e->y[i]) * v->val[i] - FlowField<Real>::dot_grad(v,e,i) );

  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::VolumetricResidual::vector_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  static ReactionTerm mu(DG);
  
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * u_ext[0]->val[i] * ( mu.value(e->x[i], e->y[i]) * v->val[i] - FlowField<Real>::dot_grad(v,e,i) );
  
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::InterfaceJacobian::matrix_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* u, Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  static FlowField<Real> beta;
  
  Scalar result = 0;
  for (int i = 0; i < n; i++) 
    result += wt[i] * ( AVG(u) * beta.dot(e->x[i], e->y[i], JUMP(v)) + 
                        theta * magn(beta.dot_n(e,i)) * dot2<Real>(JUMP(u), JUMP(v)) );
                        
  return result;
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormDiscontinuousGalerkin::InterfaceResidual::vector_form(int n, double* wt, Func< Scalar >* u_ext[], Func< Real >* v, Geom< Real >* e, ExtData< Scalar >* ext) const
{
  static FlowField<Real> beta;
  
  Scalar result = 0;
  for (int i = 0; i < n; i++) 
    result += wt[i] * ( AVG(u_ext[0]) * beta.dot(e->x[i], e->y[i], JUMP(v)) + 
                        theta * magn(beta.dot_n(e,i)) * dot2<Real>(JUMP(u_ext[0]), JUMP(v)) );
  
  return result;
}
