#define HERMES_REPORT_ALL

////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "definitions.h"

CustomWeakForm::CustomWeakForm( const MaterialPropertyMaps& matprop,
                                Hermes::vector<Solution*>& iterates,
                                double init_keff, std::string bdy_vacuum )
  : DefaultWeakFormSourceIteration(matprop, iterates, init_keff, HERMES_AXISYM_Y)
{
  for (unsigned int g = 0; g < matprop.get_G(); g++)
  {
    add_matrix_form_surf(new VacuumBoundaryCondition::Jacobian(g, bdy_vacuum, HERMES_AXISYM_Y));
    add_vector_form_surf(new VacuumBoundaryCondition::Residual(g, bdy_vacuum, HERMES_AXISYM_Y));
  }
}

scalar ErrorForm::value(int n, double *wt, Func<scalar> *u_ext[],
                        Func<scalar> *u, Func<scalar> *v, Geom<double> *e,
                        ExtData<scalar> *ext) const
{
  switch (projNormType)
  {
    case HERMES_L2_NORM:
      return l2_error_form_axisym<double, scalar>(n, wt, u_ext, u, v, e, ext);
    case HERMES_H1_NORM:
      return h1_error_form_axisym<double, scalar>(n, wt, u_ext, u, v, e, ext);
    default:
      error("Only the H1 and L2 norms are currently implemented.");
      return 0.0;
  }
}

Ord ErrorForm::ord(int n, double *wt, Func<Ord> *u_ext[],
                   Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                   ExtData<Ord> *ext) const
{
  switch (projNormType)
  {
    case HERMES_L2_NORM:
      return l2_error_form_axisym<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    case HERMES_H1_NORM:
      return h1_error_form_axisym<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    default:
      error("Only the H1 and L2 norms are currently implemented.");
      return Ord();
  }
}

// Integral over the active core.
double integrate(MeshFunction* sln, std::string area)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  
  double integral = 0.0;
  Element* e;
  Mesh *mesh = sln->get_mesh();
  int marker = mesh->get_element_markers_conversion().get_internal_marker(area);
  
  for_all_active_elements(e, mesh)
  {
    if (e->marker == marker)
    {
      update_limit_table(e->get_mode());
      sln->set_active_element(e);
      RefMap* ru = sln->get_refmap();
      int o = sln->get_fn_order() + ru->get_inv_ref_order();
      limit_order(o);
      sln->set_quad_order(o, H2D_FN_VAL);
      scalar *uval = sln->get_fn_values();
      double* x = ru->get_phys_x(o);
      double result = 0.0;
      h1_integrate_expression(x[i] * uval[i]);
      integral += result;
    }
  }
  
  return 2.0 * M_PI * integral;
}

// Calculate number of negative solution values.
int get_num_of_neg(MeshFunction *sln)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  Element* e;
  Mesh* mesh = sln->get_mesh();
  
  int n = 0;
  
  for_all_active_elements(e, mesh)
  {
    update_limit_table(e->get_mode());
    sln->set_active_element(e);
    RefMap* ru = sln->get_refmap();
    int o = sln->get_fn_order() + ru->get_inv_ref_order();
    limit_order(o);
    sln->set_quad_order(o, H2D_FN_VAL);
    scalar *uval = sln->get_fn_values();
    int np = quad->get_num_points(o);
    
    for (int i = 0; i < np; i++)
      if (uval[i] < -1e-12)
        n++;
  }
  
  return n;
}

int power_iteration(const Hermes2D& hermes2d, const MaterialPropertyMaps& matprop, 
                    const Hermes::vector<Space *>& spaces, DefaultWeakFormSourceIteration* wf, 
                    const Hermes::vector<Solution *>& solutions, const std::string& fission_region, 
                    double tol, SparseMatrix *mat, Vector* rhs, Solver *solver)
{
  // Sanity checks.
  if (spaces.size() != solutions.size()) 
    error("Spaces and solutions supplied to power_iteration do not match."); 
 
  // Number of energy groups.
  int G = spaces.size();
  
  // Initialize the discrete problem.
  DiscreteProblem dp(wf, spaces);
  int ndof = Space::get_num_dofs(spaces);
    
  // The following variables will store pointers to solutions obtained at each iteration and will be needed for 
  // updating the eigenvalue. 
  Hermes::vector<Solution*> new_solutions;
  for (int g = 0; g < G; g++) 
    new_solutions.push_back(new Solution(solutions[g]->get_mesh()));
  
  // This power iteration will most probably run on a different mesh than the previous one and so will be different
  // the corresponding algebraic system. We will need to factorize it anew (but then, the L and U factors may be 
  // reused until the next adaptation changes the mesh again).
  // TODO: This could be solved more elegantly by defining a function Solver::reinit().
  solver->set_factorization_scheme(HERMES_FACTORIZE_FROM_SCRATCH);
  
  // Initial coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[ndof];
  
  // Force the Jacobian assembling in the first iteration.
  bool Jacobian_changed = true;
  
  bool eigen_done = false; int it = 0;
  do 
  {
    memset(coeff_vec, 0.0, ndof*sizeof(scalar));

    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, mat, rhs, Jacobian_changed)) 
      error("Newton's iteration failed.");
    
    // The matrix doesn't change within the power iteration loop, so it does not need to be reassembled again and 
    // the first computed LU factorization may be completely reused in following iterations.
    Jacobian_changed = false;
    solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);
    
    // Convert coefficients vector into a set of Solution pointers.
    Solution::vector_to_solutions(solver->get_solution(), spaces, new_solutions);

    // Update fission sources.
    using WeakFormsNeutronics::Multigroup::SupportClasses::Common::SourceFilter;
    SourceFilter new_source(new_solutions, &matprop, fission_region);
    SourceFilter old_source(solutions, &matprop, fission_region);

    // Compute the eigenvalue for current iteration.
    double k_new = wf->get_keff() * (integrate(&new_source, fission_region) / integrate(&old_source, fission_region));

    info("      dominant eigenvalue (est): %g, rel. difference: %g", k_new, fabs((wf->get_keff() - k_new) / k_new));

    // Stopping criterion.
    if (fabs((wf->get_keff() - k_new) / k_new) < tol) eigen_done = true;

    // Update the final eigenvalue.
    wf->update_keff(k_new);

    it++;
        
    // Store the new eigenvector approximation in the result.
    for (int g = 0; g < G; g++) 
      solutions[g]->copy(new_solutions[g]); 
  }
  while (!eigen_done);
  
  // Free memory.
  for (int g = 0; g < G; g++) 
    delete new_solutions[g];
  
  return it;
}
