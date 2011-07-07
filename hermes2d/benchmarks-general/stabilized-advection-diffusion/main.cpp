//  This benchmark solves a linear advection-diffusion problem, using either the
//  standard continuous Galerkin method, one of its stabilized variants, or the 
//  discontinuous Galerkin method. 
//
//  PDE: div(\beta u - \epsilon \nabla u) = 0 where \beta = (B1, B2) is a constant vector.
//
//  Domain: Square (0, 1)x(0, 1).
//
//  BC:  Dirichlet, see the function scalar essential_bc_values() below.
//

#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace RefinementSelectors;

GalerkinMethod method = CG_STAB_SUPG;

const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements. 
const int INIT_BDY_REF_NUM = 2;                   // Number of initial refinements towards boundary. If INIT_BDY_REF_NUM == 0, 
                                                  // the first solution will be performed on a mesh (INIT_REF_NUM + 1) times 
                                                  // globally refined.
const int ORDER_INCREASE = 0;                     // Order increase for the refined space. If no change of order is allowed 
                                                  // (not even for computing the reference solution), set to 0.
const double THRESHOLD = 0.2;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = -1... do not refine.
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_H_ANISO;           // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // exact and the coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK. 

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;
  
  info("Discretization method: %s", method_names[method].c_str());
  
  if (method != CG && method != DG)
  {
    if (STRATEGY > -1 && CAND_LIST != H2D_H_ISO && CAND_LIST != H2D_H_ANISO)
      warning("The %s method should be used only with h-refinement.", method_names[method].c_str());
    
    int eff_order = (STRATEGY == -1) ? P_INIT : P_INIT + ORDER_INCREASE;
    
    if (method != CG_STAB_GLS_2)
    {
      if (eff_order != 1)
        warning("The %s method should be used only with 1st order elements.", method_names[method].c_str());
    }
    else
    {
      if (eff_order != 2)
        warning("The %s method should be used only with 2nd order elements.", method_names[method].c_str());
    }
  }
  
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square_quad.mesh", &mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("outflow", INIT_BDY_REF_NUM);
  //mesh.refine_towards_boundary("nonzero_Dirichlet", INIT_BDY_REF_NUM/2);

  // Create a space and refinement selector appropriate for the selected discretization method.
  Space *space;
  ProjBasedSelector *selector;
  ProjNormType norm;
  
  WeaklyImposableBC bc_fn(Hermes::vector<std::string>("nonzero_Dirichlet"), new NonzeroBoundaryValues(&mesh));
  WeaklyImposableBC bc_zero(Hermes::vector<std::string>("zero_Dirichlet", "outflow"), 0.0);
  EssentialBCs bcs(Hermes::vector<EssentialBoundaryCondition*>(&bc_fn, &bc_zero));
  
  // Initialize the weak formulation.
  WeakForm *wf;
  
  if (method != DG)
  {    
    space = new H1Space(&mesh, &bcs, P_INIT);
    selector = new L2ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
    norm = HERMES_L2_NORM;  // WARNING: In order to compare the errors with DG, L2 norm should be here.
    
    wf = new CustomWeakFormContinuousGalerkin(method, EPSILON);
  }
  else
  {
    space = new L2Space(&mesh, P_INIT);
    selector = new L2ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
    norm = HERMES_L2_NORM;
    // Disable weighting of refinement candidates.
    selector->set_error_weights(1, 1, 1);
    
    wf = new CustomWeakFormDiscontinuousGalerkin(bcs, EPSILON);
  }
  
  // Initialize coarse and reference mesh solution.
  Solution sln, ref_sln;
  
  // Set exact solution.
  CustomExactSolution exact(&mesh, EPSILON);
  
  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 440, 350));
  sview.fix_scale_width(50);
  sview.show_mesh(false);
  OrderView  oview("Polynomial orders", new WinGeom(450, 0, 400, 350));
  
  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu, graph_dof_exact, graph_cpu_exact;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  
  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  Space* actual_sln_space;
  do
  {
    info("---- Adaptivity step %d:", as);

    if (STRATEGY == -1)
      actual_sln_space = space;
    else
      // Construct globally refined reference mesh and setup reference space.
      actual_sln_space = Space::construct_refined_space(space, ORDER_INCREASE);

    int ndof_fine = Space::get_num_dofs(actual_sln_space);
    int ndof_coarse = Space::get_num_dofs(space);
    
    // Solve the linear system. If successful, obtain the solution.
    info("Solving on the refined mesh (%d NDOF).", ndof_fine);
    
    DiscreteProblem dp(wf, actual_sln_space);
    
    // Initial coefficient vector for the Newton's method.  
    scalar* coeff_vec = new scalar[ndof_fine];
    memset(coeff_vec, 0, ndof_fine * sizeof(scalar));
    
    // Perform Newton's iteration.
    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs)) 
      error("Newton's iteration failed.");
    Solution::vector_to_solution(solver->get_solution(), actual_sln_space, &ref_sln);
    
    // Calculate exact error.
    double err_exact_rel = hermes2d.calc_rel_error(&ref_sln, &exact, norm) * 100;
    info("ndof_fine: %d, err_exact_rel: %g%%", ndof_fine, err_exact_rel);

    // Time measurement.
    cpu_time.tick();
    
    // View the fine mesh solution and polynomial orders.
    sview.show(&ref_sln);
    oview.show(actual_sln_space);
    
    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    if (STRATEGY == -1) done = true;  // Do not adapt.
    else
    { 
      Adapt* adaptivity = new Adapt(space, norm);
      
      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection::project_global(space, &ref_sln, &sln, matrix_solver, norm); 
      
      // Calculate element errors and total error estimate.
      info("Calculating error estimate."); 
      bool solutions_for_adapt = true;
      double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt) * 100;

      // Report results.
      info("ndof_coarse: %d, err_est_rel: %g%%", ndof_coarse, err_est_rel);

      // Time measurement.
      cpu_time.tick();
      
      // Add entry to DOF and CPU convergence graphs.
      graph_dof.add_values(ndof_coarse, err_est_rel);
      graph_dof.save("conv_dof_est.dat");
      graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
      graph_cpu.save("conv_cpu_est.dat");
      graph_dof_exact.add_values(ndof_coarse, err_exact_rel);
      graph_dof_exact.save("conv_dof_exact.dat");
      graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel);
      graph_cpu_exact.save("conv_cpu_exact.dat");

      // Skip graphing time.
      cpu_time.tick(HERMES_SKIP);
      
      // If err_est too large, adapt the mesh.
      if (err_exact_rel < ERR_STOP) done = true;
      else 
      {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
  
        // Increase the counter of performed adaptivity steps.
        if (done == false)  as++;
      }
      if (Space::get_num_dofs(space) >= NDOF_STOP) done = true;
      
      // Clean up.
      if(done == false) 
      {
        delete actual_sln_space->get_mesh();
        delete actual_sln_space;
      }
      
      delete adaptivity;
    }
  }
  while (done == false);
  
  if (space != actual_sln_space) 
  {
    delete space;
    delete actual_sln_space->get_mesh();
  }
  delete actual_sln_space;
  delete solver;
  delete matrix;
  delete rhs;
  delete selector;
  
  verbose("Total running time: %g s", cpu_time.accumulated());
  
  // Wait for all views to be closed.
  View::wait();
  
  return 0;
}
