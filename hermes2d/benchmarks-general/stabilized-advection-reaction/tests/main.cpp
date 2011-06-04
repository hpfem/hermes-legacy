/** \addtogroup t_bench_stabilized-advection-reaction benchmarks/stabilized-advection-reaction
 *  \{
 *  \brief This test makes sure that the benchmark "stabilized-advection-reaction" works correctly.
 *
 *  \section s_params Parameters
 *   - P_INIT=0
 *   - INIT_REF=0
 *   - ORDER_INCREASE=1
 *   - THRESHOLD=0.2
 *   - STRATEGY=0
 *   - CAND_LIST=H2D_HP_ISO;
 *   - MESH_REGULARITY=-1
 *   - CONV_EXP=1.0
 *   - ERR_STOP=5.0
 *   - NDOF_STOP=8000
 *   - matrix_solver = SOLVER_UMFPACK
 *   - method = SUPG/DG
 *
 *  \section s_res Results
 *  - DOFs: 174 (SUPG) / 1229 (DG)
 *  - Iterations: 9 (SUPG) / 16 (DG)
 */

#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

// For writing results into file.
#include <fstream>
#include <iterator>

// The following two parameters will be set according to command line arguments.
DiscretizationMethod method;
int P_INIT;                                       // Initial polynomial degrees of mesh elements.

using namespace RefinementSelectors;

const bool SAVE_FINAL_SOLUTION = true;            // Save the final solution at specified points for comparison with the
                                                  // semi-analytic solution in Mathematica?
                                                  
const int INIT_REF = 0;                           // Number of initial uniform mesh refinements.
const double THRESHOLD = 0.20;                    // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = -1 ... do not refine.
                                                  // STRATEGY =  0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY =  1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY =  2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ISO;            // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                                  // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const int ORDER_INCREASE = 1;                     // Difference in polynomial orders of the coarse and the reference spaces.
const double ERR_STOP = 5.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 8000;                       // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).

// Construct a string representing the program options.
void make_str_from_program_options(std::stringstream& str);

int main(int argc, char* args[])
{
  std::string arg = (argc == 1) ? "SUPG" : args[1];
 
  if (arg == "SUPG")
  {
    method = SUPG;
    P_INIT = 1;
  }
  else if (arg == "DG")
  {
    method = DG;
    P_INIT = 0;
    error("DG option not yet supported.");
  }
  else
  {
    error("Invalid discretization method.");
  }
  
  // Instantiate a class with global functions.
  Hermes2D hermes2d;
  
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square.mesh", &mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF; i++) mesh.refine_all_elements();
  
  // Create a space and refinement selector appropriate for the selected discretization method.
  Space *space;
  ProjBasedSelector *selector;
  ProjNormType norm;
  
  if (method != DG)
  { 
    space = new H1Space(&mesh, P_INIT);
    norm = HERMES_L2_NORM;
    //norm = HERMES_H1_NORM;
    
    if (STRATEGY > -1)
    {
      selector = new H1ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
      selector->set_error_weights(1.15, 1.0, 1.414); // Prefer h-refinement over p-refinement.
      //selector = new L2ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
      //selector->set_error_weights(0.5, 2.0, 1.414); // Prefer h-refinement over p-refinement.
    }
  }
  else
  {
    space = new L2Space(&mesh, P_INIT);
    norm = HERMES_L2_NORM;
    
    if (STRATEGY > -1)
    {
      selector = new L2ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);   
      selector->set_error_weights(0.5, 2.0, 1.414); // Prefer h-refinement over p-refinemet.
    }
  }
  
  // Initialize the weak formulation.
  info("Discretization method: %s", method_names[method].c_str());
  
  WeakForm* wf;
  switch(method)
  {
    case CG:
      wf = new CustomWeakFormContinuousGalerkin(false);
      break;
    case SUPG:
      wf = new CustomWeakFormContinuousGalerkin(true);
      break; 
    case DG:
      wf = new CustomWeakFormDiscontinuousGalerkin();
      break;
  }      
  
  Solution sln;
  Solution ref_sln;
  
  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_ex, graph_cpu_ex, graph_dof_outfl, graph_cpu_outfl;
  
  // Time measurement.
  TimePeriod cpu_time, clk_time;
  cpu_time.tick();
  clk_time.tick();
  
  // Setup data structures for solving the discrete algebraic problem.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  
  // Initialize the preconditioner in the case of SOLVER_AZTECOO.
  if (matrix_solver == SOLVER_AZTECOO) 
  {
    ((AztecOOSolver*) solver)->set_solver(iterative_method);
    ((AztecOOSolver*) solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }
  
  // Load the exact solution evaluated at the Gauss-Kronrod quadrature points.
  SemiAnalyticSolution exact_sln("../exact/sol_GaussKronrod50.map");
  // Create a class for calculating a functional with known exact value (used for benchmarking).
  BoundaryIntegral benchmark_functional("outflow");
  
  
  int as = 1; bool done = false;
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
    
    // Initialize the FE problem.
    DiscreteProblem* dp = new DiscreteProblem(wf, actual_sln_space);
    
    // Speed-up assembly when 0-th order elements are used.
    bool is_fvm = true;
    Element *e;
    for_all_active_elements(e, actual_sln_space->get_mesh()) 
    {
      int order2d = actual_sln_space->get_element_order(e->id);
      if (H2D_GET_H_ORDER(order2d) > 0 || H2D_GET_V_ORDER(order2d) > 0)
      {
        is_fvm = false;
        break;
      }
    }
    if (is_fvm)
      dp->set_fvm();
    
    // Solve the linear system. If successful, obtain the solution.
    info("Solving on the refined mesh (%d NDOF).", ndof_fine);
    
    // Initial coefficient vector for the Newton's method.  
    scalar* coeff_vec = new scalar[ndof_fine];
    memset(coeff_vec, 0, ndof_fine * sizeof(scalar));
    
    // Perform Newton's iteration.
    if (!hermes2d.solve_newton(coeff_vec, dp, solver, matrix, rhs)) 
      error("Newton's iteration failed.");
    Solution::vector_to_solution(solver->get_solution(), actual_sln_space, &ref_sln);
    
    // Process the intermediate solution, but don't accumulate cpu time.
    cpu_time.tick();
    
    // Calculate relative error w.r.t. exact solution.
    info("Calculating relative L2 error w.r.t. exact solution.");
    double err_exact_rel = exact_sln.get_l2_rel_err(&ref_sln);
    info("Relative L2 error w.r.t. exact solution: %g%%", err_exact_rel*100);
    
    double outflow = benchmark_functional.value(&ref_sln);
    double exact_outflow = benchmark_functional.exact();
    double err_outflow = std::abs(outflow - exact_outflow)/exact_outflow;
    info("Integrated outflow: %f, relative error w.r.t. the ref. value (~%1.4f): %g%%", 
         outflow, exact_outflow, 100*err_outflow);
         
    cpu_time.tick(HERMES_SKIP);
         
    if (STRATEGY == -1) done = true;  // Do not adapt.
    else
    {
      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection::project_global(space, &ref_sln, &sln, matrix_solver, norm);  
          
      // Calculate element errors and total error estimate.
      info("Calculating error estimate."); 
      Adapt* adaptivity = new Adapt(space, norm);
      bool solutions_for_adapt = true;
      double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt, 
                          HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL);
                          
      int ndof_fine = Space::get_num_dofs(actual_sln_space);
      int ndof_coarse = Space::get_num_dofs(space);
  
      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
           ndof_coarse, ndof_fine, err_est_rel*100);

      // Get current time, use it in convergence graphs, and skip the graphing time.
      cpu_time.tick();
      double t = cpu_time.accumulated();

      // Add entry to DOF and CPU convergence graphs.
      graph_dof_est.add_values(ndof_fine, err_est_rel);
      graph_cpu_est.add_values(t, err_est_rel);
      graph_dof_ex.add_values(ndof_fine, err_exact_rel);
      graph_cpu_ex.add_values(t, err_exact_rel);
      graph_dof_outfl.add_values(ndof_fine, err_outflow);
      graph_cpu_outfl.add_values(t, err_outflow);
      
      cpu_time.tick(HERMES_SKIP);

      // If err_est_rel too large, adapt the mesh.
      if (err_est_rel*100 < ERR_STOP) done = true;
      else 
      {
        info("Adapting the coarse mesh.");
        done = adaptivity->adapt(selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

        if (Space::get_num_dofs(space) >= NDOF_STOP) 
        {
          done = true;
          break;
        }
      }

      // Clean up current adaptation step.
      delete adaptivity;
      if(done == false)
      {
        delete actual_sln_space->get_mesh();
        delete actual_sln_space;
      }
      
      as++;
    }
    
    delete dp;
  }
  while (done == false);

  // Solving finished - destroy the matrix solver and weak form objects.
  delete solver;
  delete matrix;
  delete rhs;
  delete wf;
  
  clk_time.tick();
  info("Total running time: %g s", clk_time.accumulated());
  info("Total cpu time: %g s", cpu_time.accumulated());
  
  // Wait for keyboard or mouse input.
  // View::wait();
  
  // Get the string summarizing the currently selected options.
  std::stringstream str;
  make_str_from_program_options(str);
  
  if (STRATEGY > -1)
  {
    // Save convergence graphs.   
    std::stringstream fdest, fcest, fdex, fcex, fdoutfl, fcoutfl;
    fdest << "conv_dof_est_" << str.str() << ".dat";
    fcest << "conv_cpu_est_" << str.str() << ".dat";
    fdex << "conv_dof_ex_" << str.str() << ".dat";
    fcex << "conv_cpu_ex_" << str.str() << ".dat";
    fdoutfl << "conv_dof_outfl_" << str.str() << ".dat";
    fcoutfl << "conv_cpu_outfl_" << str.str() << ".dat";
    
    graph_dof_est.save(fdest.str().c_str());
    graph_cpu_est.save(fcest.str().c_str());
    graph_dof_ex.save(fdex.str().c_str());
    graph_cpu_ex.save(fcex.str().c_str());
    graph_dof_outfl.save(fdoutfl.str().c_str());
    graph_cpu_outfl.save(fcoutfl.str().c_str());
  }
  
  if (SAVE_FINAL_SOLUTION)
  {
    // Save the final solution at 101x101 points of the domain for comparison with the semi-analytic solution.
    double y = 0.0, step = 0.01;
    int npts = int(1./step+0.5);
    double *res = new double [npts*npts];
    double *p = res;
    std::stringstream sssln;
    sssln << "sln_" << str.str() << ".dat";
    std::ofstream fs(sssln.str().c_str());
    info("Saving final solution to %s", sssln.str().c_str());
    for (int i = 0; i < npts; i++, y+=step)
    {
      std::cout << "."; std::cout.flush(); 
      double x = 0.0;
      for (int j = 0; j < npts; j++, x+=step)
        *p++ = ref_sln.get_pt_value(x, y);    
    }
    std::copy(res, res+npts*npts, std::ostream_iterator<double>(fs, "\n"));
    
    fs.close();
    delete [] res;
  }
  
  int ndof = Space::get_num_dofs(space);

  // Destroy spaces and refinement selector object.
  if (STRATEGY > -1)
  {
    delete space;
    delete actual_sln_space->get_mesh();
    delete selector;
  }
  delete actual_sln_space; 
  
  int n_dof_allowed = (method == DG) ? 1235 : 174;
  printf("\n n_dof_actual = %d\n", ndof);
  printf("n_dof_allowed = %d\n", n_dof_allowed);
  if (ndof <= n_dof_allowed) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

void make_str_from_program_options(std::stringstream& str)
{ 
  switch (method)
  {
    case CG:
      str << "cg";
      break;
    case SUPG:
      str << "supg";
      break;
    case DG:
      str << "dg";
      break;
  }
  
  if (STRATEGY > -1)
  {
    switch (CAND_LIST) 
    {
      case H2D_H_ANISO:
      case H2D_H_ISO:
        str << "_h" << P_INIT;
        break;
      case H2D_P_ANISO:
      case H2D_P_ISO:
        str << "_p" << INIT_REF;
        break;
      default:
        str << "_hp";
        break;
    }
    switch (CAND_LIST) 
    {
      case H2D_H_ANISO:
      case H2D_P_ANISO:
      case H2D_HP_ANISO:
        str << "_aniso";
        break;
      case H2D_H_ISO:
      case H2D_P_ISO:
      case H2D_HP_ISO:
        str << "_iso";
        break;
      case H2D_HP_ANISO_H:
        str << "_anisoh";
        break;
      case H2D_HP_ANISO_P:
        str << "_anisop";
        break;
    }
  }
  else
  {
    str << "_h-" << INIT_REF << "_p-" << P_INIT;
  }
}
