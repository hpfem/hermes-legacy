#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

using namespace RefinementSelectors;

/** \addtogroup t_bench_stabilized-advection-diffusion benchmarks-general/stabilized-advection-diffusion
 *  \{
 *  \brief This test makes sure that the benchmark "stabilized-advection-diffusion" works correctly.
 *
 *  \section s_params Common parameters
 *   - INIT_REF_NUM = 1
 *   - INIT_BDY_REF_NUM = 2
 *   - THRESHOLD=0.2
 *   - STRATEGY=0
 *   - MESH_REGULARITY=-1
 *   - CONV_EXP=1.0
 *   - ERR_STOP=5.0
 *   - NDOF_STOP=60000
 *   - matrix_solver = SOLVER_UMFPACK
 *
 *  \section s_res Results
 *   - DOFs:
 *      - CG:       5
 *      - SUPG:     184
 *      - GLS_1:    187
 *      - GLS_2:    75
 *      - SGS:      417
 *      - SGS_ALT:  178
 *      - DG:       TODO
 *   - Adaptivity steps: 
 *      - CG:       2
 *      - SUPG:     10
 *      - GLS_1:    10
 *      - GLS_2:    3
 *      - SGS:      13
 *      - SGS_ALT:  10  
 *      - DG:       TODO
 */

const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements. 
const int INIT_BDY_REF_NUM = 2;                   // Number of initial refinements towards boundary. If INIT_BDY_REF_NUM == 0, 
                                                  // the first solution will be performed on a mesh (INIT_REF_NUM + 1) times 
                                                  // globally refined.
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
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 3.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // exact and the coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,

int main(int argc, char* argv[])
{ 
  GalerkinMethod method;
  CandList CAND_LIST;
  int P_INIT;                             
  int ORDER_INCREASE;
  int n_dof_allowed;
  int iter_allowed;
  
  std::string arg = (argc == 1) ? "CG" : argv[1];
  
  if (arg == "CG")
  {
    method = CG;
    CAND_LIST = H2D_HP_ANISO;
    P_INIT = 1;
    ORDER_INCREASE = 1;
    n_dof_allowed = 7;
    iter_allowed = 3;
  }
  else if (arg == "SUPG")
  {
    method = CG_STAB_SUPG;
    CAND_LIST = H2D_H_ANISO;
    P_INIT = 1;
    ORDER_INCREASE = 0;
    n_dof_allowed = 190;
    iter_allowed = 11;
  }
  else if (arg == "GLS_1")
  {
    method = CG_STAB_GLS_1;
    CAND_LIST = H2D_H_ANISO;
    P_INIT = 1;
    ORDER_INCREASE = 0;
    n_dof_allowed = 195;
    iter_allowed = 12;
  }
  else if (arg == "GLS_2")
  {
    method = CG_STAB_GLS_2;
    CAND_LIST = H2D_H_ANISO;
    P_INIT = 2;
    ORDER_INCREASE = 0;
    n_dof_allowed = 80;
    iter_allowed = 4;
  }
  else if (arg == "SGS")
  {
    method = CG_STAB_SGS;
    CAND_LIST = H2D_H_ANISO;
    P_INIT = 1;
    ORDER_INCREASE = 0;
    n_dof_allowed = 430;
    iter_allowed = 14;
  }
  else if (arg == "SGS_ALT")
  {
    method = CG_STAB_SGS_ALT;
    CAND_LIST = H2D_H_ANISO;
    P_INIT = 1;
    ORDER_INCREASE = 0;
    n_dof_allowed = 190;
    iter_allowed = 11;
  }
  else if (arg == "DG")
  {
    method = DG;
    CAND_LIST = H2D_HP_ANISO;
    P_INIT = 0;
    ORDER_INCREASE = 1;
    n_dof_allowed = 0;//TODO
    iter_allowed = 0;//TODO
    error("DG option not yet supported.");
  }
  else
  {
    error("Invalid discretization method.");
  }
  
  // Instantiate a class with global functions.
  Hermes2D hermes2d;
  
  info("Discretization method: %s", method_names[method].c_str());
  
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square_quad.mesh", &mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("outflow", INIT_BDY_REF_NUM);
  //mesh.refine_towards_boundary(NONZERO_DIRICHLET, INIT_BDY_REF_NUM/2);

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
  
  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu, graph_dof_exact, graph_cpu_exact;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Setup data structures for solving the discrete algebraic problem.
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
  
  int ndof = Space::get_num_dofs(space);

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

  // Test number of DOF.
  printf("n_dof_actual = %d\n", ndof);
  printf("n_dof_allowed = %d\n", n_dof_allowed);
  if (ndof <= n_dof_allowed) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
  
  // Test number of iterations.
  printf("iterations = %d\n", as);
  printf("iterations allowed = %d\n", iter_allowed);
  if (as <= iter_allowed) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
