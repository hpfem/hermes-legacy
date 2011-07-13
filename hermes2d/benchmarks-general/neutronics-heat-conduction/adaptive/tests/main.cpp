#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../../definitions.h"

using namespace RefinementSelectors;

// This test makes sure that the benchmark "neutronics-heat-conduction-adapt" works correctly.

using namespace RefinementSelectors;

const bool SOLVE_ON_COARSE_MESH = false;   // true... The discrete problem is solved on coarse meshes in every adaptivity step 
                                           //         to obtain the coarse solution for error computation.
                                           // false...The discrete problem is always solved on the refined meshes; second
                                           //         solution for error computation is obtained by orthogonal projection onto 
                                           //         the coarse meshes.
const int INIT_GLOB_REF_NUM = 1;           // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;            // Number of initial refinements towards boundary.
const int P_INIT = 1;                      // Initial polynomial degree of all mesh elements.

// Time-stepping:
const double TIME_STEP = 0.1;              // Time step in seconds.
const double T_FINAL = 1;                  // Time interval length (secs.).

// Adaptivity:
const int UNREF_FREQ = 5;                  // Every UNREF_FREQ time step the mesh is unrefined.
const int UNREF_METHOD = 3;                // 1... mesh reset to basemesh and poly degrees to P_INIT.   
                                           // 2... one ref. layer shaved off, poly degrees reset to P_INIT.
                                           // 3... one ref. layer shaved off, poly degrees decreased by one. 
const double THRESHOLD = 0.3;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ISO;     // Predefined list of element refinement candidates. Possible values are
                                           // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                           // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                           // See User Documentation for details.
const int MESH_REGULARITY = -1;            // Maximum allowed level of hanging nodes:
                                           // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                           // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                           // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                           // Note that regular meshes are not supported, this is due to
                                           // their notoriously bad performance.
const double CONV_EXP = 1.0;               // Default value is 1.0. This parameter influences the selection of
                                           // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1;                 // Stopping criterion for hp-adaptivity
                                           // (relative error between reference and coarse solution in percent).
const int NDOF_STOP = 100000;              // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.
const int ORDER_INCREASE = 1;              // Increase in approximation order associated with the global refinement.


// Newton's method:
const double NEWTON_TOL_COARSE = 1.0e-6;   // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 5.0e-6;     // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 100;           // Maximum allowed number of Newton iterations.

/* 
  Time-stepping method. Possibilities: 
    Explicit methods:
      Explicit_RK_1, Explicit_RK_2, Explicit_RK_3, Explicit_RK_4.
    Implicit methods: 
      Implicit_RK_1, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, Implicit_ESIRK_2_2, Implicit_SDIRK_2_2, 
      Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Implicit_Lobatto_IIIA_3_4, 
      Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_5_4.
    Embedded explicit methods:
      Explicit_HEUN_EULER_2_12_embedded, Explicit_BOGACKI_SHAMPINE_4_23_embedded, Explicit_FEHLBERG_6_45_embedded,
      Explicit_CASH_KARP_6_45_embedded, Explicit_DORMAND_PRINCE_7_45_embedded.
    Embedded implicit methods:
      Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded, 
      Implicit_SDIRK_BILLINGTON_3_23_embedded, Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded, 
      Implicit_DIRK_ISMAIL_7_45_embedded.
*/
ButcherTableType butcher_table_type = Implicit_Crank_Nicolson_2_2;

/* 
  Linear algebraic system solver for the coarse and refined problems, respectively. Possibilities: 
     SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK. 
*/
MatrixSolverType matrix_solver_coarse = SOLVER_UMFPACK;  
MatrixSolverType matrix_solver_fine = SOLVER_UMFPACK;

int main(int argc, char* argv[])
{
  Hermes2D hermes2d;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
  
  info("Simulation length = %f s, number of time steps: %d", T_FINAL, int(T_FINAL / TIME_STEP + 0.5));
  
  // Load the mesh.
  Mesh basemesh, mesh_T, mesh_phi;
  H2DReader mloader;
  mloader.load("../../domain.mesh", &basemesh);

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_GLOB_REF_NUM; i++)
    basemesh.refine_all_elements(0, true);
  basemesh.refine_towards_boundary("zero_Dirichlet", INIT_BDY_REF_NUM, true, true);
  
  const int LX = 100, LY = 100; // Domain extents in x- and y- directions.
  
  // Create a special mesh for each physical field.
  mesh_T.copy(&basemesh);
  mesh_phi.copy(&basemesh);
 
  // Solutions in the previous time step (converging within the time stepping loop).
  // Initially set to the exact distributions at time t = 0.
  TemperatureField::ExactDistribution T_prev_time(&mesh_T, LX, LY);
  NeutronField::ExactDistribution phi_prev_time(&mesh_phi, LX, LY);
  Hermes::vector<Solution*> prev_time_solutions(&T_prev_time, &phi_prev_time);
  
  // Solutions on the coarse and refined meshes in current time step (converging within the Newton's loop).
  Solution T_coarse, phi_coarse, T_fine, phi_fine;
  Hermes::vector<Solution*> coarse_mesh_solutions(&T_coarse, &phi_coarse);
  Hermes::vector<Solution*> fine_mesh_solutions(&T_fine, &phi_fine);
  
  // Exact solutions for error evaluation.
  TemperatureField::ExactDistribution T_exact_solution(&mesh_T, LX, LY);
  NeutronField::ExactDistribution phi_exact_solution(&mesh_phi, LX, LY);
  Hermes::vector<Solution*> exact_solutions(&T_exact_solution, &phi_exact_solution);
  
  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_essential("zero_Dirichlet", 0.0);
  EssentialBCs bcs(&bc_essential);
  
  // Create H1 spaces with default shapesets.
  H1Space space_T(&mesh_T, &bcs, P_INIT);
  H1Space space_phi(&mesh_phi, &bcs, P_INIT);
  Hermes::vector<Space*> spaces(&space_T, &space_phi);
  
  // Initialize the weak formulation.
  CustomWeakForm wf(LX, LY);
  
  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
  
  // Initialize Runge-Kutta time stepping.
  
  double current_time = 0;
    
  // Choose the Butcher's table for to the selected method.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());
   
  // Run the time loop.  
  int ts = 1;
  do 
  {
    // Periodic global derefinement.
    
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      switch (UNREF_METHOD) 
      {
        case 1: 
          mesh_T.copy(&basemesh);
          mesh_phi.copy(&basemesh);
          space_T.set_uniform_order(P_INIT);
          space_phi.set_uniform_order(P_INIT);
          break;
        case 2: 
          mesh_T.unrefine_all_elements();
          mesh_phi.unrefine_all_elements();
          space_T.set_uniform_order(P_INIT);
          space_phi.set_uniform_order(P_INIT);
          break;
        case 3: 
          mesh_T.unrefine_all_elements();
          mesh_phi.unrefine_all_elements();
          //space_T.adjust_element_order(-1, P_INIT);
          space_T.adjust_element_order(-1, -1, P_INIT, P_INIT);
          //space_phi.adjust_element_order(-1, P_INIT);
          space_phi.adjust_element_order(-1, -1, P_INIT, P_INIT);
          break;
        default: 
          error("Wrong global derefinement method.");
      }
    }

    // Adaptivity loop:
    bool done = false;
    int as = 0;
    do 
    {
      as++;
      
      info("---- Time step %d, adaptivity step %d:", ts, as);
      
      // Construct globally refined reference mesh and setup reference space.
      Hermes::vector<Space *> *ref_spaces = Space::construct_refined_spaces(spaces, ORDER_INCREASE);

      // Initialize discrete problem on reference mesh.
      DiscreteProblem dp_fine(&wf, *ref_spaces);
      
      // Initialize a Runge-Kutta solver to obtain the approximation of next time level solution
      // on the refined mesh.
      RungeKutta runge_kutta_fine(&dp_fine, &bt, matrix_solver_fine);
      
      // Perform one Runge-Kutta time step according to the selected Butcher's table.
      info("Runge-Kutta time step on refined meshes (t = %g s, tau = %g s, stages: %d).", current_time, TIME_STEP, bt.get_size());
      
      bool verbose = true;
      bool reassemble_jacobian = true;
      if (!runge_kutta_fine.rk_time_step(current_time, TIME_STEP, prev_time_solutions, fine_mesh_solutions, 
          reassemble_jacobian, verbose, NEWTON_TOL_FINE, NEWTON_MAX_ITER)) 
      {
        error("Runge-Kutta time step failed, try to decrease time step size.");
      }
     
      if (SOLVE_ON_COARSE_MESH) 
      {        
        // Initialize discrete problem on the coarse meshes.
        DiscreteProblem dp_coarse(&wf, spaces);

        // Initialize a Runge-Kutta solver to obtain the approximation of next time level solution
        // on the coarse mesh prior refinement.
        RungeKutta runge_kutta_coarse(&dp_coarse, &bt, matrix_solver_coarse);
        
        // Perform one Runge-Kutta time step according to the selected Butcher's table.
        info("Runge-Kutta time step on coarse meshes (t = %g s, tau = %g s, stages: %d).", current_time, TIME_STEP, bt.get_size());
        
        bool verbose = true;
        bool reassemble_jacobian = true;
        if (!runge_kutta_coarse.rk_time_step(current_time, TIME_STEP, prev_time_solutions, coarse_mesh_solutions, 
            reassemble_jacobian, verbose, NEWTON_TOL_COARSE, NEWTON_MAX_ITER)) 
        {
          error("Runge-Kutta time step failed, try to decrease time step size.");
        }
      } 
      else 
      {
        // Projection onto the new coarse meshes.
        info("Projecting fine mesh solutions back onto coarse meshes.");
        OGProjection::project_global(spaces, fine_mesh_solutions, coarse_mesh_solutions, matrix_solver_coarse); 
      }

      // Calculate element errors.
      info("Calculating error estimate and exact error."); 
      Adapt* adaptivity = new Adapt(spaces);

      // Calculate error estimate for each solution component and the total error estimate.
      Hermes::vector<double> err_est_rel;
      double err_est_rel_total = adaptivity->calc_err_est(coarse_mesh_solutions, fine_mesh_solutions, &err_est_rel) * 100;

      info("T: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%", 
            space_T.get_num_dofs(), (*ref_spaces)[0]->get_num_dofs(), err_est_rel[0]*100);
      info("phi: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%", 
            space_phi.get_num_dofs(), (*ref_spaces)[1]->get_num_dofs(), err_est_rel[1]*100);
       
      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) done = true;
      else 
      {
        info("Adapting the coarse meshes.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector*> (&selector, &selector), THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (Space::get_num_dofs(spaces) >= NDOF_STOP) done = true; 
      }
      
      delete adaptivity;
      
      for (unsigned int i = 0; i < ref_spaces->size(); i++)
      {
        delete ref_spaces->at(i); // Mesh dynamically allocated for each space is held by the corresponding reference solution.
        if (!done) delete fine_mesh_solutions[i]->get_mesh();
      }
      delete ref_spaces;
      
    }
    while (!done);
    
    // Set the time indicator to that for which rk_time_step returned the new solution.
    current_time += TIME_STEP;
    
    // Update time-dependent exact solutions.
    T_exact_solution.update(current_time);
    phi_exact_solution.update(current_time);
    
    // Calculate exact error.
    info("Calculating H1 error between the latest fine mesh solution and the exact solution in time t = %3.2f s.", current_time);
    
    Hermes::vector<double> exact_errors;
    exact_errors.push_back(hermes2d.calc_rel_error(&T_fine, &T_exact_solution, HERMES_H1_NORM));
    exact_errors.push_back(hermes2d.calc_rel_error(&phi_fine, &phi_exact_solution, HERMES_H1_NORM));
    
    double maxerr = std::max(exact_errors[0], exact_errors[1])*100;
    info("Exact solution error for T (H1 norm): %g %%", exact_errors[0]*100);
    info("Exact solution error for phi (H1 norm): %g %%", exact_errors[1]*100);
    info("Exact solution error (maximum): %g %%", maxerr);
          
    // Make the fine mesh solution at current time level the previous time level solution in the following time step.
    T_prev_time.copy(&T_fine);
    phi_prev_time.copy(&phi_fine);
    
    // Delete the refined meshes associated with fine_mesh_solutions. Although this makes
    // the meshes in prev_time_solutions invalid, these solutions will be then projected
    // on the new meshes by rk_time_step and take the new meshes. Without this deallocation,
    // several MBs of memory would be consumed every time-step, which would eventually become
    // a major memory leak.
    delete T_fine.get_mesh();
    delete phi_fine.get_mesh();
    
    ts++;
  }
  while (fabs(current_time - T_FINAL) > 1e-12);
    
  cpu_time.tick();
  verbose("Total running time: %g s", cpu_time.accumulated());

  info("Coordinate ( 0.0,  0.0) T value = %lf", T_prev_time.get_pt_value(0.0, 0.0));
  info("Coordinate ( 25.0, 25.0) T value = %lf", T_prev_time.get_pt_value(25.0, 25.0));
  info("Coordinate ( 25.0, 75.0) T value = %lf", T_prev_time.get_pt_value(25.0, 75.0));
  info("Coordinate ( 75.0, 25.0) T value = %lf", T_prev_time.get_pt_value(75.0, 25.0));
  info("Coordinate ( 75.0, 75.0) T value = %lf", T_prev_time.get_pt_value(75.0, 75.0));

  info("Coordinate ( 0.0,  0.0) phi value = %lf", phi_prev_time.get_pt_value(0.0, 0.0));
  info("Coordinate ( 25.0, 25.0) phi value = %lf", phi_prev_time.get_pt_value(25.0, 25.0));
  info("Coordinate ( 25.0, 75.0) phi value = %lf", phi_prev_time.get_pt_value(25.0, 75.0));
  info("Coordinate ( 75.0, 25.0) phi value = %lf", phi_prev_time.get_pt_value(75.0, 25.0));
  info("Coordinate ( 75.0, 75.0) phi value = %lf", phi_prev_time.get_pt_value(75.0, 75.0));
  
  int success = 1;
  
  // Test convergence to correct results.  
  double eps = 1e-5;
  if (fabs(T_prev_time.get_pt_value(0.0, 0.0) - 0.000000) > eps) {
    info("Coordinate (  0,  0) T value = %lf\n", T_prev_time.get_pt_value(0.0, 0.0));
    success = 0;
  }
  
  if (fabs(T_prev_time.get_pt_value(25.0, 25.0) - 0.880524) > eps) {
    printf("Coordinate ( 25, 25) T value = %lf\n", T_prev_time.get_pt_value(25.0, 25.0));
    success = 0;
  }
  
  if (fabs(T_prev_time.get_pt_value(25.0, 75.0) - 0.880524) > eps) {
    printf("Coordinate ( 25, 75) T value = %lf\n", T_prev_time.get_pt_value(25.0, 75.0));
    success = 0;
  }
  
  if (fabs(T_prev_time.get_pt_value(75.0, 25.0) - 0.880524) > eps) {
    printf("Coordinate ( 75, 25) T value = %lf\n", T_prev_time.get_pt_value(75.0, 25.0));
    success = 0;
  }
  
  if (fabs(T_prev_time.get_pt_value(75.0, 75.0) - 0.880524) > eps) {
    printf("Coordinate ( 75, 75) T value = %lf\n", T_prev_time.get_pt_value(75.0, 75.0));
    success = 0;
  }
  
  if (fabs(phi_prev_time.get_pt_value(0.0, 0.0) - 0.000000) > eps) {
    printf("Coordinate (  0,  0) phi value = %lf\n", phi_prev_time.get_pt_value(0.0, 0.0));
    success = 0;
  }
  
  if (fabs(phi_prev_time.get_pt_value(25.0, 25.0) - 0.071353) > eps) {
    printf("Coordinate ( 25, 25) phi value = %lf\n", phi_prev_time.get_pt_value(25.0, 25.0));
    success = 0;
  }
  
  if (fabs(phi_prev_time.get_pt_value(25.0, 75.0) - 0.214066) > eps) {
    printf("Coordinate ( 25, 75) phi value = %lf\n", phi_prev_time.get_pt_value(25.0, 75.0));
    success = 0;
  }
  
  if (fabs(phi_prev_time.get_pt_value(75.0, 25.0) - 0.214066) > eps) {
    printf("Coordinate ( 75, 25) phi value = %lf\n", phi_prev_time.get_pt_value(75.0, 25.0));
    success = 0;
  }
  
  if (fabs(phi_prev_time.get_pt_value(75.0, 75.0) - 0.642227) > eps) {
    printf("Coordinate ( 75, 75) phi value = %lf\n", phi_prev_time.get_pt_value(75.0, 75.0));
    success = 0;
  }
  
  // Test adaptivity.
  int ndof_allowed_T = 30;
  int ndof_allowed_phi = 50;
  int ndof_T = Space::get_num_dofs(&space_T);
  int ndof_phi = Space::get_num_dofs(&space_phi);
  printf("ndof_actual_T = %d\n", ndof_T);
  printf("ndof_actual_phi = %d\n", ndof_phi);
  printf("ndof_allowed_T = %d\n", ndof_allowed_T);
  printf("ndof_allowed_phi = %d\n", ndof_allowed_phi);
  if ((ndof_T > ndof_allowed_T) || (ndof_phi > ndof_allowed_phi)) {
    // ndofs_T was 25 and ndof_phi was 46 at the time this test was created
    printf("Adaptivity failed.");
    success = 0;
  }
  
  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
