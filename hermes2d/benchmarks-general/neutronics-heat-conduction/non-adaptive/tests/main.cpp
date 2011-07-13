#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../../definitions.h"

const int INIT_GLOB_REF_NUM = 2;                    // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;                     // Number of initial refinements towards boundary.
const int P_INIT = 2;                               // Initial polynomial degree of all mesh elements.
const double NEWTON_TOL = 1e-6;                     // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                    // Maximum allowed number of Newton iterations.
const double TIME_STEP = 0.1;                       // Time step in seconds.
const double T_FINAL = 10.0;                        // Time interval length (secs.).

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
ButcherTableType butcher_table_type = Implicit_SDIRK_2_2; 

/* 
Linear algebraic system solver. Possibilities: 
SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK. 
*/
MatrixSolverType matrix_solver = SOLVER_UMFPACK;    

int main(int argc, char* argv[])
{
  Hermes2D hermes2d;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
  
  info("Simulation length = %d s, number of time steps: %f", T_FINAL, int(T_FINAL / TIME_STEP + 0.5));
  
  // Load the mesh file.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../../domain.mesh", &mesh);
  
  const int LX = 100, LY = 100; // Domain extents in x- and y- directions.
  
  // Perform initial mesh refinements.
  for (int i=0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("zero_Dirichlet", INIT_BDY_REF_NUM);
  
  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_essential("zero_Dirichlet", 0.0);
  EssentialBCs bcs(&bc_essential);
  
  // Create H1 spaces with default shapesets.
  H1Space space_T(&mesh, &bcs, P_INIT);
  H1Space space_phi(&mesh, &bcs, P_INIT);
  Hermes::vector<Space*> spaces(&space_T, &space_phi);
  
  // Exact solutions for error evaluation.
  TemperatureField::ExactDistribution T_exact_solution(&mesh, LX, LY);
  NeutronField::ExactDistribution phi_exact_solution(&mesh, LX, LY);
  Hermes::vector<Solution*> sln_exact(&T_exact_solution, &phi_exact_solution);
  
  // Solutions in the previous time step.
  Solution T_prev_time, phi_prev_time;
  Hermes::vector<Solution*> sln_time_prev(&T_prev_time, &phi_prev_time);
  
  // Solutions in the current time.
  Solution T_current_time, phi_current_time;
  Hermes::vector<Solution*> sln_time_new(&T_current_time, &phi_current_time);
    
  // Initialize the weak formulation.
  CustomWeakForm wf(LX, LY);
  
  // Initialize the FE problem.
  DiscreteProblem dp(&wf, spaces);
  
  // Initialize Runge-Kutta time stepping.
  
  double current_time = 0;
  T_exact_solution.update(current_time);
  phi_exact_solution.update(current_time);
  
  // Use the exact solution in time t = 0s as initial condition (could represent
  // e.g. the results of initial steady state calculation).
  info("Projecting the initial condition on computation mesh.");
  OGProjection::project_global(spaces, sln_exact, sln_time_prev, matrix_solver);
  
  // Choose the Butcher's table for to the selected method.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());
  
  RungeKutta runge_kutta(&dp, &bt, matrix_solver);
  // Although the Jacobian changes every iteration due to the nonlinearity in Sigma_r and lambda,
  // the change is not too big and some information from its first factorization may be reused.
  bool reassemble_jacobian = true;
  runge_kutta.get_matrix_solver()->set_factorization_scheme(HERMES_REUSE_MATRIX_REORDERING_AND_SCALING);
  
  // Run the time loop.  
  do {    
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, tau = %g s, stages: %d).", current_time, TIME_STEP, bt.get_size());
    
    bool verbose = false;
    if (!runge_kutta.rk_time_step(current_time, TIME_STEP, sln_time_prev, sln_time_new, 
      reassemble_jacobian, verbose, NEWTON_TOL, NEWTON_MAX_ITER)) 
    {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }
    
    // Set the time indicator to that for which rk_time_step returned the new solution.
    current_time += TIME_STEP;
    
    // Update the exact solution for comparing with the computed approximation.
    T_exact_solution.update(current_time);
    phi_exact_solution.update(current_time);
    
    // Calculate exact error.
    info("Calculating error (exact).");
    
    Hermes::vector<double> exact_errors;
    exact_errors.push_back(hermes2d.calc_rel_error(&T_current_time, &T_exact_solution, HERMES_H1_NORM));
    exact_errors.push_back(hermes2d.calc_rel_error(&phi_current_time, &phi_exact_solution, HERMES_H1_NORM));
    
    double maxerr = std::max(exact_errors[0], exact_errors[1])*100;
    info("Exact solution error for T (H1 norm): %g %%", exact_errors[0]*100);
    info("Exact solution error for phi (H1 norm): %g %%", exact_errors[1]*100);
    info("Exact solution error (maximum): %g %%", maxerr);
    
    // Prepare previous time level solution for the next time step.
    T_prev_time.copy(&T_current_time);
    phi_prev_time.copy(&phi_current_time);
  }
  while (fabs(current_time - T_FINAL) > 1e-12);
  
  cpu_time.tick();
  verbose("Total running time: %g s", cpu_time.accumulated());
  
  info("Coordinate (  0,  0) T value = %lf", T_current_time.get_pt_value(0.0, 0.0));
  info("Coordinate ( 25, 25) T value = %lf", T_current_time.get_pt_value(25.0, 25.0));
  info("Coordinate ( 75, 25) T value = %lf", T_current_time.get_pt_value(75.0, 25.0));
  info("Coordinate ( 25, 75) T value = %lf", T_current_time.get_pt_value(25.0, 75.0));
  info("Coordinate ( 75, 75) T value = %lf", T_current_time.get_pt_value(75.0, 75.0));

  info("Coordinate (  0,  0) phi value = %lf", phi_current_time.get_pt_value(0.0, 0.0));
  info("Coordinate ( 25, 25) phi value = %lf", phi_current_time.get_pt_value(25.0, 25.0));
  info("Coordinate ( 75, 25) phi value = %lf", phi_current_time.get_pt_value(75.0, 25.0));
  info("Coordinate ( 25, 75) phi value = %lf", phi_current_time.get_pt_value(25.0, 75.0));
  info("Coordinate ( 75, 75) phi value = %lf", phi_current_time.get_pt_value(75.0, 75.0));
 
  int success = 1;
  double eps = 1e-5;
  if (fabs(T_current_time.get_pt_value(0.0, 0.0) - 0.000000) > eps) {
    printf("Coordinate (  0,  0) T value = %lf\n", T_current_time.get_pt_value(0.0, 0.0));
    success = 0;
  }
  
  if (fabs(T_current_time.get_pt_value(25.0, 25.0) - 1.004828) > eps) {
    printf("Coordinate ( 25, 25) T value = %lf\n", T_current_time.get_pt_value(25.0, 25.0));
    success = 0;
  }
  
  if (fabs(T_current_time.get_pt_value(75.0, 25.0) - 1.004828) > eps) {
    printf("Coordinate ( 75, 25) T value = %lf\n", T_current_time.get_pt_value(75.0, 25.0));
    success = 0;
  }
  
  if (fabs(T_current_time.get_pt_value(25.0, 75.0) - 1.004828) > eps) {
    printf("Coordinate ( 25, 75) T value = %lf\n", T_current_time.get_pt_value(25.0, 75.0));
    success = 0;
  }
  
  if (fabs(T_current_time.get_pt_value(75.0, 75.0) - 1.004828) > eps) {
    printf("Coordinate ( 75, 75) T value = %lf\n", T_current_time.get_pt_value(75.0, 75.0));
    success = 0;
  }
  
  if (fabs(phi_current_time.get_pt_value(0.0, 0.0) - 0.000000) > eps) {
    printf("Coordinate (  0,  0) phi value = %lf\n", phi_current_time.get_pt_value(0.0, 0.0));
    success = 0;
  }
  
  if (fabs(phi_current_time.get_pt_value(25.0, 25.0) - 0.409032) > eps) {
    printf("Coordinate ( 25, 25) phi value = %lf\n", phi_current_time.get_pt_value(25.0, 25.0));
    success = 0;
  }
  
  if (fabs(phi_current_time.get_pt_value(75.0, 25.0) - 1.233987) > eps) {
    printf("Coordinate ( 75, 25) phi value = %lf\n", phi_current_time.get_pt_value(75.0, 25.0));
    success = 0;
  }
  
  if (fabs(phi_current_time.get_pt_value(25.0, 75.0) - 1.233987) > eps) {
    printf("Coordinate ( 25, 75) phi value = %lf\n", phi_current_time.get_pt_value(25.0, 75.0));
    success = 0;
  }
  
  if (fabs(phi_current_time.get_pt_value(75.0, 75.0) - 3.724930) > eps) {
    printf("Coordinate ( 75, 75) phi value = %lf\n", phi_current_time.get_pt_value(75.0, 75.0));
    success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
