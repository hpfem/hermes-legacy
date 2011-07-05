#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

using namespace RefinementSelectors;

// This test makes sure that example 22-newton-timedep-gp works correctly.

const int INIT_REF_NUM = 3;                       // Number of initial uniform refinements.
const int P_INIT = 4;                             // Initial polynomial degree.
double time_step = 0.005;                         // Time step.
const double T_FINAL = 4*time_step + 1e-4;        // Time interval length.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
// Choose one of the following time-integration methods, or define your own Butcher's table. The last number 
// in the name of each method is its order. The one before last, if present, is the number of stages.
// Explicit methods:
//   Explicit_RK_1, Explicit_RK_2, Explicit_RK_3, Explicit_RK_4.
// Implicit methods: 
//   Implicit_RK_1, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, Implicit_ESIRK_2_2, Implicit_SDIRK_2_2, 
//   Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Implicit_Lobatto_IIIA_3_4, 
//   Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_5_4.
// Embedded explicit methods:
//   Explicit_HEUN_EULER_2_12_embedded, Explicit_BOGACKI_SHAMPINE_4_23_embedded, Explicit_FEHLBERG_6_45_embedded,
//   Explicit_CASH_KARP_6_45_embedded, Explicit_DORMAND_PRINCE_7_45_embedded.
// Embedded implicit methods:
//   Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded, 
//   Implicit_SDIRK_BILLINGTON_3_23_embedded, Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded, 
//   Implicit_DIRK_ISMAIL_7_45_embedded. 

ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

// Problem constants
const double h = 1;                               // Planck constant 6.626068e-34.
const double m = 1;                               // Mass of boson.
const double g = 1;                               // Coupling constant.
const double omega = 1;                           // Frequency.

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Convert initial condition into a Solution.
  CustomInitialCondition psi_time_prev(&mesh);
  Solution psi_time_new(&mesh);

  // Initialize the weak formulation.
  double current_time = 0;

  CustomWeakFormGPRK wf(h, m, g, omega);
  
  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_essential("Bdy", 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);
 
  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Initialize Runge-Kutta time stepping.
  RungeKutta runge_kutta(&dp, &bt, matrix_solver);
  
  // Time stepping:
  int ts = 1;
  int nstep = (int)(T_FINAL/time_step + 0.5);
  for(int ts = 1; ts <= nstep; ts++)
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, time step = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool jacobian_changed = false;
    bool verbose = true;
    if (!runge_kutta.rk_time_step(current_time, time_step, &psi_time_prev, 
                                  &psi_time_new, jacobian_changed, verbose)) {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

    // Copy solution for the new time step.
    psi_time_prev.copy(&psi_time_new);

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;
  }

  AbsFilter mag2(&psi_time_prev);
  int success = 1;
  double eps = 1e-5;
  double val = std::abs(mag2.get_pt_value(0.1, 0.1));
  info("Coordinate ( 0.1, 0.1) psi value = %lf", val);
  if (fabs(val - (0.807026)) > eps) {
    printf("Coordinate ( 0.1, 0.1) psi value = %lf\n", val);
    success = 0;
  }

  val = std::abs(mag2.get_pt_value(0.1, -0.1));
  info("Coordinate ( 0.1, -0.1) psi value = %lf", val);
  if (fabs(val - (0.807026)) > eps) {
    printf("Coordinate ( 0.1, -0.1) psi value = %lf\n", val);
    success = 0;
  }

  val = std::abs(mag2.get_pt_value(0.2, 0.1));
  info("Coordinate ( 0.2, 0.1) psi value = %lf", val);
  if (fabs(val - (0.606902)) > eps) {
    printf("Coordinate ( 0.2, 0.1) psi value = %lf\n", val);
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
