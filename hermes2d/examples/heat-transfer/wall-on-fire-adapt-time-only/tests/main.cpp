#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"
#include "runge_kutta.h"

using namespace RefinementSelectors;

// This test makes sure that the example "heat-transfer/wall-on-fire-adapt-time-only" works correctly.

// The following parameters can be changed:
const int P_INIT = 3;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 3;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 3;                   // Number of initial uniform mesh refinements towards the boundary.
double time_step = 20;                            // Time step in seconds.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
const double TIME_TOL_UPPER = 1.0;                // If rel. temporal error is greater than this threshold, decrease time 
                                                  // step size and repeat time step.
const double TIME_TOL_LOWER = 0.5;                // If rel. temporal error is less than this threshold, increase time step
                                                  // but do not repeat time step (this might need further research).
const double TIME_STEP_INC_RATIO = 1.1;           // Time step increase ratio (applied when rel. temporal error is too small).
const double TIME_STEP_DEC_RATIO = 0.8;           // Time step decrease ratio (applied when rel. temporal error is too large).
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
ButcherTableType butcher_table_type = Implicit_SDIRK_CASH_3_23_embedded;

// Boundary markers.
const std::string BDY_FIRE = "Bottom";
const std::string BDY_LEFT = "Left";
const std::string BDY_RIGHT = "Right";
const std::string BDY_AIR  = "Top";

// Problem parameters.
const double TEMP_INIT = 20;         // Initial temperature.
const double TEMP_EXT_AIR = 20;      // Exterior temperature top;

const double ALPHA_FIRE = 25;        // Heat flux coefficient on the bottom edge.
const double ALPHA_AIR = 8;          // Heat flux coefficient on the top edge.
const double HEATCAP = 1020;         // Heat capacity.
const double RHO = 2200;             // Material density.
const double T_FINAL = time_step*15; // Length of time interval in seconds.

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../wall.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_RIGHT, 2);
  mesh.refine_towards_boundary(BDY_FIRE, INIT_REF_NUM_BDY);

  // Initialize essential boundary conditions (none).
  EssentialBCs bcs;

  // Initialize an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);
 
  // Previous and next time level solutions.
  Solution* sln_time_prev = new Solution(&mesh, TEMP_INIT);
  Solution* sln_time_new = new Solution(&mesh);
  Solution* time_error_fn = new Solution(&mesh, 0.0);

  // Initialize the weak formulation.
  double current_time = 0;
  CustomWeakFormHeatRK wf(BDY_FIRE, BDY_AIR, ALPHA_FIRE, ALPHA_AIR,
                          RHO, HEATCAP, TEMP_EXT_AIR, TEMP_INIT, &current_time);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Graph for time step history.
  SimpleGraph time_step_graph;
  info("Time step history will be saved to file time_step_history.dat.");

  // Initialize Runge-Kutta time stepping.
  RungeKutta runge_kutta(&dp, &bt, matrix_solver);

  // Time stepping loop:
  int ts = 1;
  do 
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, tau = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool jacobian_changed = false;
    bool verbose = true;
    if (!runge_kutta.rk_time_step(current_time, time_step, sln_time_prev, sln_time_new, 
                                  time_error_fn, jacobian_changed, verbose)) {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

    // Calculate relative time stepping error and decide whether the 
    // time step can be accepted. If not, then the time step size is 
    // reduced and the entire time step repeated. If yes, then another
    // check is run, and if the relative error is very low, time step 
    // is increased.
    double rel_err_time = hermes2d.calc_norm(time_error_fn, HERMES_H1_NORM) 
                          / hermes2d.calc_norm(sln_time_new, HERMES_H1_NORM) * 100;
    info("rel_err_time = %g%%", rel_err_time);
    if (rel_err_time > TIME_TOL_UPPER) {
      info("rel_err_time above upper limit %g%% -> decreasing time step from %g to %g and restarting time step.", 
           TIME_TOL_UPPER, time_step, time_step * TIME_STEP_DEC_RATIO);
      time_step *= TIME_STEP_DEC_RATIO;
      continue;
    }
    if (rel_err_time < TIME_TOL_LOWER) {
      info("rel_err_time = below lower limit %g%% -> increasing time step from %g to %g", 
           TIME_TOL_UPPER, time_step, time_step * TIME_STEP_INC_RATIO);
      time_step *= TIME_STEP_INC_RATIO;
    }

    // Add entry to the timestep graph.
    time_step_graph.add_values(current_time, time_step);
    time_step_graph.save("time_step_history.dat");

    // Update time.
    current_time += time_step;

    // Copy solution for the new time step.
    sln_time_prev->copy(sln_time_new);

    // Increase counter of time steps.
    ts++;
  } 
  while (current_time < T_FINAL);

  info("Coordinate (1.0, 0.0) value = %lf", sln_time_prev->get_pt_value(1.0, 0.0));
  info("Coordinate (1.5, 0.0) value = %lf", sln_time_prev->get_pt_value(1.5, 0.0));
  info("Coordinate (2.0, 0.0) value = %lf", sln_time_prev->get_pt_value(2.0, 0.0));
  info("Coordinate (2.5, 0.0) value = %lf", sln_time_prev->get_pt_value(2.5, 0.0));
  info("Coordinate (3.0, 0.0) value = %lf", sln_time_prev->get_pt_value(3.0, 0.0));
  info("Coordinate (3.5, 0.0) value = %lf", sln_time_prev->get_pt_value(3.5, 0.0));

  double coor_x[6] = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5};
  double value[6] = {28.942554, 38.108420, 48.615642, 59.122873, 68.288669, 74.772370};
  bool success = true;
  for (int i = 0; i < 6; i++)
  {
    if (fabs(value[i] - sln_time_prev->get_pt_value(coor_x[i], 0)) > 1E-6) success = false;
  }

  // Cleanup.
  delete sln_time_prev;
  delete sln_time_new;
  delete time_error_fn;

  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
