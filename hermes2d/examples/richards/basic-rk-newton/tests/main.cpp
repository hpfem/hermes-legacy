#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"
#include "function/function.h"

using namespace RefinementSelectors;

// This test makes sure that example "richards/basic-rk-newton" works correctly.

const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 5;                   // Number of initial refinements towards boundary.
const int P_INIT = 2;                             // Initial polynomial degree.
double time_step = 5e-4;                          // Time step.
const double T_FINAL = 4*time_step + 1e-4;        // Time interval length.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Newton's method.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 1000;                 // Maximum allowed number of Newton iterations.
const double DAMPING_COEFF = 1.0;

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

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("Top", INIT_REF_NUM_BDY);

  // Initialize boundary conditions.
  CustomEssentialBCNonConst bc_essential(Hermes::vector<std::string>("Bottom", "Right", "Top", "Left"));
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Initial condition vector is the zero vector. This is why we
  // use the H_OFFSET. 
  scalar* coeff_vec = new scalar[ndof];
  memset(coeff_vec, 0, ndof*sizeof(double));

  // Convert initial condition into a Solution.
  Solution h_time_prev, h_time_new;
  Solution::vector_to_solution(coeff_vec, &space, &h_time_prev);
  delete [] coeff_vec;

  // Initialize the weak formulation.
  CustomWeakFormRichardsRK wf;

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Initialize Runge-Kutta time stepping.
  RungeKutta runge_kutta(&dp, &bt, matrix_solver);

  // Time stepping:
  double current_time = 0;
  int ts = 1;
  do 
  {
    info("---- Time step %d, time %3.5f s", ts, current_time);

    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, time step = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool jacobian_changed = true;
    bool verbose = true;
    double damping_coeff = 1.0;
    double max_allowed_residual_norm = 1e10;
    if (!runge_kutta.rk_time_step(current_time, time_step, &h_time_prev, 
                                  &h_time_new, jacobian_changed, verbose,
                                  NEWTON_TOL, NEWTON_MAX_ITER, damping_coeff,
                                  max_allowed_residual_norm)) 
    {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

    // Copy solution for the new time step.
    h_time_prev.copy(&h_time_new);

    // Increase current time.
    current_time += time_step;

    // Increase time step counter.
    ts++;
  }
  while (current_time < T_FINAL);

  info("Coordinate (  0,   0) value = %lf", h_time_prev.get_pt_value(0.0, 0.0));
  info("Coordinate ( 25,  25) value = %lf", h_time_prev.get_pt_value(25.0, 25.0));
  info("Coordinate ( 50,  50) value = %lf", h_time_prev.get_pt_value(50.0, 50.0));
  info("Coordinate ( 75,  75) value = %lf", h_time_prev.get_pt_value(75.0, 75.0));
  info("Coordinate (100, 100) value = %lf", h_time_prev.get_pt_value(100.0, 100.0));

  double coor_x_y[5] = {0.0, 25.0, 50.0, 75.0, 100.0};
  double value[5] = {0.000000, 0.005529, 1.694636, 98.825517, 0.000000};
  bool success = true;
  for (int i = 0; i < 5; i++)
  {
    if (fabs(value[i] - h_time_prev.get_pt_value(coor_x_y[i], coor_x_y[i])) > 1E-6) 
      success = false;
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

