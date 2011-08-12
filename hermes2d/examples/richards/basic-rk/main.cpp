#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  This example solves a simple version of the time-dependent
//  Richard's equation using the backward Euler method in time 
//  combined with the Newton's method in each time step.
//
//  PDE: C(h)dh/dt - div(K(h)grad(h)) - (dK/dh)*(dh/dy) = 0
//  where K(h) = K_S*exp(alpha*h)                          for h < 0,
//        K(h) = K_S                                       for h >= 0,
//        C(h) = alpha*(theta_s - theta_r)*exp(alpha*h)    for h < 0,
//        C(h) = alpha*(theta_s - theta_r)                 for h >= 0.
//
//  Domain: square (0, 100)^2.
//
//  BC: Dirichlet, given by the initial condition.
//  IC: See the function init_cond().
//
//  The following parameters can be changed:

// If this is defined, use van Genuchten's constitutive relations, otherwise use Gardner's.
//#define CONSTITUTIVE_GENUCHTEN

const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;                   // Number of initial refinements towards boundary.
const int P_INIT = 3;                             // Initial polynomial degree.
double time_step = 5e-3;                          // Time step.
const double T_FINAL = 0.4;                       // Time interval length.
const int TIME_INTEGRATION = 1;                   // 1... implicit Euler, 2... Crank-Nicolson.
const double NEWTON_TOL = 1e-6;                   // Stopping criterion for the Newton's method.
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

// For the definition of initial condition.
int y_power = 10;

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
  mloader.load("square.mesh", &mesh);

  // Convert initial condition into a Solution.
  CustomInitialCondition u_time_prev(&mesh, y_power);
  Solution u_time_new(&mesh);
  
  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  CustomDirichletCondition bc_essential(Hermes::vector<std::string>("Bdy"), y_power);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Initialize the weak formulation.
  double current_time = 0;

  CustomWeakFormRichardsRK wf;

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Initialize views.
  ScalarView view("", new WinGeom(0, 0, 600, 500));
  view.fix_scale_width(80);

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
    if (!runge_kutta.rk_time_step(current_time, time_step, &u_time_prev, 
                                  &u_time_new, jacobian_changed, verbose)) 
    {
      error("Runge-Kutta time step failed, try to decrease time step size.");
    }

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Time %3.2f s", current_time);
    view.set_title(title);
    view.show(&u_time_new);

    // Copy solution for the new time step.
    u_time_prev.copy(&u_time_new);

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;
  }

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

