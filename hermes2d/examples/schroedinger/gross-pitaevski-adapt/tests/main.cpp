#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

using namespace RefinementSelectors;

// This test makes sure that example timedep-gp-adapt works correctly.

const int INIT_REF_NUM = 3;                       // Number of initial uniform refinements.
const int P_INIT = 2;                             // Initial polynomial degree.
double time_step = 0.005;                         // Time step.
const double T_FINAL = 2*time_step + 1e-4;              // Time interval length.

// Adaptivity.
const int UNREF_FREQ = 1;                         // Every UNREF_FREQ time step the mesh is unrefined.
const int UNREF_METHOD = 3;                       // 1... mesh reset to basemesh and poly degrees to P_INIT.   
                                                  // 2... one ref. layer shaved off, poly degrees reset to P_INIT.
                                                  // 3... one ref. layer shaved off, poly degrees decreased by one. 
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See the User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 5.0;                      // Stopping criterion for hp-adaptivity
                                                  // (relative error between reference and coarse solution in percent)
const double SPACE_ERR_TOL = 1.0;                 // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Temporal adaptivity.
bool ADAPTIVE_TIME_STEP_ON = true;                // This flag decides whether adaptive time stepping will be done.
                                                  // The methods for the adaptive and fixed-step versions are set
                                                  // below. An embedded method must be used with adaptive time stepping. 
const double TIME_ERR_TOL_UPPER = 1.0;            // If rel. temporal error is greater than this threshold, decrease time 
                                                  // step size and repeat time step.
const double TIME_ERR_TOL_LOWER = 0.1;            // If rel. temporal error is less than this threshold, increase time step
                                                  // but do not repeat time step (this might need further research).
const double TIME_STEP_INC_RATIO = 1.1;           // Time step increase ratio (applied when rel. temporal error is too small).
const double TIME_STEP_DEC_RATIO = 0.8;           // Time step decrease ratio (applied when rel. temporal error is too large).

// Newton's method.
const double NEWTON_TOL_COARSE = 0.01;            // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;              // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 50;                   // Maximum allowed number of Newton iterations.

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

//ButcherTableType butcher_table_type = Implicit_RK_1;
ButcherTableType butcher_table_type = Implicit_SDIRK_CASH_3_23_embedded;

// Problem parameters.
const double h = 1;                               // Planck constant 6.626068e-34.
const double m = 1;                               // Mass of boson.
const double g = 1;                               // Coupling constant.
const double omega = 1;                           // Frequency.

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Turn off adaptive time stepping if R-K method is not embedded.
  if (bt.is_embedded() == false && ADAPTIVE_TIME_STEP_ON == true) {
    warn("R-K method not embedded, turning off adaptive time stepping.");
    ADAPTIVE_TIME_STEP_ON = false;
  }

  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("../square.mesh", &basemesh);
  mesh.copy(&basemesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Convert initial condition into a Solution.
  CustomInitialCondition psi_time_prev(&mesh);

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

  // Create a refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Graph for time step history.
  SimpleGraph time_step_graph;
  if (ADAPTIVE_TIME_STEP_ON) info("Time step history will be saved to file time_step_history.dat.");
  
  // Time stepping:
  int num_time_steps = (int)(T_FINAL/time_step + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  // Time stepping loop.
  double current_time = 0.0; int ts = 1;
  do 
  {
    info("Begin time step %d.", ts);
    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
        case 1: mesh.copy(&basemesh);
                space.set_uniform_order(P_INIT);
                break;
        case 2: mesh.unrefine_all_elements();
                space.set_uniform_order(P_INIT);
                break;
        case 3: mesh.unrefine_all_elements();
                //space.adjust_element_order(-1, P_INIT);
                space.adjust_element_order(-1, -1, P_INIT, P_INIT);
                break;
        default: error("Wrong global derefinement method.");
      }

      ndof = Space::get_num_dofs(&space);
    }
    info("ndof: %d", ndof);

    // Spatial adaptivity loop. Note: psi_time_prev must not be 
    // changed during spatial adaptivity. 
    Solution ref_sln;
    Solution* time_error_fn;
    if (bt.is_embedded() == true) time_error_fn = new Solution(&mesh);
    else time_error_fn = NULL;
    bool done = false; int as = 1;
    double err_est;
    do {
      // Construct globally refined reference mesh and setup reference space.
      Space* ref_space = Space::construct_refined_space(&space);

      // Initialize discrete problem on reference mesh.
      DiscreteProblem* ref_dp = new DiscreteProblem(&wf, ref_space);
      
      RungeKutta runge_kutta(ref_dp, &bt, matrix_solver);

      // Runge-Kutta step on the fine mesh.
      info("Runge-Kutta time step on fine mesh (t = %g s, time step = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
      bool verbose = true;
      bool jacobian_changed = true;
      if (!runge_kutta.rk_time_step(current_time, time_step, &psi_time_prev, 
                                    &ref_sln, time_error_fn, jacobian_changed, verbose, 
                                    NEWTON_TOL_FINE, NEWTON_MAX_ITER)) {
    
        error("Runge-Kutta time step failed, try to decrease time step size.");
      }

      /* If ADAPTIVE_TIME_STEP_ON == true, estimate temporal error. 
         If too large or too small, then adjust it and restart the time step. */

      double rel_err_time = 0;
      if (bt.is_embedded() == true) {
        info("Calculating temporal error estimate.");

        rel_err_time = hermes2d.calc_norm(time_error_fn, HERMES_H1_NORM) / 
                       hermes2d.calc_norm(&ref_sln, HERMES_H1_NORM) * 100;
        if (ADAPTIVE_TIME_STEP_ON == false) info("rel_err_time: %g%%", rel_err_time);
      }

      if (ADAPTIVE_TIME_STEP_ON) {
        if (rel_err_time > TIME_ERR_TOL_UPPER) {
          info("rel_err_time %g%% is above upper limit %g%%", rel_err_time, TIME_ERR_TOL_UPPER);
          info("Decreasing time step from %g to %g s and restarting time step.", 
               time_step, time_step * TIME_STEP_DEC_RATIO);
          time_step *= TIME_STEP_DEC_RATIO;
          delete ref_space;
          delete ref_dp;
          continue;
        }
        else if (rel_err_time < TIME_ERR_TOL_LOWER) {
          info("rel_err_time = %g%% is below lower limit %g%%", rel_err_time, TIME_ERR_TOL_LOWER);
          info("Increasing time step from %g to %g s.", time_step, time_step * TIME_STEP_INC_RATIO);
          time_step *= TIME_STEP_INC_RATIO;
          delete ref_space;
          delete ref_dp;
          continue;
        }
        else {
          info("rel_err_time = %g%% is in acceptable interval (%g%%, %g%%)", 
            rel_err_time, TIME_ERR_TOL_LOWER, TIME_ERR_TOL_UPPER);
        }

        // Add entry to time step history graph.
        time_step_graph.add_values(current_time, time_step);
        time_step_graph.save("time_step_history.dat");
      }

      /* Estimate spatial errors and perform mesh refinement */

      info("Spatial adaptivity step %d.", as);

      // Project the fine mesh solution onto the coarse mesh.
      Solution sln;
      info("Projecting fine mesh solution on coarse mesh for error estimation.");
      OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 

      // Calculate element errors and spatial error estimate.
      info("Calculating spatial error estimate.");
      Adapt* adaptivity = new Adapt(&space);
      double err_rel_space = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

      // Report results.
      info("ndof: %d, ref_ndof: %d, err_rel_space: %g%%", 
           Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_rel_space);

      // If err_est too large, adapt the mesh.
      if (err_rel_space < SPACE_ERR_TOL) done = true;
      else 
      {
        info("Adapting the coarse mesh.");
        done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

        if (Space::get_num_dofs(&space) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }
      
      // Clean up.
      delete adaptivity;
      delete ref_space;
      delete ref_dp;
    }
    while (done == false);

    // Copy last reference solution into psi_time_prev.
    psi_time_prev.copy(&ref_sln);

    // Increase current time and counter of time steps.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  AbsFilter mag2(&psi_time_prev);
  int success = 1;
  double eps = 1e-5;
  double val = std::abs(mag2.get_pt_value(0.1, 0.1));
  info("Coordinate ( 0.1, 0.1) xvel value = %lf", val);
  if (fabs(val - (0.636525)) > eps) {
    printf("Coordinate ( 0.1, 0.1) xvel value = %lf\n", val);
    success = 0;
  }

  val = std::abs(mag2.get_pt_value(0.1, -0.1));
  info("Coordinate ( 0.1, -0.1) xvel value = %lf", val);
  if (fabs(val - (0.636525)) > eps) {
    printf("Coordinate ( 0.1, -0.1) xvel value = %lf\n", val);
    success = 0;
  }

  val = std::abs(mag2.get_pt_value(0.2, 0.1));
  info("Coordinate ( 0.2, 0.1) xvel value = %lf", val);
  if (fabs(val - (0.413034)) > eps) {
    printf("Coordinate ( 0.2, 0.1) xvel value = %lf\n", val);
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
