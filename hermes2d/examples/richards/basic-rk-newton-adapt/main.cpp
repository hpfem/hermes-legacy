#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This example uses adaptivity with dynamical meshes to solve
//  the Tracy problem with arbitrary Runge-Kutta methods in time. 
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
//  IC: Flat in all elements except the top layer, within this 
//      layer the solution rises linearly to match the Dirichlet condition.
//
//  NOTE: The pressure head 'h' is between -1000 and 0. For convenience, we
//        increase it by an offset H_OFFSET = 1000. In this way we can start
//        from a zero coefficient vector.
//
//  The following parameters can be changed:

// If this is defined, use van Genuchten's constitutive relations, otherwise use Gardner's.
//#define CONSTITUTIVE_GENUCHTEN

const int INIT_GLOB_REF_NUM = 1;                  // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 6;                   // Number of initial refinements towards boundary.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
double time_step = 5e-4;                          // Time step.
const double T_FINAL = 0.4;                       // Time interval length.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Adaptivity
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_METHOD = 3;                       // 1... mesh reset to basemesh and poly degrees to P_INIT.   
                                                  // 2... one ref. layer shaved off, poly degrees reset to P_INIT.
                                                  // 3... one ref. layer shaved off, poly degrees decreased by one. 
                                                  // and just one polynomial degree subtracted.
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
                                                  // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.5;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// Newton's method
const double NEWTON_TOL = 5e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.

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
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);
  mesh.copy(&basemesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("Top", INIT_REF_NUM_BDY);

  // Initialize boundary conditions.
  CustomEssentialBCNonConst bc_essential(Hermes::vector<std::string>("Bottom", "Right", "Top", "Left"));
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof_coarse = Space::get_num_dofs(&space);
  info("ndof_coarse = %d.", ndof_coarse);

  // Initial condition vector is the zero vector. This is why we
  // use the H_OFFSET. 
  scalar* coeff_vec = new scalar[ndof_coarse];
  memset(coeff_vec, 0, ndof_coarse*sizeof(double));

  // Convert initial condition into a Solution.
  Solution h_time_prev, h_time_new;
  Solution::vector_to_solution(coeff_vec, &space, &h_time_prev);
  delete [] coeff_vec;

  // Initialize the weak formulation.
  CustomWeakFormRichardsRK wf;

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Create a refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Visualize initial condition.
  char title[100];
  ScalarView view("Initial condition", new WinGeom(0, 0, 440, 350));
  OrderView ordview("Initial mesh", new WinGeom(445, 0, 440, 350));
  view.show(&h_time_prev);
  ordview.show(&space);

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
  
  // Time stepping loop.
  double current_time = 0; int ts = 1;
  do 
  {
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
                space.adjust_element_order(-1, -1, P_INIT, P_INIT);
                break;
        default: error("Wrong global derefinement method.");
      }

      ndof_coarse = Space::get_num_dofs(&space);
    }

    // Spatial adaptivity loop. Note: h_time_prev must not be changed 
    // during spatial adaptivity. 
    bool done = false; int as = 1;
    double err_est;
    do {
      info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh and setup reference space.
      Space* ref_space = Space::construct_refined_space(&space);
      int ndof_ref = Space::get_num_dofs(ref_space);

      // Initialize discrete problem on reference mesh.
      DiscreteProblem dp(&wf, ref_space);

      // Time measurement.
      cpu_time.tick();

      // Initialize Runge-Kutta time stepping.
      RungeKutta runge_kutta(&dp, &bt, matrix_solver);

      // Perform one Runge-Kutta time step according to the selected Butcher's table.
      info("Runge-Kutta time step (t = %g s, tau = %g s, stages: %d).",
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

      // Project the fine mesh solution onto the coarse mesh.
      Solution sln_coarse;
      info("Projecting fine mesh solution on coarse mesh for error estimation.");
      OGProjection::project_global(&space, &h_time_new, &sln_coarse, matrix_solver); 

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt* adaptivity = new Adapt(&space);
      double err_est_rel_total = adaptivity->calc_err_est(&sln_coarse, &h_time_new) * 100;

      // Report results.
      info("ndof_coarse: %d, ndof_ref: %d, err_est_rel: %g%%", 
           Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_est_rel_total);

      // Time measurement.
      cpu_time.tick();

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) done = true;
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
      if(!done)
        delete h_time_new.get_mesh();
    }
    while (done == false);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(current_time, Space::get_num_dofs(&space));
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(current_time, cpu_time.accumulated());
    graph_cpu.save("conv_cpu_est.dat");

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution, time %g", current_time);
    view.set_title(title);
    view.show_mesh(false);
    view.show(&h_time_new);
    sprintf(title, "Mesh, time %g", current_time);
    ordview.set_title(title);
    ordview.show(&space);

    // Copy last reference solution into h_time_prev.
    h_time_prev.copy(&h_time_new);

    // Increase current time and counter of time steps.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
