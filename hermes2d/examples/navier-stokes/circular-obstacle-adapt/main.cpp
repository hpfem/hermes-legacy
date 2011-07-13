#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace RefinementSelectors;

// The time-dependent laminar incompressible Navier-Stokes equations are
// discretized in time via the implicit Euler method. The Newton's method 
// is used to solve the nonlinear problem at each time step. We show how
// to use discontinuous ($L^2$) elements for pressure and thus make the
// velocity discreetely divergence free. Comparison to approximating the
// pressure with the standard (continuous) Taylor-Hood elements is enabled.
// The Reynolds number Re = 200 which is embarrassingly low. You
// can increase it but then you will need to make the mesh finer, and the
// computation will take more time.
//
// PDE: incompressible Navier-Stokes equations in the form
// \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
// div v = 0
//
// BC: u_1 is a time-dependent constant and u_2 = 0 on Gamma_4 (inlet)
//     u_1 = u_2 = 0 on Gamma_1 (bottom), Gamma_3 (top) and Gamma_5 (obstacle)
//     "do nothing" on Gamma_2 (outlet)
//
// Geometry: Rectangular channel containing an off-axis circular obstacle. The
//           radius and position of the circle, as well as other geometry
//           parameters can be changed in the mesh file "domain.mesh".
//
// The following parameters can be changed:

const bool STOKES = false;                        // For application of Stokes flow (creeping flow).
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 3;                   // Number of initial mesh refinements towards boundary.
#define PRESSURE_IN_L2                            // If this is defined, the pressure is approximated using
                                                  // discontinuous L2 elements (making the velocity discreetely
                                                  // divergence-free, more accurate than using a continuous
                                                  // pressure approximation). Otherwise the standard continuous
                                                  // elements are used. The results are striking - check the
                                                  // tutorial for comparisons.
const int P_INIT_VEL = 2;                         // Initial polynomial degree for velocity components
const int P_INIT_PRESSURE = 1;                    // Initial polynomial degree for pressure
                                                  // Note: P_INIT_VEL should always be greater than
                                                  // P_INIT_PRESSURE because of the inf-sup condition

// Adaptivity
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is unrefined.
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
const CandList CAND_LIST = H2D_H_ANISO;           // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See the Used Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 5.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows over
                                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
const bool NEWTON = true;                         // If NEWTON == true then the Newton's iteration is performed.
                                                  // in every time step. Otherwise the convective term is linearized
                                                  // using the velocities from the previous time step.
// Problem parameters
const double RE = 200.0;                          // Reynolds number.
const double VEL_INLET = 1.0;                     // Inlet velocity (reached after STARTUP_TIME).
const double STARTUP_TIME = 1.0;                  // During this time, inlet velocity increases gradually
                                                  // from 0 to VEL_INLET, then it stays constant.
const double TAU = 0.01;                          // Time step.
const double T_FINAL = 30000.0;                   // Time interval length.
const double NEWTON_TOL = 0.05;                   // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;                   // Maximum allowed number of Newton iterations.
const double H = 5;                               // Domain height (necessary to define the parabolic
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.                                                  // velocity profile at inlet)

// Current time (defined as global since needed in weak forms)
double TIME = 0;

// Boundary markers.
const std::string BDY_BOTTOM = "b1";
const std::string BDY_RIGHT = "b2";
const std::string BDY_TOP = "b3";
const std::string BDY_LEFT = "b4";
const std::string BDY_OBSTACLE = "b5";

// Current time (used in weak forms).
double current_time = 0;

/*// Boundary condition values for x-velocity
scalar essential_bc_values_xvel(double x, double y, double time) {
  // time-dependent inlet velocity (parabolic profile)
  double val_y = VEL_INLET * y*(H-y) / (H/2.)/(H/2.); //parabolic profile with peak VEL_INLET at y = H/2
  if (time <= STARTUP_TIME) return val_y * time/STARTUP_TIME;
  else return val_y;
}
*/

/*
void mag(int n, scalar* a, scalar* dadx, scalar* dady,
                scalar* b, scalar* dbdx, scalar* dbdy,
                scalar* out, scalar* outdx, scalar* outdy)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = sqrt(sqr(a[i]) + sqr(b[i]));
    outdx[i] = (0.5 / out[i]) * (2.0 * a[i] * dadx[i] + 2.0 * b[i] * dbdx[i]);
    outdy[i] = (0.5 / out[i]) * (2.0 * a[i] * dady[i] + 2.0 * b[i] * dbdy[i]);
  }
}
*/
int main(int argc, char* argv[])
{
  Hermes2D hermes_2D;

  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_OBSTACLE, 4, false);
  mesh.refine_towards_boundary(BDY_TOP, 4, true);     // '4' is the number of levels,
  mesh.refine_towards_boundary(BDY_BOTTOM, 4, true);  // 'true' stands for anisotropic refinements.

  // Initialize boundary conditions.
  EssentialBCNonConst bc_left_vel_x(BDY_LEFT, VEL_INLET, H, STARTUP_TIME);
  DefaultEssentialBCConst bc_other_vel_x(Hermes::vector<std::string>(BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE), 0.0);
  EssentialBCs bcs_vel_x(Hermes::vector<EssentialBoundaryCondition *>(&bc_left_vel_x, &bc_other_vel_x));
  DefaultEssentialBCConst bc_vel_y(Hermes::vector<std::string>(BDY_LEFT, BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE), 0.0);
  EssentialBCs bcs_vel_y(&bc_vel_y);
  EssentialBCs bcs_pressure;

  // Spaces for velocity components and pressure.
  H1Space xvel_space(&mesh, &bcs_vel_x, P_INIT_VEL);
  H1Space yvel_space(&mesh, &bcs_vel_y, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space p_space(&mesh, &bcs_pressure, P_INIT_PRESSURE);
#else
  H1Space p_space(&mesh, &bcs_pressure, P_INIT_PRESSURE);
#endif
  Hermes::vector<Space *> spaces = Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space);

  // Calculate and report the number of degrees of freedom.
  int ndof = Space::get_num_dofs(spaces);
  info("ndof = %d.", ndof);

  // Define projection norms.
  ProjNormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  ProjNormType p_proj_norm = HERMES_L2_NORM;
#else
  ProjNormType p_proj_norm = HERMES_H1_NORM;
#endif

  // Solutions for the Newton's iteration and time stepping.
  info("Setting initial conditions.");
  Solution xvel_sln, yvel_sln, p_sln;
  Solution xvel_ref_sln, yvel_ref_sln, p_ref_sln;
  Solution xvel_prev_time, yvel_prev_time, p_prev_time;

  // Define initial conditions on the coarse mesh.
  xvel_prev_time.set_zero(&mesh);
  yvel_prev_time.set_zero(&mesh);
  p_prev_time.set_zero(&mesh);
 
  xvel_sln.copy(&xvel_prev_time);
  yvel_sln.copy(&yvel_prev_time);
  p_sln.copy(&p_prev_time);

  // Initialize weak formulation.
  WeakForm* wf;
  if (NEWTON)
    wf = new WeakFormNSNewton(STOKES, RE, TAU, &xvel_prev_time, &yvel_prev_time);
  else
    wf = new WeakFormNSSimpleLinearization(STOKES, RE, TAU, &xvel_prev_time, &yvel_prev_time);

  // Initialize the FE problem.
  DiscreteProblem dp(wf, spaces);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  VectorView vview("velocity [m/s]", new WinGeom(0, 0, 750, 240));
  ScalarView pview("pressure [Pa]", new WinGeom(0, 290, 750, 240));
  vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  pview.fix_scale_width(80);
  pview.show_mesh(true);

  // Time-stepping loop:
  char title[100];
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    current_time += TAU;
    info("---- Time step %d, time = %g:", ts, current_time);

    // Update time-dependent essential BCs.
    info("Updating time-dependent essential BC.");
    Space::update_essential_bc_values(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space), current_time);

    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      xvel_space.set_uniform_order(P_INIT_VEL);
      yvel_space.set_uniform_order(P_INIT_VEL);
      p_space.set_uniform_order(P_INIT_PRESSURE);
    }

    // Spatial adaptivity loop. Note: xvel_prev_time, yvel_prev_time and pvel_prev_time
    // must not be changed during spatial adaptivity. 
    bool done = false; int as = 1;
    double err_est;
    do {
      info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh
      // and setup reference space.
      Hermes::vector<Space *>* ref_spaces = Space::construct_refined_spaces(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space));

      // Initialize discrete problem on the reference mesh.
      DiscreteProblem dp(wf, *ref_spaces);

      // Initialize matrix solver.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Calculate initial coefficient vector for Newton on the fine mesh.
      scalar* coeff_vec = new scalar[Space::get_num_dofs(*ref_spaces)];

      if (as == 1) {
        info("Projecting coarse mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection::project_global(*ref_spaces, Hermes::vector<MeshFunction *>(&xvel_sln, &yvel_sln, &p_sln), 
                      coeff_vec, matrix_solver);
      }
      else {
        info("Projecting previous fine mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection::project_global(*ref_spaces, Hermes::vector<MeshFunction *>(&xvel_ref_sln, &yvel_ref_sln, &p_ref_sln), 
                      coeff_vec, matrix_solver);
        delete xvel_ref_sln.get_mesh();
        delete yvel_ref_sln.get_mesh();
        delete p_ref_sln.get_mesh();
      }

      // Perform Newton's iteration.
      info("Solving nonlinear problem:");
      bool verbose = true;
      bool jacobian_changed = true;
      if (!hermes_2D.solve_newton(coeff_vec, &dp, solver, matrix, rhs, jacobian_changed,
          NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

      // Update previous time level solutions.
      Solution::vector_to_solutions(coeff_vec, Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space),
                                    Hermes::vector<Solution *>(&xvel_prev_time, &yvel_prev_time, &p_prev_time));
      if (as > 1) 
      {
        // Project the fine mesh solution onto the coarse mesh.
        info("Projecting reference solution on coarse mesh.");
        OGProjection::project_global(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space), 
                      Hermes::vector<Solution *>(&xvel_ref_sln, &yvel_ref_sln, &p_ref_sln), 
                      Hermes::vector<Solution *>(&xvel_sln, &yvel_sln, &p_sln), matrix_solver, 
                      Hermes::vector<ProjNormType>(vel_proj_norm, vel_proj_norm, p_proj_norm) );
      }

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      //Adapt* adaptivity = new Adapt(*ref_spaces);
      Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space));

      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&xvel_sln, &yvel_sln, &p_sln), 
                                 Hermes::vector<Solution *>(&xvel_ref_sln, &yvel_ref_sln, &p_ref_sln)) * 100.;

      // Report results.
      info("ndof: %d, ref_ndof: %d, err_est_rel: %g%%", 
           Space::get_num_dofs(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space)), 
           Space::get_num_dofs(*ref_spaces), err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) done = true;
      else 
      {
        info("Adapting the coarse mesh.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *>(&selector, &selector, &selector), 
                                 THRESHOLD, STRATEGY, MESH_REGULARITY);

        if (Space::get_num_dofs(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space)) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
      delete adaptivity;
      delete ref_spaces;
      delete [] coeff_vec;
    }
    while (done == false);

    // Copy new time level reference solution into prev_time.
    xvel_prev_time.copy(&xvel_ref_sln);
    yvel_prev_time.copy(&yvel_ref_sln);
    p_prev_time.copy(&p_ref_sln);

    // Show the solution at the end of time step.
    sprintf(title, "Velocity, time %g", TIME);
    vview.set_title(title);
    vview.show(&xvel_prev_time, &yvel_prev_time, HERMES_EPS_LOW);
    sprintf(title, "Pressure, time %g", TIME);
    pview.set_title(title);
    pview.show(&p_prev_time);
  }

  ndof = Space::get_num_dofs(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space));
  info("ndof = %d", ndof);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
