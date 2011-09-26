#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace RefinementSelectors;

//  This example is a very simple flame propagation model (laminar flame,
//  zero flow velocity), and its purpose is to show how the Newton's method
//  is applied to a time-dependent two-equation system.
//
//  PDEs:
//
//  dT/dt - laplace T = omega(T,Y),
//  dY/dt - 1/Le * laplace Y = - omega(T,Y).
//
//  Domain: rectangle with cooled rods.
//
//  BC:  T = 1, Y = 0 on the inlet,
//       dT/dn = - kappa T on cooled rods,
//       dT/dn = 0, dY/dn = 0 elsewhere.
//
//  Time-stepping: a second order BDF formula.

const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int P_INIT = 2;                             // Initial polynomial degree.

// Newton's method.
const double DAMPING_COEFF = 1.0;
const double NEWTON_TOL = 1e-4;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 50;                   // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
const bool NEWTON = true;                         // If NEWTON == true then the Newton's iteration is performed.
                                                  // in every time step. Otherwise the convective term is linearized
                                                  // using the velocities from the previous time step.
const int PRECOND = 2;                            // Preconditioning by jacobian (1) (less GMRES iterations, more time to create precond)
                                                  // or by approximation of jacobian (2) (less time for precond creation, more GMRES iters).
                                                  // in case of jfnk,
                                                  // default Ifpack proconditioner in case of Newton.

// Problem constants.
const double TAU   = 0.05;                        // Time step.
const double T_FINAL = 60.0;                      // Time interval length.
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;


int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes_2D;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst left_t("Left", 1.0);
  EssentialBCs bcs_t(&left_t);

  DefaultEssentialBCConst left_c("Left", 0.0);
  EssentialBCs bcs_c(&left_c);

  // Create H1 spaces with default shapesets.
  H1Space* t_space = new H1Space(&mesh, &bcs_t, P_INIT);
  H1Space* c_space = new H1Space(&mesh, &bcs_c, P_INIT);
  int ndof = Space::get_num_dofs(Hermes::vector<Space*>(t_space, c_space));
  info("ndof = %d.", ndof);

  // Define initial conditions.
  InitialSolutionTemperature t_prev_time_1(&mesh, x1);
  InitialSolutionConcentration c_prev_time_1(&mesh, x1, Le);
  InitialSolutionTemperature t_prev_time_2(&mesh, x1);
  InitialSolutionConcentration c_prev_time_2(&mesh, x1, Le);
  Solution t_prev_newton;
  Solution c_prev_newton;

  // Filters for the reaction rate omega and its derivatives.
  CustomFilter omega(Hermes::vector<MeshFunction*>(&t_prev_time_1, &c_prev_time_1), Le, alpha, beta, kappa, x1, TAU);
  CustomFilterDt omega_dt(Hermes::vector<MeshFunction*>(&t_prev_time_1, &c_prev_time_1), Le, alpha, beta, kappa, x1, TAU);
  CustomFilterDc omega_dc(Hermes::vector<MeshFunction*>(&t_prev_time_1, &c_prev_time_1), Le, alpha, beta, kappa, x1, TAU);

  // Initialize visualization.
  ScalarView rview("Reaction rate", new WinGeom(0, 0, 800, 230));

  scalar* coeff_vec = new scalar[Space::get_num_dofs(Hermes::vector<Space*>(t_space, c_space))];
  memset(coeff_vec, 0, ndof * sizeof(double));
  Solution::vector_to_solutions(coeff_vec, Hermes::vector<Space*>(t_space, c_space), 
                                Hermes::vector<Solution *>(&t_prev_time_1, &c_prev_time_1));

  // Initialize the weak formulation.
  double current_time = 0;
  CustomWeakForm wf(Le, alpha, beta, kappa, x1, TAU, false, PRECOND, &omega, &omega_dt, 
                    &omega_dc, &t_prev_time_1, &c_prev_time_1, &t_prev_time_2, &c_prev_time_2);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, Hermes::vector<Space*>(t_space, c_space));

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Time stepping:
  int ts = 1;
  bool jacobian_changed = true;
  do 
  {
    info("---- Time step %d, time %3.5f s", ts, current_time);

    if (NEWTON)
    {
      // Perform Newton's iteration.
      info("Solving nonlinear problem:");
      bool verbose = true;
      bool jacobian_changed = true;
      if (!hermes_2D.solve_newton(coeff_vec, &dp, solver, matrix, rhs, jacobian_changed,
                                  NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

      Solution::vector_to_solutions(coeff_vec, Hermes::vector<Space*>(t_space, c_space), 
                                    Hermes::vector<Solution*>(&t_prev_newton, &c_prev_newton));
 
      // Saving solutions for the next time step.
      if(ts > 1)
      {
        t_prev_time_2.copy(&t_prev_time_1);
        c_prev_time_2.copy(&c_prev_time_1);
      }

      t_prev_time_1.copy(&t_prev_newton);
      c_prev_time_1.copy(&c_prev_newton);

    }

    // Visualization.
    rview.set_min_max_range(0.0,2.0);
    rview.show(&omega);

    // Increase current time and time step counter.
    current_time += TAU;
    ts++;
  }
  while (current_time < T_FINAL);

  // Clean up.
  delete [] coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
