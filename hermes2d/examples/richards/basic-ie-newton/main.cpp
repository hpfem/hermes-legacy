#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This example solves a simple version of the time-dependent
//  Richard's equation using the backward Euler method in time 
//  combined with the Newton's method in each time step. It describes
//  infiltration into an initially dry soil. The example has a exact 
//  solution that is given in terms of a Fourier series (see a paper 
//  by Tracy). The exact solution is not used here.
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

const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 5;                   // Number of initial refinements towards boundary.
const int P_INIT = 4;                             // Initial polynomial degree.
double time_step = 5e-4;                          // Time step.
const double T_FINAL = 0.4;                       // Time interval length.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Newton's method.
const double NEWTON_TOL = 1e-6;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                 // Maximum allowed number of Newton iterations.
const double DAMPING_COEFF = 1.0;

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

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
  Solution h_time_prev;
  Solution::vector_to_solution(coeff_vec, &space, &h_time_prev);

  // Initialize views.
  ScalarView view("", new WinGeom(0, 0, 600, 500));
  view.fix_scale_width(80);

  // Visualize the initial condition.
  view.show(&h_time_prev);

  // Initialize the weak formulation.
  double current_time = 0;
  CustomWeakFormRichardsIE wf(time_step, &h_time_prev);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Time stepping:
  int ts = 1;
  do 
  {
    info("---- Time step %d, time %3.5f s", ts, current_time);

    // Perform Newton's iteration.
    bool verbose = true;
    bool jacobian_changed = true;
    bool residual_as_function = false;
    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, jacobian_changed,
			       NEWTON_TOL, NEWTON_MAX_ITER, verbose, residual_as_function, 
                               DAMPING_COEFF)) error("Newton's iteration failed.");

    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solution(coeff_vec, &space, &h_time_prev);

    // Visualize the solution.
    char title[100];
    sprintf(title, "Time %3.2f s", current_time);
    view.set_title(title);
    view.show(&h_time_prev);

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  // Cleaning up.
  delete [] coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

