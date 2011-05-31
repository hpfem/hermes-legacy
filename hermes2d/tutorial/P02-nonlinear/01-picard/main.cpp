#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This example uses the Picard's method to solve a nonlinear problem.
//
//  PDE: Stationary heat transfer equation with nonlinear thermal
//  conductivity, -div[lambda(u) grad u] + src(x, y) = 0.
//
//  Nonlinearity: lambda(u) = 1 + pow(u, 4).
//
//  Picard's linearization: -div[lambda(u^n) grad u^{n+1}] + src(x, y) = 0.
//
//  Domain: square (-10, 10)^2.
//
//  BC: Nonconstant Dirichlet.
//
//  The following parameters can be changed:

const int P_INIT = 2;                             // Initial polynomial degree.
const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 5;                   // Number of initial refinements towards boundary.
const double PICARD_TOL = 1e-6;                   // Stopping criterion for the Picard's method.
const int PICARD_MAX_ITER = 1000;                 // Maximum allowed number of Picard's iterations.
const double INIT_COND_CONST = 3.0;               // Constant initial condition.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
double heat_src = 1.0;
double alpha = 4.0;

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("Bdy", INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  CustomEssentialBCNonConst bc_essential("Bdy");
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);

  // Initialize previous iteration solution for the Picard's method.
  Solution sln_prev_iter(&mesh, INIT_COND_CONST);

  // Initialize the weak formulation.
  CustomNonlinearity lambda(alpha);
  HermesFunction src(-heat_src);
  CustomWeakFormPicard wf(&sln_prev_iter, &lambda, &src);

  // Perform the Picard's iteration.
  bool verbose = true;
  hermes2d.solve_picard(&wf, &space, &sln_prev_iter, matrix_solver, PICARD_TOL, 
	                PICARD_MAX_ITER, verbose);

  // Visualise the solution and mesh.
  ScalarView s_view("Solution", new WinGeom(0, 0, 440, 350));
  s_view.show_mesh(false);
  s_view.show(&sln_prev_iter);
  OrderView o_view("Mesh", new WinGeom(450, 0, 420, 350));
  o_view.show(&space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

