#define HERMES_REPORT_ALL
#include "definitions.h"

//  This example solves a simple linear first-order equation using 
//  the Discontinuous Galerkin (DG) method. The exact solution is 
//  u(x, y) = x + y.
//
//  PDE:  du/dx - du/dy = 0.
//
//  Boundary conditions: Dirichlet matching exact solution u(x, y) = x + y.
//  Prescribed in the weak sense inside the weak formulation.
//
//  Domain: Unit square.
//
//  The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const int P_INIT = 1;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 3;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double BDY_A_PARAM = 1.0;
const double BDY_B_PARAM = 1.0;
const double BDY_C_PARAM = 0.0;

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize the weak formulation.
  CustomWeakForm wf;

  // Create an L2 space with default shapeset.
  L2Space space(&mesh, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initial coefficient vector for the Newton's method.  
  scalar* coeff_vec = new scalar[ndof];
  memset(coeff_vec, 0, ndof*sizeof(scalar));

  // Perform Newton's iteration.
  bool jacobian_changed = true;
  double newton_tol = 1e-8;
  int newton_max_iter = 100; 
  bool verbose = true; 
  if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, jacobian_changed, 
                             newton_tol, newton_max_iter, verbose)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution sln;
  Solution::vector_to_solution(coeff_vec, &space, &sln);

  // VTK output.
  if (VTK_VISUALIZATION) {
    // Output solution in VTK format.
    Linearizer lin;
    bool mode_3D = true;
    lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D);
    info("Solution in VTK format saved to file %s.", "sln.vtk");

    // Output mesh and element orders in VTK format.
    Orderizer ord;
    ord.save_orders_vtk(&space, "ord.vtk");
    info("Element orders in VTK format saved to file %s.", "ord.vtk");
  }

  // Visualize the solution.
  if (HERMES_VISUALIZATION) {
    ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
    view.show(&sln, HERMES_EPS_HIGH);
    View::wait();
  }

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;
  delete [] coeff_vec;

  return 0;
}

