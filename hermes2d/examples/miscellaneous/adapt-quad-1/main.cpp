#define HERMES_REPORT_ALL
#include "definitions.h"

// This example shows how to turn on adaptive quadrature in the 
// evaluation of weak forms. It is derived from example P01-linear/03-poisson.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// Geometry: L-Shape domain (see file domain.mesh). 
//
// The following parameters can be changed:

const bool ADAPTIVE_QUADRATURE = true;            // Evaluate weak forms using adaptive quadrature.
const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;             // Set to "true" to enable VTK output.
const int P_INIT = 5;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5e3;        // Volume heat sources generated (for example) by electric current.        
const double FIXED_BDY_TEMP = 20.0;        // Fixed temperature on the boundary.

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

  // Initialize the weak formulation. Not providing the order determination form 
  // (or callback) turns on adaptive numerical quadrature. The quadrature begins 
  // with using a first-order rule in the entire element. Then the element is split 
  // uniformly in space and the quadrature order is increased by "adapt_order_increase".
  // Then the form is calculated again by employing the new quadrature in subelements. 
  // This provides a more accurate result. If relative error is less than 
  // "adapt_rel_error_tol", the computation stops, otherwise the same procedure is 
  // applied recursively to all four subelements. 
  int adapt_order_increase = 1;
  double adapt_rel_error_tol = 1e1;
  CustomWeakFormPoisson wf("Aluminum", new HermesFunction(LAMBDA_AL), "Copper", 
                           new HermesFunction(LAMBDA_CU), new HermesFunction(-VOLUME_HEAT_SRC),
                           ADAPTIVE_QUADRATURE, adapt_order_increase, adapt_rel_error_tol);
  
  if (ADAPTIVE_QUADRATURE)
    info("Adaptive quadrature ON.");    
  else
    info("Adaptive quadrature OFF.");    

  // Initialize essential boundary conditions.
  DefaultEssentialBCConst bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"), 
                                       FIXED_BDY_TEMP);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
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
  if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs)) 
    error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into a Solution.
  Solution sln;
  Solution::vector_to_solution(coeff_vec, &space, &sln);

  // VTK output.
  if (VTK_VISUALIZATION) 
  {
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
  if (HERMES_VISUALIZATION) 
  {
    ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
    // Hermes uses adaptive FEM to approximate higher-order FE solutions with linear
    // triangles for OpenGL. The second parameter of View::show() sets the error 
    // tolerance for that. Options are HERMES_EPS_LOW, HERMES_EPS_NORMAL (default), 
    // HERMES_EPS_HIGH and HERMES_EPS_VERYHIGH. The size of the graphics file grows 
    // considerably with more accurate representation, so use it wisely.
    view.show(&sln, HERMES_EPS_HIGH);
    View::wait();
  }

  // Clean up.
  delete [] coeff_vec;
  delete solver;
  delete matrix;
  delete rhs;

  return 0;
}

