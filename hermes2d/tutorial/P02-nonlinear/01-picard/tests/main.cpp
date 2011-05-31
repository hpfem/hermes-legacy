#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This test makes sure that example 15-picard works correctly.

const int P_INIT = 2;                             // Initial polynomial degree.
const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 5;                   // Number of initial refinements towards boundary.
const double PICARD_TOL = 1e-6;                   // Stopping criterion for the Picard's method.
const int PICARD_MAX_ITER = 1000;                 // Maximum allowed number of Picard iterations.
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
  mloader.load("../square.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("Bdy", INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  CustomEssentialBCNonConst bc_essential("Bdy");
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);

  // Initialize previous iteration solution for the Picard method.
  Solution sln_prev_iter(&mesh, INIT_COND_CONST);

  // Initialize the weak formulation.
  CustomNonlinearity lambda(alpha);
  HermesFunction src(-heat_src);
  CustomWeakFormPicard wf(&sln_prev_iter, &lambda, &src);

  // Perform the Picard's iteration.
  bool verbose = true;
  hermes2d.solve_picard(&wf, &space, &sln_prev_iter, matrix_solver, PICARD_TOL,
                        PICARD_MAX_ITER, verbose);

  ndof = Space::get_num_dofs(&space);
  info("Coordinate (-10, -10) value = %lf", sln_prev_iter.get_pt_value(-10.0, -10.0));
  info("Coordinate ( -6,  -6) value = %lf", sln_prev_iter.get_pt_value(-6.0, -6.0));
  info("Coordinate ( -2,  -2) value = %lf", sln_prev_iter.get_pt_value(-2.0, -2.0));
  info("Coordinate (  2,   2) value = %lf", sln_prev_iter.get_pt_value(2.0, 2.0));
  info("Coordinate (  6,   6) value = %lf", sln_prev_iter.get_pt_value(6.0, 6.0));
  info("Coordinate ( 10,  10) value = %lf", sln_prev_iter.get_pt_value(10.0, 10.0));

  double coor_x_y[6] = {-10.0, -6.0, -2.0, 2.0, 6.0, 10.0};
  double value[6] = {0.000000, 2.311879, 2.748742, 2.920216, 3.146187, 4.000000};
  bool success = true;
  for (int i = 0; i < 6; i++)
  {
    if (abs(value[i] - sln_prev_iter.get_pt_value(coor_x_y[i], coor_x_y[i])) > 1E-6) success = false;
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

