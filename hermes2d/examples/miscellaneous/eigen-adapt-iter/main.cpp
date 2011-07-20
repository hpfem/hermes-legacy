#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace RefinementSelectors;
using Teuchos::RCP;
using Teuchos::rcp;
using Hermes::EigenSolver;

//  This example shows how one can perform adaptivity to a selected eigenfunction
//  without calling the eigensolver again in each adaptivity step. The eigensolver 
//  is only called once at the beginning.
//
//  PDE: -Laplace u + V*u = lambda_k u,
//  where lambda_0, lambda_1, ... are the eigenvalues.
//
//  Domain: Square (-pi/2, pi/2)^2 ... file domain_square.mesh,
//          L-Shape domain         ... file domain_lshape.mesh.
//
//  BC:  Homogeneous Dirichlet.
//
//  The following parameters can be changed:

// Select one of the mesh files below.
//const char* mesh_file = "domain_square_quad_1_sym.mesh";     // Square domain with one single element (symmetric).
//const char* mesh_file = "domain_square_quad_2_sym.mesh";     // Square domain with four quad elements (symmetric).
//const char* mesh_file = "domain_lshape_quad_sym.mesh";       // L-Shape domain with quadrilateral mesh (symmetric). 
//const char* mesh_file = "domain_square_quad_2_nonsym.mesh";  // Square domain with four quad elements (non-symmetric).
//const char* mesh_file = "domain_square_tria_nonsym.mesh";    // Square domain with triangular mesh    (non-symmetric).
const char* mesh_file = "domain_lshape_tria_nonsym.mesh";    // L-Shape domain with triangular mesh   (non-symmetric).  

int TARGET_EIGENFUNCTION = 5;                     // Desired eigenfunction: 1 for the first, 2 for the second, etc.

int ITERATIVE_METHOD = 2;                         // 1 = Newton, 2 = Picard.

int P_INIT = 2;                                   // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial mesh refinements.
const double THRESHOLD = 0.2;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                                  // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double ERR_STOP = 1e-1;                     // Stopping criterion for adaptivity (rel. error tolerance between the
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Pysparse parameters.
const double PYSPARSE_TARGET_VALUE = 2.0;         // PySparse parameter: Eigenvalues in the vicinity of this number will be computed. 
const double PYSPARSE_TOL = 1e-10;                // PySparse parameter: Error tolerance.
const int PYSPARSE_MAX_ITER = 1000;               // PySparse parameter: Maximum number of iterations.

// Parameters for the Newton's and Picard's methods.
const double NEWTON_TOL = 1e-3;
const double NEWTON_ABSTOL = 1e-10;
const int NEWTON_MAX_ITER = 100;
const double PICARD_TOL = 1e-3;
const double PICARD_ABSTOL = 1e-10;
const int PICARD_MAX_ITER = 100;
const bool USE_SHIFT = false;

// ORTHOGONALITY:

// The orthogonality technologies are used to converge to the eigenvalues and eigenfunctions in the target 
// eigenspace. There are three possible settings:
// USE_ORTHO == false: No orthogonality - the method could converges not to the target eigenfucntions
// USE_ORTHO == true AND USE_IMPROVED_ORTHO == false: Standard orthogonality -  
//     all eigenpairs with eigenvalues smaller than the target eigenvalues are 
//     also computed in order to make sure that the method converges to the correct ones, forcing, in the 
//     iterative method, the computed eigenfunctions to be orthogonal to all eigenfunctions of smaller 
//     eigenvalues. This method is robust but quite expensive, because a lot of unwanted eigenpairs are 
//     computed
// SE_ORTHO == true AND USE_IMPROVED_ORTHO == true: Improved orthogonality - 
//     the method try to compute only the target eigenfunctions, in case that
//     unwanted eigenfunctions are computed, the method automatically discard them and keep all future
//     computed eigenfunctions orthogonal to the unwanted. The value THRESHOLD_ORTHO is used to decide 
//     wether a compute eigenfucntion is to keep or to discard. This method is also very robust and much more
//     efficient.  
const bool USE_ORTHO = true;
const bool USE_IMPROVED_ORTHO = true;
const double THRESHOLD_ORTHO = 0.5; 

// RECONSTRUCTION TECHNOLOGY:

// The reconstruction technology makes possible to follow the target eigenfunction, no matter how the continuous
// eigenspace is splitted by the discretization. This is possible by computing an approximation of all discrete
// eigenfunctions, corresponding to a basis of eigenfunctions for the continuous eigenspace that contains 
// the eigenfunction of index TARGET_EIGENFUNCTION, and thencomputing an linear interpolation of those.
// So in practise if the continuous eigenspace has dimensions 2, then the eigenfunction of index TARGET_EIGENFUNCTION
// and another eigenfunction are computed on each adapted mesh.
bool RECONSTRUCTION_ON = true;                    // Use eigenfunction reconstruction.

// Dimension of the continuous eigenspace that contains the eigenfunction of index TARGET_EIGENFUNCTION. 
// If the actual number of dimensions of the continuous eigenspace is unknown, an upperbound of it is also enough.
int DIMENSION_TARGET_EIGENSPACE = 1;  

int FIRST_INDEX_EIGENSPACE = 5;                   // Index of the first eigenfunction in the continuous eigenspace 
                                                  // that contains the eigenfunction of index TARGET_EIGENFUNCTION            

// ORTHOGONALIZATION TECHNOLOGY:

// The orthogonalization technology ensures that the iterative methods converge to the eigenfunction
// of index TARGET_EIGENFUNCTION. This is possible by computing an approximation of all discrete
// eigenfunctions of index less or equal to TARGET_EIGENFUNCTION. So in practise if TARGET_EIGENFUNCTION = 3, then
// also the eigenfunctions of indeces 1 and 2 are computed on each adapted mesh.
// The value of DIMENSION_SUBSPACE is set automatically in the code.
int DIMENSION_SUBSPACE = 0;              	  

// Main function.
int main(int argc, char* argv[])
{
  info("Desired eigenfunction to calculate: %d.", TARGET_EIGENFUNCTION);

  // Check the parameters
  if (FIRST_INDEX_EIGENSPACE > TARGET_EIGENFUNCTION) 
  {
    error("ERROR: FIRST_INDEX_EIGENSPACE should be less or equal to TARGET_EIGENFUNCTION");
  }
  
  // Disable reconstruction if no orthogonalizations are used
  if (USE_ORTHO == false)
    RECONSTRUCTION_ON = false; 

  // Set the value of DIMENSION_SUBSPACE
  if (RECONSTRUCTION_ON) 
  {
    DIMENSION_SUBSPACE = FIRST_INDEX_EIGENSPACE + DIMENSION_TARGET_EIGENSPACE - 1;
  }
  else{
    DIMENSION_SUBSPACE = TARGET_EIGENFUNCTION;
    DIMENSION_TARGET_EIGENSPACE = 1;
    FIRST_INDEX_EIGENSPACE = TARGET_EIGENFUNCTION;
  }

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load(mesh_file, &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) 
    mesh.refine_all_elements();
  
  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_essential("Bdy", 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);

  // Initialize the weak formulation for the left hand side.
  WeakFormS wf_S;
  WeakFormM wf_M;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView sview("", new WinGeom(0, 0, 440, 350));
  sview.fix_scale_width(50);
  OrderView oview("", new WinGeom(450, 0, 410, 350));

  // DOF convergence graph.
  SimpleGraph graph_dof;

  // Initialize matrices and matrix solver.
  SparseMatrix* matrix_S = create_matrix(matrix_solver);
  SparseMatrix* matrix_M = create_matrix(matrix_solver);

  // Assemble the matrices.
  DiscreteProblem dp_S(&wf_S, &space);
  dp_S.assemble(matrix_S);
  DiscreteProblem dp_M(&wf_M, &space);
  dp_M.assemble(matrix_M);

  // Initialize matrices.
  RCP<SparseMatrix> matrix_rcp_S = rcp(matrix_S);
  RCP<SparseMatrix> matrix_rcp_M = rcp(matrix_M);

  EigenSolver es(matrix_rcp_S, matrix_rcp_M);
  info("Calling Pysparse...");
  es.solve(DIMENSION_SUBSPACE, PYSPARSE_TARGET_VALUE, PYSPARSE_TOL, PYSPARSE_MAX_ITER);
  info("Pysparse finished.");
  es.print_eigenvalues();

  
  // Initialize subspace - coefficients for all computed eigenfunctions
  double* coeff_vec = new double[ndof];
  // Space for the orthogonal eigenfunctions
  double** coeff_space = new double*[DIMENSION_SUBSPACE];
  for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
  { 
    coeff_space[i] = new double[ndof];
  }
  

  // Read solution vectors from file and visualize it.
  double* eigenval =new double[DIMENSION_SUBSPACE];
  
  int neig = es.get_n_eigs();  
  for (int ieig = 0; ieig < neig; ieig++) 
  {
    // Get next eigenvalue from the file
    eigenval[ieig] = es.get_eigenvalue(ieig);         
    int n;
    es.get_eigenvector(ieig, &coeff_vec, &n);
    for (int i = 0; i < ndof; i++)
    {
      coeff_space[ieig][i] = coeff_vec[i];
    }
    // Normalize the eigenvector.
    normalize((UMFPackMatrix*)matrix_M, coeff_space[ieig], ndof);

  }

  // Retrieve desired eigenvalue.
  double lambda = eigenval[TARGET_EIGENFUNCTION-1];
  info("Eigenvalue on coarse mesh: %g", lambda);

  // Convert eigenvector into eigenfunction. After this, the 
  // eigenvector on the coarse mesh will not be needed anymore.
  Solution sln;
  Solution::vector_to_solution(coeff_space[TARGET_EIGENFUNCTION-1], &space, &sln);
  Solution* sln_space = new Solution[DIMENSION_SUBSPACE];
  for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
  {  
    Solution::vector_to_solution(coeff_space[i], &space, &sln_space[i]);
  }
  for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
  { 
    delete [] coeff_space[i];
  }
  delete [] coeff_vec;

  // Visualize the eigenfunction.
  info("Plotting initial eigenfunction on coarse mesh.");
  char title[100];
  sprintf(title, "Eigenfunction %d on initial mesh", neig);
  sview.set_title(title);
  sview.show_mesh(false);
  sview.show(&sln);
  sprintf(title, "Initial mesh");
  oview.set_title(title);
  oview.show(&space);
  View::wait(HERMES_WAIT_KEYPRESS);

  /*** Begin adaptivity ***/

  // Adaptivity loop:
  Solution ref_sln;
  Solution* ref_sln_space = new Solution[DIMENSION_SUBSPACE];
  Space* ref_space = NULL;  
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    ref_space = Space::construct_refined_space(&space);
    int ndof_ref = Space::get_num_dofs(ref_space);
    info("ndof: %d, ndof_ref: %d", ndof, ndof_ref);

    // Obtain initial approximation on new reference mesh.
    double* coeff_vec_ref = new double[ndof_ref];
    double** coeff_space_ref = new double*[DIMENSION_SUBSPACE];
    for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
    { 
      coeff_space_ref[i] = new double[ndof_ref];
    }
    if (as == 1) 
    {
      // Project the coarse mesh eigenfunction to the reference mesh.
      info("Projecting coarse mesh solution to reference mesh.");
      OGProjection::project_global(ref_space, &sln, coeff_vec_ref, matrix_solver);     
      for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
      {  
        OGProjection::project_global(ref_space, &sln_space[i], coeff_space_ref[i], matrix_solver);
      }
    }
    else {
      // Project the last reference mesh solution to the reference mesh.
      info("Projecting last reference mesh solution to new reference mesh.");
      OGProjection::project_global(ref_space, &ref_sln, coeff_vec_ref, matrix_solver);     
      for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
      {  
        OGProjection::project_global(ref_space, &ref_sln_space[i], coeff_space_ref[i], matrix_solver);
      }
    }
    Solution::vector_to_solution(coeff_vec_ref, ref_space, &ref_sln); 
    for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
    {  
      Solution::vector_to_solution(coeff_space_ref[i], ref_space, &ref_sln_space[i]); 
    }     

    // Initialize matrices and matrix solver on reference mesh.
    SparseMatrix* matrix_S_ref = create_matrix(matrix_solver);
    SparseMatrix* matrix_M_ref = create_matrix(matrix_solver);

    // Assemble matrices S and M on reference mesh.
    info("Assembling matrices S and M on reference mesh.");
    DiscreteProblem dp_S_ref(&wf_S, ref_space);
    dp_S_ref.assemble(matrix_S_ref);
    DiscreteProblem dp_M_ref(&wf_M, ref_space);
    dp_M_ref.assemble(matrix_M_ref);

    // Calculate eigenvalue corresponding to the new reference solution.
    lambda = calc_mass_product((UMFPackMatrix*)matrix_S_ref, coeff_space_ref[TARGET_EIGENFUNCTION-1], ndof_ref)
             / calc_mass_product((UMFPackMatrix*)matrix_M_ref, coeff_space_ref[TARGET_EIGENFUNCTION-1], ndof_ref);
    info("Initial guess for eigenvalue on reference mesh: %.12f", lambda);

    if (ITERATIVE_METHOD == 1) 
    {
      if (USE_ORTHO==false)
      {
        // No orthogonalization
        lambda = calc_mass_product((UMFPackMatrix*)matrix_S_ref, coeff_space_ref[TARGET_EIGENFUNCTION-1], ndof_ref)
               / calc_mass_product((UMFPackMatrix*)matrix_M_ref, coeff_space_ref[TARGET_EIGENFUNCTION-1], ndof_ref);
        if(!solve_newton_eigen(ref_space, (UMFPackMatrix*)matrix_S_ref, (UMFPackMatrix*)matrix_M_ref, 
	  		       coeff_space_ref[TARGET_EIGENFUNCTION-1], lambda, matrix_solver, NEWTON_TOL,
                               NEWTON_ABSTOL, NEWTON_MAX_ITER))
        error("Newton's method failed.");
      } 
      else if (USE_IMPROVED_ORTHO) 
      {
        // Improved orthogonality
        // Discarded eigenfunctions
        double** coeff_space_discard = new double*[DIMENSION_SUBSPACE];
        for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
        { 
      	   coeff_space_discard[i] = new double[ndof_ref];
        }
        // copy of the subspace from the previous mesh
        double** coeff_space_prev = new double*[DIMENSION_SUBSPACE];
        for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
        { 
      	   coeff_space_prev[i] = new double[ndof_ref];
        }
        for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
        { 
      	  for (int k = 0; k < ndof_ref; k++)
          {
            coeff_space_prev[i][k] = coeff_space_ref[i][k];
          }
        }
        // counter for number of eigenfunctions discarded
        int counter = 0;
        // counter for number of eigenfunction to keep
        int counter_target = 0;
        for (int i = DIMENSION_SUBSPACE-1; i >= 0; i--) 
        {  
          // create a temporary copy of the vector
          double* coeff_ref_tmp = new double[ndof_ref]; 
          info("Initial guess index: %d",DIMENSION_SUBSPACE-1-counter_target);
          for (int k = 0; k < ndof_ref; k++)
          {
            coeff_ref_tmp[k] = coeff_space_ref[DIMENSION_SUBSPACE-1-counter_target][k];
          }
          lambda = calc_mass_product((UMFPackMatrix*)matrix_S_ref, coeff_ref_tmp, ndof_ref)
             / calc_mass_product((UMFPackMatrix*)matrix_M_ref, coeff_ref_tmp, ndof_ref);
          if(!solve_newton_eigen_ortho(ref_space, (UMFPackMatrix*)matrix_S_ref, (UMFPackMatrix*)matrix_M_ref, 
	  		     coeff_ref_tmp, lambda, matrix_solver, NEWTON_TOL, NEWTON_ABSTOL, NEWTON_MAX_ITER,USE_ORTHO,
                             coeff_space_discard, counter+counter_target, DIMENSION_SUBSPACE))
          /*if(!solve_picard_eigen_ortho(ref_space, (UMFPackMatrix*)matrix_S_ref, (UMFPackMatrix*)matrix_M_ref, 
	  		     coeff_ref_tmp, lambda, matrix_solver, PICARD_TOL, PICARD_ABSTOL, PICARD_MAX_ITER,
                             USE_ORTHO, USE_SHIFT,
                             coeff_space_discard,counter+counter_target,DIMENSION_SUBSPACE))*/
          error("Newton's method failed.");
          // compute orthogonality between the computed eigenfunction and the target eigenspace from the previous mesh
          double inner = 0.0;
          for (int k = 0; k < DIMENSION_TARGET_EIGENSPACE; k++)
          {
            inner = inner + fabs(calc_inner_product((UMFPackMatrix*)matrix_M_ref, coeff_ref_tmp, 
               coeff_space_prev[DIMENSION_SUBSPACE-1-k], ndof_ref));
          }
          info("Orthogonality computed eigenfunction:%d %g",i,inner);
          if (inner < THRESHOLD_ORTHO)
          {
            info("Eigenfunction to discard");
            for (int k = 0; k < ndof_ref; k++)
            {
              coeff_space_ref[counter][k] = coeff_ref_tmp[k];
            }
            counter = counter + 1;
          }
          else{
            info("eigenfunction to keep: %d",DIMENSION_SUBSPACE-1-counter_target);
            for (int k = 0; k < ndof_ref; k++)
            {
              coeff_space_ref[DIMENSION_SUBSPACE-1-counter_target][k] = coeff_ref_tmp[k];
            }
            counter_target = counter_target + 1;
          }
          for (int k = 0; k < ndof_ref; k++)
          {
            coeff_space_discard[counter_target+counter-1][k] = coeff_ref_tmp[k];           
          }
          delete coeff_ref_tmp;
          // found all approximations of eigenfunctions in the target space
          if (counter_target==DIMENSION_TARGET_EIGENSPACE)
          {
            info("Number of computed eigenfunctions: %d",counter_target + counter );
            break;          
          }
          // computed all eigenfunctions
          else if((counter + counter_target) == DIMENSION_SUBSPACE)
          {
            info("Number of computed eigenfunctions: %d",counter_target + counter );
            break;
          }
        }
        for(int i = 0; i < DIMENSION_SUBSPACE; i++) 
        { 
          delete [] coeff_space_discard[i];
        }
        delete [] coeff_space_discard;
        for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
        { 
          delete [] coeff_space_prev[i];
        }
        delete [] coeff_space_prev;
      }
      else { // standard orthogonalization
        // Try to compute an approximation for the first eigenvalue in the target eigenspace
        for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
        {  
          lambda = calc_mass_product((UMFPackMatrix*)matrix_S_ref, coeff_space_ref[i], ndof_ref)
               / calc_mass_product((UMFPackMatrix*)matrix_M_ref, coeff_space_ref[i], ndof_ref);
          if(!solve_newton_eigen_ortho(ref_space, (UMFPackMatrix*)matrix_S_ref, (UMFPackMatrix*)matrix_M_ref, 
	  		     coeff_space_ref[i], lambda, matrix_solver, NEWTON_TOL,
                             NEWTON_ABSTOL, NEWTON_MAX_ITER,USE_ORTHO,
                             coeff_space_ref,i,DIMENSION_SUBSPACE))
          error("Newton's method failed.");
        }
      }
    }
    else if (ITERATIVE_METHOD == 2) 
    {
      if (USE_ORTHO==false)
      {
        // No orthogonalization
        lambda = calc_mass_product((UMFPackMatrix*)matrix_S_ref, coeff_space_ref[TARGET_EIGENFUNCTION-1], ndof_ref)
               / calc_mass_product((UMFPackMatrix*)matrix_M_ref, coeff_space_ref[TARGET_EIGENFUNCTION-1], ndof_ref);
        if(!solve_picard_eigen(ref_space, (UMFPackMatrix*)matrix_S_ref, (UMFPackMatrix*)matrix_M_ref, 
	  		     coeff_space_ref[TARGET_EIGENFUNCTION-1], lambda, matrix_solver, PICARD_TOL, 
                             PICARD_ABSTOL, PICARD_MAX_ITER, USE_SHIFT))
          error("Picard's method failed.");
      }
      else if (USE_IMPROVED_ORTHO) 
      {
        // Improved orthogonality
        // set of discarded eigenfunctions
        double** coeff_space_discard = new double*[DIMENSION_SUBSPACE];
        for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
        { 
      	   coeff_space_discard[i] = new double[ndof_ref];
        }
        // copy of the subspace from the previous maeh
        double** coeff_space_prev = new double*[DIMENSION_SUBSPACE];
        for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
        { 
      	   coeff_space_prev[i] = new double[ndof_ref];
        }
        for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
        { 
      	   for (int k = 0; k < ndof_ref; k++)
           {
              coeff_space_prev[i][k] = coeff_space_ref[i][k];
           }
        }
        // counter for number of eigenfunctions discarded
        int counter = 0;
        // counter for number of eigenfunction to keep
        int counter_target = 0;
        for (int i = DIMENSION_SUBSPACE-1; i >= 0; i--) 
        {  
          // create a temporary copy of the vector
          double* coeff_ref_tmp = new double[ndof_ref]; 
          info("Initial guess index: %d",DIMENSION_SUBSPACE-1-counter_target);
          for (int k = 0; k < ndof_ref; k++)
          {
            coeff_ref_tmp[k] = coeff_space_ref[DIMENSION_SUBSPACE-1-counter_target][k];
          }
          lambda = calc_mass_product((UMFPackMatrix*)matrix_S_ref, coeff_ref_tmp, ndof_ref)
             / calc_mass_product((UMFPackMatrix*)matrix_M_ref, coeff_ref_tmp, ndof_ref);
          if(!solve_picard_eigen_ortho(ref_space, (UMFPackMatrix*)matrix_S_ref, (UMFPackMatrix*)matrix_M_ref, 
	  		     coeff_ref_tmp, lambda, matrix_solver, PICARD_TOL, PICARD_ABSTOL, PICARD_MAX_ITER,
                             USE_ORTHO, USE_SHIFT,
                             coeff_space_discard,counter+counter_target,DIMENSION_SUBSPACE))
            error("Picard's method failed.");
          // compute orthogonality between the computed eigenfunction and the target eigenspace from the previous mesh
          double inner = 0.0;
          for (int k = 0; k < DIMENSION_TARGET_EIGENSPACE; k++)
          {
             inner = inner + fabs(calc_inner_product((UMFPackMatrix*)matrix_M_ref, coeff_ref_tmp, 
               coeff_space_prev[DIMENSION_SUBSPACE-1-k], ndof_ref));
          }
          info("Orthogonality computed eigenfunction:%d %g",i,inner);
          if (inner < THRESHOLD_ORTHO)
          {
            info("Eigenfunction to discard");
            for (int k = 0; k < ndof_ref; k++)
            {
              coeff_space_ref[counter][k] = coeff_ref_tmp[k];
            }
            counter = counter + 1;
          }
          else{
            info("eigenfunction to keep: %d",DIMENSION_SUBSPACE-1-counter_target);
            for (int k = 0; k < ndof_ref; k++)
            {
              coeff_space_ref[DIMENSION_SUBSPACE-1-counter_target][k] = coeff_ref_tmp[k];
            }
            counter_target = counter_target + 1;
          }
          for (int k = 0; k < ndof_ref; k++)
          {
            coeff_space_discard[counter_target+counter-1][k] = coeff_ref_tmp[k];  
          }
          delete coeff_ref_tmp;
          // found all approximations of eigenfunctions in the target space
          if (counter_target==DIMENSION_TARGET_EIGENSPACE)
          {
            info("Number of computed eigenfunctions: %d",counter_target + counter );
            break;          
          }
          // computed all eigenfunctions
          else if ((counter + counter_target) == DIMENSION_SUBSPACE)
          {
            info("Number of computed eigenfunctions: %d",counter_target + counter );
            break;
          }
        }
        for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
        { 
          delete [] coeff_space_discard[i];
        }
        delete [] coeff_space_discard;
        for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
        { 
          delete [] coeff_space_prev[i];
        }
        delete [] coeff_space_prev;
      }
      else { // standard orthogonalization
        // Try to compute an approximation for the first eigenvalue in the target eigenspace
        for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
        {  
          lambda = calc_mass_product((UMFPackMatrix*)matrix_S_ref, coeff_space_ref[i], ndof_ref)
               / calc_mass_product((UMFPackMatrix*)matrix_M_ref, coeff_space_ref[i], ndof_ref);
          if(!solve_picard_eigen_ortho(ref_space, (UMFPackMatrix*)matrix_S_ref, (UMFPackMatrix*)matrix_M_ref, 
	  		     coeff_space_ref[i], lambda, matrix_solver, PICARD_TOL,PICARD_ABSTOL,
                             PICARD_MAX_ITER,USE_ORTHO, USE_SHIFT,
                             coeff_space_ref,i,DIMENSION_SUBSPACE))
            error("Picard's method failed.");
        }
      }
    }
    else {
         error("Solver unknown.");
    }

    for (int i = 0; i < DIMENSION_SUBSPACE; i++)
      Solution::vector_to_solution(coeff_space_ref[i], ref_space, &ref_sln_space[i]);

    // Perform eigenfunction reconstruction.
    if (RECONSTRUCTION_ON == false)
      Solution::vector_to_solution(coeff_space_ref[TARGET_EIGENFUNCTION-1], ref_space, &ref_sln);

    else {
      double* inners = new double[DIMENSION_TARGET_EIGENSPACE];
      double* coeff_vec_rec = new double[ndof_ref];
      for (int i = 0; i < DIMENSION_TARGET_EIGENSPACE; i++)
      {
        inners[i] = calc_inner_product((UMFPackMatrix*)matrix_M_ref, coeff_space_ref[FIRST_INDEX_EIGENSPACE-1+i], coeff_vec_ref, ndof_ref);
        info("inners:%d %g",i,inners[i]);
      }
      
      for (int j = 0; j < ndof_ref; j++) 
      {
        coeff_vec_rec[j] = 0.0;
        for (int i = 0; i < DIMENSION_TARGET_EIGENSPACE; i++)
          coeff_vec_rec[j] += inners[i] * coeff_space_ref[FIRST_INDEX_EIGENSPACE-1+i][j];
          
      }
      normalize((UMFPackMatrix*)matrix_M_ref, coeff_vec_rec, ndof_ref);
      lambda = calc_mass_product((UMFPackMatrix*)matrix_S_ref, coeff_vec_rec, ndof_ref)
             / calc_mass_product((UMFPackMatrix*)matrix_M_ref, coeff_vec_rec, ndof_ref);
      info("Reconstructed Eigenvalue: %.12f",lambda);
      Solution::vector_to_solution(coeff_vec_rec, ref_space, &ref_sln);

      delete [] coeff_vec_rec;
      delete [] inners;
    }

    // Clean up.
    delete matrix_S_ref;
    delete matrix_M_ref;
    delete [] coeff_vec_ref;
    for (int i = 0; i < DIMENSION_SUBSPACE; i++) 
    { 
      delete [] coeff_space_ref[i];
    }
    delete [] coeff_space_ref;


    // Project reference solution to coarse mesh for error estimation.
    if (as > 1) 
    {
      // Project reference solution to coarse mesh.
      info("Projecting reference solution to coarse mesh for error calculation.");
      OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 
    }

    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    Adapt* adaptivity = new Adapt(&space);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
      Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_est_rel);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    ndof = Space::get_num_dofs(&space);
    if (ndof >= NDOF_STOP) done = true;

    // Clean up.
    delete adaptivity;

    //delete ref_space->get_mesh();
    delete ref_space;

    // Visualize the projection.
    info("Plotting projection of reference solution to new coarse mesh.");
    char title[100];
    sprintf(title, "Coarse mesh projection");
    sview.set_title(title);
    sview.show_mesh(false);
    sview.show(&sln);
    sprintf(title, "Coarse mesh, step %d", as);
    oview.set_title(title);
    oview.show(&space);

    // Increase the counter of performed adaptivity steps.
    if (done == false) as++;

    // Wait for keypress.
    View::wait(HERMES_WAIT_KEYPRESS);
  }
  while (done == false);

  // Wait for all views to be closed.
  info("Computation finished.");
  View::wait();

  return 0; 
};

