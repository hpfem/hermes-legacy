#define HERMES_REPORT_ALL
#include "hermes2d.h"
#include <stdio.h>

using namespace WeakFormsH1;

/* Weak forms */

class WeakFormS : public WeakForm
{
public:
  WeakFormS();
};

class WeakFormM : public WeakForm
{
public:
  WeakFormM();
};

/* Extras */

#define HERMES_REPORT_ALL

// Write the matrix in Matrix Market format.
void write_matrix_mm(const char* filename, Matrix* mat);

// Calculate mass norm vec^T*mat*vec.
double calc_mass_product(UMFPackMatrix* mat, double* vec, int length);

// Normalizes vector so that vec^T*mat*vec = 1. 
void normalize(UMFPackMatrix* mat, double* vec, int length);

// Multiply two vectors.
double scalar_product(double* vec1, double* vec2, int length);

// Calculate inner product u^T*mat*vec.
double calc_inner_product(UMFPackMatrix* mat, double* u, double* vec, int length);

void create_augmented_linear_system(SparseMatrix* matrix_S_ref, SparseMatrix* matrix_M_ref, 
                                    double* coeff_vec_ref, double lambda, UMFPackMatrix* new_matrix, 
                                    UMFPackVector* new_vector);

bool solve_newton_eigen(Space* ref_space, UMFPackMatrix* matrix_S_ref, UMFPackMatrix* matrix_M_ref, 
                        double* coeff_vec_ref, double &lambda, MatrixSolverType matrix_solver,
                        double newton_tol, double newton_abstol, int newton_max_iter);

bool solve_newton_eigen_ortho(Space* ref_space, UMFPackMatrix* matrix_S_ref, UMFPackMatrix* matrix_M_ref, 
                              double* coeff_vec_ref, double &lambda, MatrixSolverType matrix_solver,
                              double newton_tol, double newton_abstol, int newton_max_iter, bool use_ortho, 
                              double** coeff_space_ortho_ref, int index, int dim_space);

// This method always converges to the eigenvalue closest to the value of the argument lambda. 
// This is possible because the spectrum of the problem is shifted in such a way that the sought 
// eigenvalue comes to be very close to the origin where the method tends to converge.
bool solve_picard_eigen(Space* ref_space, UMFPackMatrix* matrix_S_ref, UMFPackMatrix* matrix_M_ref, 
                        double* coeff_vec_ref, double &lambda, MatrixSolverType matrix_solver,
                        double picard_tol, double picard_abstol, int picard_max_iter, bool use_shift);

bool solve_picard_eigen_ortho(Space* ref_space, UMFPackMatrix* matrix_S_ref, UMFPackMatrix* matrix_M_ref, 
                              double* coeff_vec_ref, double &lambda, MatrixSolverType matrix_solver,
                              double picard_tol, double picard_abstol, int picard_max_iter, bool use_ortho, bool use_shift, 
                              double** coeff_space_ortho_ref, int index, int dim_space);


