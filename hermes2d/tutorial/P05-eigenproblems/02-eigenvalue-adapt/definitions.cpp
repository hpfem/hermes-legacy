#include "definitions.h"

/* Weak forms */

WeakFormEigenLeft::WeakFormEigenLeft() : WeakForm(1) 
{
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0));
}

WeakFormEigenRight::WeakFormEigenRight() : WeakForm(1) 
{
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 0));
}

