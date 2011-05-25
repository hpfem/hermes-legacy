#include "hermes2d.h"

/* Exact solution */

class CustomExactFunction
{
public:
  CustomExactFunction(double K) : K(K) {};

  double uhat(double x) 
  {
    return 1. - (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
  }

  double duhat_dx(double x) 
  {
    return -K * (exp(K*x) - exp(-K*x)) / (exp(K) + exp(-K));
  }

  double dduhat_dxx(double x) 
  {
    return -K*K * (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
  }
 
  protected:
    double K;
};

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double K) : ExactSolutionScalar(mesh)
  {
    cef = new CustomExactFunction(K);
  };

  virtual scalar value (double x, double y) const 
  {
    return cef->uhat(x) * cef->uhat(y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const 
  {
    dx = cef->duhat_dx(x) * cef->uhat(y);
    dy = cef->uhat(x) * cef->duhat_dx(y);
  }

  virtual Ord ord(Ord x, Ord y) const 
  {
    return Ord(20);
  };

  ~CustomExactSolution() 
  { 
    delete cef;
  }

  CustomExactFunction* cef;
};

/* Custom function */

class CustomFunction: public HermesFunction
{
public:
  CustomFunction(double coeff1) : HermesFunction(), coeff1(coeff1) 
  {
    cef = new CustomExactFunction(coeff1);
  };

  virtual scalar value(double x, double y) const 
  {
    return -(-(cef->dduhat_dxx(x) * cef->uhat(y) + cef->uhat(x) * cef->dduhat_dxx(y))
	    + coeff1 * coeff1 * cef->uhat(x) * cef->uhat(y));
  }

  virtual Ord value(Ord x, Ord y) const 
  {
    return Ord(5);
  }

  ~CustomFunction() 
  { 
    delete cef;
  }

  CustomExactFunction* cef;
  double coeff1;
};

/* Weak forms */

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(CustomFunction* f) : WeakForm(1) 
  {
    // Jacobian.
    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0));
    add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 0, HERMES_ANY, 
                                                          new HermesFunction(f->coeff1*f->coeff1)));
    // Residual.
    add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0));
    add_vector_form(new WeakFormsH1::DefaultResidualVol(0, HERMES_ANY, 
                                                        new HermesFunction(f->coeff1*f->coeff1)));
    add_vector_form(new WeakFormsH1::DefaultVectorFormVol(0, HERMES_ANY, f));
  };
};
