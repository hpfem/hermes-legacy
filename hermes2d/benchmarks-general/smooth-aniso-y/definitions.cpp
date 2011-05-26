#include "hermes2d.h"

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) 
  { 
  };

  virtual void derivatives(double x, double y, scalar& dx, scalar& dy) const 
  {
    dx = 0;
    dy = cos(y);
  };

  virtual double value(double x, double y) const 
  {
    return sin(y);
  }

  virtual Ord ord(Ord x, Ord y) const 
  {
    return Ord(7);
  }
};

/* Custom function f */

class CustomFunction: public HermesFunction
{
public:
  CustomFunction() : HermesFunction() 
  { 
  }

  virtual scalar value(double x, double y) const 
  {
    return -sin(y);
  };

  virtual Ord value(Ord x, Ord y) const 
  {
    return Ord(7);
  }
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(std::string bdy_marker_top) : WeakForm(1) 
  {
    // Jacobian.
    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0));

    // Residual.
    add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0));
    add_vector_form(new WeakFormsH1::DefaultVectorFormVol(0, HERMES_ANY, new CustomFunction));
    add_vector_form_surf(new WeakFormsH1::DefaultVectorFormSurf(0, bdy_marker_top, new HermesFunction(1.0)));
  }
};
