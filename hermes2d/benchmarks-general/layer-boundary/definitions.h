#include "hermes2d.h"

/* Exact solution */

class CustomExactFunction
{
public:
  CustomExactFunction(double K) : K(K) {};

  double uhat(double x);

  double duhat_dx(double x);

  double dduhat_dxx(double x);
 
protected:
  double K;
};

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double K);

  virtual scalar value (double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  ~CustomExactSolution();

  CustomExactFunction* cef;
};

/* Custom function */

class CustomFunction: public HermesFunction
{
public:
  CustomFunction(double coeff1);

  virtual scalar value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

  ~CustomFunction();

  CustomExactFunction* cef;
  double coeff1;
};

/* Weak forms */

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(CustomFunction* f);
};
