#include "hermes2d.h"

// Problem parameters.
const double e_0  = 8.8541878176 * 1e-12;
const double mu_0 = 1.256 * 1e-6;
const double k = 1.0;

/* Fresnel integral */

extern "C" void fresnl( double xxa, double *ssa, double *cca );

/* Exact solution */

scalar Fn(double u);

scalar Fder(double u);

scalar Fder2(double u);

scalar der_Hr(double x, double y);

scalar der_Hrr(double x, double y);

scalar der_Hrt(double x, double y);

scalar der_Ht(double x, double y);

scalar der_Htr(double x, double y);

scalar der_Htt(double x, double y);

scalar exact0(double x, double y, scalar& dx, scalar& dy);

scalar exact1(double x, double y, scalar& dx, scalar& dy);

/* Exact solution */

class CustomExactSolution : public ExactSolutionVector
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionVector(mesh) {};
  ~CustomExactSolution() {};

  virtual scalar2 value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar2& dx, scalar2& dy) const;
  
  virtual Ord ord(Ord x, Ord y) const;
};

/* Weak forms */

class CustomWeakFormScreen : public WeakForm
{
public:
  CustomWeakFormScreen();
};
