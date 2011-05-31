#include "hermes2d.h"

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual void derivatives(double x, double y, scalar& dx, scalar& dy) const;

  virtual double value(double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;
};

/* Custom function f */

class CustomFunction: public HermesFunction
{
public:
  CustomFunction() : HermesFunction() {};

  virtual scalar value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(std::string bdy_marker_top);
};
