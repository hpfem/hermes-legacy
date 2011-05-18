#include "hermes2d.h"

/* Right-hand side */

class CustomRightHandSide : public DefaultFunction
{
public:
  CustomRightHandSide(double alpha)
    : DefaultFunction(), alpha(alpha) 
  {
  };

  virtual double value(double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;

  protected:
  double alpha;
};


/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double alpha)
            : ExactSolutionScalar(mesh), alpha(alpha) 
  {
  };

  double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  protected:
  double alpha;
};
