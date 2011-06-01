#include "hermes2d.h"

/* Right-hand side */

class CustomFunction : public HermesFunction
{
public:
  CustomFunction(double alpha, double x_loc, double y_loc)
    : HermesFunction(), alpha(alpha), x_loc(x_loc), y_loc(y_loc) {};

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

  protected:
  double alpha;
  double x_loc;
  double y_loc;
};


/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double alpha, double x_loc, double y_loc)
            : ExactSolutionScalar(mesh), alpha(alpha), x_loc(x_loc), y_loc(y_loc) {};

  double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  protected:
  double alpha;
  double x_loc;
  double y_loc;
};
