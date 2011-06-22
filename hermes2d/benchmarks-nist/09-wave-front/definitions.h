#include "hermes2d.h"

/* Custom function f */

class CustomFunction : public HermesFunction
{
public:
  CustomFunction(double alpha, double x_loc, double y_loc, double r_zero)
        : HermesFunction(), alpha(alpha), x_loc(x_loc), y_loc(y_loc), r_zero(r_zero) {};

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

protected:
    double alpha;
    double x_loc;
    double y_loc;
    double r_zero;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double alpha, double x_loc, double y_loc, double r_zero)
        : ExactSolutionScalar(mesh), alpha(alpha), x_loc(x_loc), y_loc(y_loc), r_zero(r_zero) {};

  double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;
  
protected:
    double alpha;
    double x_loc;
    double y_loc;
    double r_zero;
};

