#include "hermes2d.h"

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double slope)
             : ExactSolutionScalar(mesh), slope(slope) {};

  virtual double value (double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double slope;
};

/* Custom function f */

class CustomFunction: public HermesFunction
{
public:
  CustomFunction(double slope)
        : HermesFunction(), slope(slope) {};

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

  double slope;
};

