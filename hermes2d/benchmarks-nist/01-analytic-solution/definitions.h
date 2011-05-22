#include "hermes2d.h"

/* Right-hand side */

class CustomRightHandSide : public HermesFunctionXY
{
public:
  CustomRightHandSide(double poly_deg)
    : HermesFunctionXY(), poly_deg(poly_deg) {};

  virtual double value(double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;

  double poly_deg;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double poly_deg)
            : ExactSolutionScalar(mesh), poly_deg(poly_deg) {};

  double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double poly_deg;
};

