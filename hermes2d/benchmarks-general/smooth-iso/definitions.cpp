#include "hermes2d.h"

/*  Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) 
  { 
  }

  virtual void derivatives(double x, double y, scalar& dx, scalar& dy) const 
  {
    dx = cos(x)*sin(y);
    dy = sin(x)*cos(y);
  };

  virtual double value(double x, double y) const 
  {
    return sin(x)*sin(y);
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
    return -2*sin(x)*sin(y);
  };

  virtual Ord value(Ord x, Ord y) const 
  {
    return Ord(7);
  }
};
