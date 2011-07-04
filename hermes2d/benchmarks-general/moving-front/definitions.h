#include "hermes2d.h"

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
 CustomExactSolution(Mesh* mesh, double s, double c, double* t_ptr)
   : ExactSolutionScalar(mesh), s(s), c(c), t_ptr(t_ptr) {};

  virtual double value (double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double s, c, *t_ptr;
};

/* Custom function f */

class CustomFunction: public HermesFunction
{
public:
 CustomFunction(double x0, double x1, double y0, double y1, double t1, 
                double s, double c, double* t_ptr)
   : HermesFunction(), x0(x0), x1(x1), y0(y0), y1(y1), t1(t1), s(s), c(c), t_ptr(t_ptr) {};

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

  double x0, x1, y0, y1, t1, s, c, *t_ptr;
};
