#include "hermes2d.h"

/* Right-hand side */

class CustomExactFunction {
public:
  CustomExactFunction(double k, double alpha): k(k), alpha(alpha) {};

  double fn(double x, double y);

  double k, alpha;
};

/* Custom function f */

class CustomFunction : public HermesFunction
{
public:
  CustomFunction(double k, double alpha) : HermesFunction(), k(k), alpha(alpha)
  {
    cef = new CustomExactFunction(k, alpha);
  };

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

  CustomExactFunction* cef;

  double k, alpha;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double k, double alpha) : ExactSolutionScalar(mesh), k(k), alpha(alpha)
  {
    cef = new CustomExactFunction(k, alpha);
  };

  virtual double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  ~CustomExactSolution() { delete cef; }

  CustomExactFunction* cef;
  double k, alpha;
};

