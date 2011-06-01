#include "hermes2d.h"

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double sigma, double tau, double rho) 
        : ExactSolutionScalar(mesh), sigma(sigma), tau(tau), rho(rho) {};

  virtual double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double sigma, tau, rho;
};

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(std::string area_1, double r, std::string area_2);
};
