#include "hermes2d.h"
#include "runge_kutta.h"
#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

/* Nonlinearity lambda(u) = pow(u, alpha) */

class CustomNonlinearity : public HermesFunction
{
public:
  CustomNonlinearity(double alpha);

  virtual scalar value(double u) const;

  virtual Ord value(Ord u) const;

  virtual scalar derivative(double u) const;

  virtual Ord derivative(Ord u) const;

  protected:
    double alpha;
};

/* Essential boundary condition */

class EssentialBCNonConst : public EssentialBoundaryCondition 
{
public:
  EssentialBCNonConst(std::string marker);

  virtual EssentialBCValueType get_value_type() const;

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;
};

/* Initial condition */

class CustomInitialCondition : public ExactSolutionScalar
{
public:
  CustomInitialCondition(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual scalar value (double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;
};
