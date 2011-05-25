#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/weakforms_h1.h"

/* Nonlinearity lambda(u) = pow(u, alpha) */

class CustomNonlinearity : public HermesFunction
{
public:
  CustomNonlinearity(double alpha): HermesFunction()
  {
    this->is_const = false;
    this->alpha = alpha;
  }

  virtual scalar value(double u) const
  {
    return 1 + pow(u, alpha);
  }

  virtual Ord value(Ord u) const
  {
    // If alpha is not an integer, then the function
    // is non-polynomial. 
    // NOTE: Setting Ord to 10 is safe but costly,
    // one could save here by looking at special cases 
    // of alpha. 
    return Ord(10);
  }

  virtual scalar derivative(double u) const
  {
    return alpha * pow(u, alpha - 1.0);
  }

  virtual Ord derivative(Ord u) const
  {
    // Same comment as above applies.
    return Ord(10);
  }

  protected:
    double alpha;
};

/* Initial condition for the Newton's method */

class CustomInitialCondition : public ExactSolutionScalar
{
public:
  CustomInitialCondition(Mesh* mesh) : ExactSolutionScalar(mesh) 
  {
  };

  virtual scalar value(double x, double y) const 
  {
    return (x+10) * (y+10) / 100. + 2;
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const 
  {
    dx = (y+10) / 100.;
    dy = (x+10) / 100.;
  };

  virtual Ord ord(Ord x, Ord y) const 
  {
    return x*y;
  }
};

/* Essential boundary conditions */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition {
public:
  CustomEssentialBCNonConst(std::string marker) 
           : EssentialBoundaryCondition(Hermes::vector<std::string>(marker))
  {
  }

  inline EssentialBCValueType get_value_type() const 
  { 
    return EssentialBoundaryCondition::BC_FUNCTION; 
  }

  virtual scalar value(double x, double y, double n_x, double n_y, 
                       double t_x, double t_y) const
  {
    return (x+10) * (y+10) / 100.;
  }
};


