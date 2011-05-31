#include "hermes2d.h"

/* Initial condition for the Newton's method */

class CustomInitialCondition : public ExactSolutionScalar
{
public:
  CustomInitialCondition(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual double value (double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;
};

/* Essential boundary conditions */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition 
{
public:
  CustomEssentialBCNonConst(std::string marker)
           : EssentialBoundaryCondition(Hermes::vector<std::string>(marker)) {};

  inline EssentialBCValueType get_value_type() const 
  {
    return EssentialBoundaryCondition::BC_FUNCTION;
  }

  virtual double value(double x, double y, double n_x, double n_y, 
                       double t_x, double t_y) const;
};
