#include "hermes2d.h"
#include "function/function.h"

/* Initial condition */

class InitialSolutionHeatTransfer : public ExactSolutionScalar
{
public:
  InitialSolutionHeatTransfer(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;
};

/* Essential BC */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition {
public:
  CustomEssentialBCNonConst(std::string marker) : EssentialBoundaryCondition(Hermes::vector<std::string>(marker)) {};

  ~CustomEssentialBCNonConst() {};

  virtual EssentialBCValueType get_value_type() const;

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;
};
