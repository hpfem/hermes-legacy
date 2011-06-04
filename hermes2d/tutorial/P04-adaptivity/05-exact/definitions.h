#include "hermes2d.h"

class ExactSolutionCustom : public ExactSolutionScalar
{
public:
  ExactSolutionCustom(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;
};