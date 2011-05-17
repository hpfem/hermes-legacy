#include "hermes2d.h"

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double ALPHA, double OMEGA)
            : ExactSolutionScalar(mesh), ALPHA(ALPHA), OMEGA(OMEGA) {};

  double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double ALPHA;
  double OMEGA;
  
  double get_angle(double y, double x) const 
    {
    double theta = atan2(y, x);
    if (theta < 0)
      theta += 2 * M_PI;
    return theta;
    };


};

