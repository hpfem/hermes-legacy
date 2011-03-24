class ExactSolutionPoisson : public ExactSolutionScalar
{
public:
  ExactSolutionPoisson(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) {
    return x*x +y*y;
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) {
    dx = 2*x;
    dy = 2*y;
  };

  virtual Ord ord(Ord x, Ord y) {
    return x*x +y*y;
  }
};