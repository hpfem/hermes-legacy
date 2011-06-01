#include "hermes2d.h"

/* Custom function f */

class CustomFunction : public HermesFunction
{
public:
  CustomFunction(double alpha_w, double alpha_p, double x_w, double y_w, double r_0, 
	             double omega_c, double epsilon, double x_p, double y_p)
        : HermesFunction(), alpha_w(alpha_w), alpha_p(alpha_p), x_w(x_w), y_w(y_w), 
		  r_0(r_0), omega_c(omega_c), epsilon(epsilon), x_p(x_p), y_p(y_p) {};

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

protected:
  double alpha_w;
  double alpha_p;
  double x_w;
  double y_w;
  double r_0;
  double omega_c;
  double epsilon;
  double x_p;
  double y_p;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double alpha_w, double alpha_p, double x_w, double y_w, 
                      double r_0, double omega_c, double epsilon, double x_p, double y_p)
        : ExactSolutionScalar(mesh), alpha_w(alpha_w), alpha_p(alpha_p), x_w(x_w), y_w(y_w), 
          r_0(r_0), omega_c(omega_c), epsilon(epsilon), x_p(x_p), y_p(y_p) {};

  double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double get_angle(double y, double x) const;

protected:
  double alpha_w;
  double alpha_p;
  double x_w;
  double y_w;
  double r_0;
  double omega_c;
  double epsilon;
  double x_p;
  double y_p;
};
