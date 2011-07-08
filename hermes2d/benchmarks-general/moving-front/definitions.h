#include "hermes2d.h"
#include "runge_kutta.h"

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double x0, double x1, double y0, double y1, 
                      double* t_ptr, double s, double c)
    : ExactSolutionScalar(mesh), x0(x0), x1(x1), y0(y0), y1(y1), t_ptr(t_ptr), s(s), c(c) {};

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double x0, x1, y0, y1, *t_ptr, s, c;
};

/* Custom function f */

class CustomFunction: public HermesFunction
{
public:
  CustomFunction(double x0, double x1, double y0, double y1,
    double s, double c)
    : HermesFunction(), x0(x0), x1(x1), y0(y0), y1(y1), s(s), c(c) {};

  virtual double value(double x, double y, double t) const;

  virtual Ord value(Ord x, Ord y) const;

  double x0, x1, y0, y1, s, c, *t_ptr;
};

class CustomVectorFormVol : public WeakForm::VectorFormVol
{
public:
  CustomVectorFormVol(int i = 0, std::string area = HERMES_ANY,
    HermesFunction* coeff = HERMES_ONE,
    GeomType gt = HERMES_PLANAR);

  CustomVectorFormVol(int i, Hermes::vector<std::string> areas,
    HermesFunction* coeff = HERMES_ONE,
    GeomType gt = HERMES_PLANAR);

  ~CustomVectorFormVol();

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
    Geom<double> *e, ExtData<scalar> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const;

  virtual WeakForm::VectorFormVol* clone();

private:
  HermesFunction* coeff;
  GeomType gt;
};

class CustomWeakFormPoisson : public WeakFormsH1::DefaultWeakFormPoisson
{
public:
  CustomWeakFormPoisson(std::string area = HERMES_ANY, 
    HermesFunction* coeff = HERMES_ONE,
    HermesFunction* f = HERMES_ONE,
    GeomType gt = HERMES_PLANAR);
};
