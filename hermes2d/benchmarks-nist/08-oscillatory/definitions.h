#include "hermes2d.h"

/* Custom function f */

class CustomFunction: public HermesFunction
{
public:
  CustomFunction(double alpha) : HermesFunction(), alpha(alpha) {};

  virtual scalar value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

  double alpha;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double alpha)
        : ExactSolutionScalar(mesh), alpha(alpha) {};

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double alpha;
};

/* Weak forms */

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(CustomFunction* f);

private:
  class CustomMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVol(int i, int j, double alpha)
      : WeakForm::MatrixFormVol(i, j), alpha(alpha) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    double alpha;
  };

  class CustomVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol(int i, CustomFunction* f)
          : WeakForm::VectorFormVol(i), f(f) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    CustomFunction* f;
  };
};
