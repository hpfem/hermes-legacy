#include "hermes2d.h"

/* Right-hand side */

class CustomRightHandSide: public HermesFunction
{
public:
  CustomRightHandSide(double epsilon) : HermesFunction(), epsilon(epsilon) {};

  virtual double value(double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;

  double epsilon;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double epsilon)
        : ExactSolutionScalar(mesh), epsilon(epsilon) {};

  virtual double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double epsilon;
};

/* Weak forms */

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(CustomRightHandSide* rhs);

private:
  class CustomMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVol(int i, int j, double epsilon)
          : WeakForm::MatrixFormVol(i, j), epsilon(epsilon) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    double epsilon;
  };

  class CustomVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol(int i, CustomRightHandSide* rhs)
          : WeakForm::VectorFormVol(i), rhs(rhs) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    CustomRightHandSide* rhs;
  };
};

