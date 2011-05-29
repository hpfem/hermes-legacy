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
          : WeakForm::MatrixFormVol(i, j), epsilon(epsilon) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar val = 0;
      for (int i=0; i < n; i++) {
        val = val + wt[i] * epsilon * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
        val = val + wt[i] * (2*u->dx[i] + u->dy[i]) * v->val[i];
      }

      return val;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    double epsilon;
  };

  class CustomVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol(int i, CustomRightHandSide* rhs)
          : WeakForm::VectorFormVol(i), rhs(rhs) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar val = 0;
      for (int i=0; i < n; i++) {
        val += wt[i] * rhs->epsilon * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
        val += wt[i] * (2*u_ext[0]->dx[i] + u_ext[0]->dy[i]) * v->val[i];
		val -= wt[i] * rhs->value(e->x[i], e->y[i]) * v->val[i]; 
      }

      return val;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord val = 0;
      for (int i=0; i < n; i++) {
        val += wt[i] * rhs->epsilon * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
        val += wt[i] * (2*u_ext[0]->dx[i] + u_ext[0]->dy[i]) * v->val[i];
		val -= wt[i] * rhs->ord(e->x[i], e->y[i]) * v->val[i]; 
      }

      return val;
    }

    CustomRightHandSide* rhs;
  };
};

