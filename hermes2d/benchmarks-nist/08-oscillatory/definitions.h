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
      : WeakForm::MatrixFormVol(i, j), alpha(alpha) 
    { 
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      Scalar val = 0;
      for (int i=0; i < n; i++) 
      {
        Scalar x = e->x[i];
        Scalar y = e->y[i];
        Scalar r = sqrt(x*x + y*y);
        Scalar h = 1/(alpha + r);
        Scalar grad_u_grad_v = u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i];
        val += wt[i] * (grad_u_grad_v - pow(h, 4) * u->val[i] * v->val[i]);
      }

      return val;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    double alpha;
  };

  class CustomVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol(int i, CustomFunction* f)
          : WeakForm::VectorFormVol(i), f(f) 
    { 
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
    {
      scalar val = 0;
      for (int i=0; i < n; i++) {
        scalar x = e->x[i];
        scalar y = e->y[i];
        scalar r = sqrt(x*x + y*y);
        scalar h = 1/(f->alpha + r);
        scalar grad_u_grad_v = u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i];
        val += wt[i] * (grad_u_grad_v - pow(h, 4) * u_ext[0]->val[i] * v->val[i]);
		val -= wt[i] * f->value(e->x[i], e->y[i]) * v->val[i]; 
      }

      return val;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      Ord val = 0;
      for (int i=0; i < n; i++) 
      {
        Ord x = e->x[i];
        Ord y = e->y[i];
        Ord r = sqrt(x*x + y*y);
        Ord h = 1/(f->alpha + r);
        Ord grad_u_grad_v = u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i];
        val += wt[i] * (grad_u_grad_v - pow(h, 4) * u_ext[0]->val[i] * v->val[i]);
		val -= wt[i] * f->value(e->x[i], e->y[i]) * v->val[i]; 
      }

      return val;
    }

    CustomFunction* f;
  };
};

