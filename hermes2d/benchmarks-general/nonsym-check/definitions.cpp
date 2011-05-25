#include "hermes2d.h"

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) 
  { 
  };

  virtual void derivatives(double x, double y, scalar& dx, scalar& dy) const 
  {
    dx = cos(x);
    dy = 0;
  };

  virtual double value(double x, double y) const 
  {
    return sin(x);
  }

  virtual Ord ord(Ord x, Ord y) const 
  {
    return Ord(20);
  }
};

/* Weak forms */

class CustomJacobian : public WeakForm::MatrixFormVol
{
public:
  CustomJacobian(int i, int j) : WeakForm::MatrixFormVol(i, j) 
  { 
  }

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                       Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * u->dx[i] * v->val[i];
    return result;
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    return Ord(20);
  }
};

class CustomResidual : public WeakForm::VectorFormVol
{
public:
  CustomResidual(int i) : WeakForm::VectorFormVol(i) 
  { 
  }

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                       Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    for (int i = 0; i < n; i++)
      result += wt[i] * (u_ext[0]->dx[i] - sin(e->x[i]) - cos(e->x[i])) * v->val[i];
    return result;
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    return Ord(20);
  }
};

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(std::string marker_bdy_right) : WeakForm(1) 
  {
    // Jacobian.
    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0));
    add_matrix_form(new CustomJacobian(0, 0));

    // Residual.
    add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0));
    add_vector_form(new CustomResidual(0));
    add_vector_form_surf(new WeakFormsH1::DefaultVectorFormSurf(0, marker_bdy_right, 
                                                                new HermesFunction(1.0)));
  }
};
