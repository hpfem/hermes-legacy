#include "hermes2d.h"

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual void derivatives(double x, double y, scalar& dx, scalar& dy) const;
  
  virtual double value(double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const; 
};

/* Weak forms */

class CustomJacobian : public WeakForm::MatrixFormVol
{
public:
  CustomJacobian(int i, int j) : WeakForm::MatrixFormVol(i, j) {};

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                       Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const; 

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const; 
};

class CustomResidual : public WeakForm::VectorFormVol
{
public:
  CustomResidual(int i) : WeakForm::VectorFormVol(i) {};

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                       Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const; 
};

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(std::string marker_bdy_right);
};
