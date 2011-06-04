#include "hermes2d.h"

/* Custom function that is used in the exact solution and in right-hand side */

class CustomExactFunction1
{
public:
  CustomExactFunction1() {};

  double val(double x);
  
  double dx(double x);
  
  double ddxx(double x);
};

class CustomExactFunction2
{
public:
  CustomExactFunction2(double K) : K(K) {};

  double val(double x);
  
  double dx(double x);
  
  double ddxx(double x);

  double K;
};

/* Right-hand side */

class CustomRightHandSide1: public HermesFunction
{
public:
  CustomRightHandSide1(double K, double d_u, double sigma);

  virtual scalar value(double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;

  ~CustomRightHandSide1();

  CustomExactFunction1* cef1;
  CustomExactFunction2* cef2;
  double d_u, sigma;
};

class CustomRightHandSide2: public HermesFunction
{
public:
  CustomRightHandSide2(double K, double d_v);

  virtual scalar value(double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;

  ~CustomRightHandSide2();

  CustomExactFunction1* cef1;
  CustomExactFunction2* cef2;
  double d_v;
};

/* Exact solution */

class ExactSolutionFitzHughNagumo1 : public ExactSolutionScalar
{
public:
  ExactSolutionFitzHughNagumo1(Mesh* mesh);

  virtual scalar value (double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  ~ExactSolutionFitzHughNagumo1();

  CustomExactFunction1* cef1;
};

class ExactSolutionFitzHughNagumo2 : public ExactSolutionScalar
{
public:
  ExactSolutionFitzHughNagumo2(Mesh* mesh, double K);

  virtual scalar value (double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  ~ExactSolutionFitzHughNagumo2();

  CustomExactFunction2* cef2;
};

/* Weak forms */

class CustomResidual1 : public WeakForm::VectorFormVol
{
public:
  CustomResidual1(double d_u, double sigma, CustomRightHandSide1* g1)
        : WeakForm::VectorFormVol(0, HERMES_ANY), d_u(d_u), sigma(sigma), g1(g1) {};

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                       Geom<double> *e, ExtData<scalar> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const;

  virtual WeakForm::VectorFormVol* clone();

private:
  double d_u;
  double sigma;
  CustomRightHandSide1* g1;
};

class CustomResidual2 : public WeakForm::VectorFormVol
{
public:
  CustomResidual2(double d_v, CustomRightHandSide2* g2)
        : WeakForm::VectorFormVol(1, HERMES_ANY), d_v(d_v), g2(g2) {};

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                       Geom<double> *e, ExtData<scalar> *ext) const;
  
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const;
  
  virtual WeakForm::VectorFormVol* clone();
  
private:
  double d_v;
  CustomRightHandSide2* g2;
};

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(CustomRightHandSide1* g1, CustomRightHandSide2* g2);
};
