#include "hermes2d.h"
#include "runge_kutta.h"
#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

/* Initial condition */

class CustomInitialCondition : public ExactSolutionScalar
{
public:
  CustomInitialCondition(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual scalar value (double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;
};


/* Weak forms */

class CustomWeakFormGPRK : public WeakForm
{
public:
  CustomWeakFormGPRK(double h, double m, double g, double omega);

private:

  class CustomFormMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomFormMatrixFormVol(int i, int j, double h, double m, double g, double omega) 
          : WeakForm::MatrixFormVol(i, j), h(h), m(m), g(g), omega(omega) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                          Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

    // Members.
    double h, m, g, omega;
  };


  class CustomFormVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomFormVectorFormVol(int i, double h, double m, double g, double omega)
          : WeakForm::VectorFormVol(i), h(h), m(m), g(g), omega(omega) {};

    template<typename Real, typename Scalar>
    Scalar vector_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                            Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

    // Members.
    double h, m, g, omega;
  };
  double h, m, g, omega;
};

