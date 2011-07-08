#include "hermes2d.h"
#include "function/function.h"
#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

// Right-hand side for the 2D equation -Laplace u = f with Dirichlet BC.
class CustomRightHandSide
{
public:
  CustomRightHandSide(double K) : K(K) 
  {
  };

  template<typename Real>
  Real value(Real x, Real y);

  // Member.
  double K;
};

// Exact solution (needed in the Dirichlet condition).
class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double K) : ExactSolutionScalar(mesh), K(K) {};

  virtual scalar value(double x, double y) const;

  virtual void derivatives(double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  // Members.
  double K;
};

// Weak forms.
class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(CustomRightHandSide* rhs, std::string bdy_left_right, double K);

private:
  class CustomMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVol(int i, int j) : WeakForm::MatrixFormVol(i, j) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
  };

  class CustomVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol(int i, CustomRightHandSide* rhs) : WeakForm::VectorFormVol(i), rhs(rhs) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Members.
    CustomRightHandSide* rhs;
  };

  class CustomVectorFormSurfRight : public WeakForm::VectorFormSurf
  {
  public:
    CustomVectorFormSurfRight(int i, double K, std::string area) : WeakForm::VectorFormSurf(i, area), K(K) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Members.
    double K;
  };

  class CustomVectorFormSurfLeft : public WeakForm::VectorFormSurf
  {
  public:
    CustomVectorFormSurfLeft(int i, double K, std::string area) : WeakForm::VectorFormSurf(i, area), K(K) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Members.
    double K;
  };
};
