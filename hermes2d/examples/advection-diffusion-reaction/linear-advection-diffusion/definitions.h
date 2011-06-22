#include "hermes2d.h"
#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

class WeakFormLinearAdvectionDiffusion : public WeakForm
{
public:
  // Problem parameters.
  double const_f;

  WeakFormLinearAdvectionDiffusion(bool stabilization_on, bool shock_capturing_on, double b1, double b2, double epsilon);

private:
  class MatrixFormVolAdvectionDiffusion : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolAdvectionDiffusion(int i, int j, double b1, double b2, double epsilon) 
          : WeakForm::MatrixFormVol(i, j), b1(b1), b2(b2), epsilon(epsilon) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const; 

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Members.
    double b1, b2, epsilon;
  };

  // Members.
  bool stabilization_on;
  bool shock_capturing_on;
};

/* Essential BC */

class EssentialBCNonConst : public EssentialBoundaryCondition 
{
public:
  EssentialBCNonConst(std::string marker) : EssentialBoundaryCondition(Hermes::vector<std::string>(marker)) {};

  ~EssentialBCNonConst() {};

  virtual EssentialBoundaryCondition::EssentialBCValueType get_value_type() const;

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;
};
