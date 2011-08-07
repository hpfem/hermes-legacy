#include "hermes2d.h"
#include "function/function.h"
#include "runge_kutta.h"
#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

// K (Gardner).
scalar K(double h);

// dK/dh (Gardner).
scalar dKdh(double h);

// ddK/dhh (Gardner).
scalar ddKdhh(double h);

// C (Gardner).
scalar C(double h);

// dC/dh (Gardner).
scalar dCdh(double h);

/* Custom non-constant Dirichlet condition */

class CustomDirichletCondition : public EssentialBoundaryCondition 
{
public:
  CustomDirichletCondition(Hermes::vector<std::string> markers, double y_power)       
        : EssentialBoundaryCondition(markers), y_power(y_power) {};

  virtual EssentialBoundaryCondition::EssentialBCValueType get_value_type() const; 

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const; 

  protected:
    double y_power;
};

/* Initial condition */

class CustomInitialCondition : public ExactSolutionScalar
{
public:
  CustomInitialCondition(Mesh* mesh, double y_power) 
        : ExactSolutionScalar(mesh), y_power(y_power) {};

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual scalar value (double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;

  protected:
    double y_power;
};

/* Weak forms */

class CustomWeakFormRichardsRK : public WeakForm
{
public:
  CustomWeakFormRichardsRK();

private:

  class CustomFormMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomFormMatrixFormVol(int i, int j) 
          : WeakForm::MatrixFormVol(i, j) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                          Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

  };

  class CustomFormVectorFormVol1 : public WeakForm::VectorFormVol
  {
  public:
    CustomFormVectorFormVol1(int i)
          : WeakForm::VectorFormVol(i) {};

    template<typename Real, typename Scalar>
    Scalar vector_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                          Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

  };

  class CustomFormVectorFormVol2 : public WeakForm::VectorFormVol
  {
  public:
    CustomFormVectorFormVol2(int i)
          : WeakForm::VectorFormVol(i) {};

    template<typename Real, typename Scalar>
    Scalar vector_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                          Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

  };

};

