#include "hermes2d.h"
#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/weakforms_h1.h"

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

// ddC/dhh (Gardner).
scalar ddCdhh(double h);

/* Custom non-constant Dirichlet condition */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition 
{
public:
  CustomEssentialBCNonConst(std::string marker, double y_power)       
        : EssentialBoundaryCondition(Hermes::vector<std::string>(marker)), y_power(y_power) {};

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, 
                       double t_x, double t_y) const;

  double y_power;
};

/* Initial condition for the Newton's method */

class CustomInitialCondition : public ExactSolutionScalar
{
public:
  CustomInitialCondition(Mesh* mesh, double y_power) : ExactSolutionScalar(mesh), y_power(y_power) {};

  virtual double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double y_power;
};

/* Weak forms */

class CustomWeakFormRichardsIE : public WeakForm
{
public:
  CustomWeakFormRichardsIE(double time_step, Solution* u_time_prev);

private:
  // This form is custom since it contains previous time-level solution.
  class CustomVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol(int i, double time_step)
          : WeakForm::VectorFormVol(i), time_step(time_step) {};

    template<typename Real, typename Scalar>
    Scalar vector_form_ie(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                          Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

    double time_step;
  };


  class CustomFormMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomFormMatrixFormVol(int i, int j, double time_step) 
          : WeakForm::MatrixFormVol(i, j), time_step(time_step) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form_ie(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                          Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();

    double time_step;

  };

  class CustomResidualFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomResidualFormVol(int i, double time_step)
          : WeakForm::VectorFormVol(i), time_step(time_step) {};

    template<typename Real, typename Scalar>
    Scalar vector_form_ie(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                          Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();

    double time_step;

  };

};

