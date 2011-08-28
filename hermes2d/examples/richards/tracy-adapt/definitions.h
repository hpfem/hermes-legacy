#include "hermes2d.h"
#include "runge_kutta.h"
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
  CustomEssentialBCNonConst(Hermes::vector<std::string>(markers))       
        : EssentialBoundaryCondition(markers) {};

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, 
                       double t_x, double t_y) const;
};

/* Weak forms */

class CustomWeakFormRichardsRK : public WeakForm
{
public:
  CustomWeakFormRichardsRK();

private:

  class CustomJacobianFormVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomJacobianFormVol(int i, int j) 
          : WeakForm::MatrixFormVol(i, j) 
    {
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual WeakForm::MatrixFormVol* clone();
  };

  class CustomResidualFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomResidualFormVol(int i)
          : WeakForm::VectorFormVol(i) 
    {
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    virtual WeakForm::VectorFormVol* clone();
  };
};

