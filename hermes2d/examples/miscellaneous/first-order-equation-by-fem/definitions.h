#include "hermes2d.h"

/* Weak forms */

class CustomWeakForm : public WeakForm
{
  public:
    CustomWeakForm();

  private:
    class CustomJacobian : public WeakForm::MatrixFormVol
    {
      public:
        CustomJacobian() : WeakForm::MatrixFormVol(0, 0)
        {
        }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                             Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const;
    };

    class CustomResidual : public WeakForm::VectorFormVol
    {
      public:
        CustomResidual() : WeakForm::VectorFormVol(0)
        {
        }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                             Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const;
    };
};

/* Custom non-constant Dirichlet condition */

class CustomDirichletCondition : public EssentialBoundaryCondition 
{
public:
  CustomDirichletCondition(Hermes::vector<std::string> markers, double A, double B, double C);

  virtual EssentialBoundaryCondition::EssentialBCValueType get_value_type() const; 

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const; 

  protected:
    double A, B, C;
};

