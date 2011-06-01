#include "hermes2d.h"

/* Nonlinearity lambda(u) = pow(u, alpha) */

class CustomNonlinearity : public HermesFunction
{
public:
  CustomNonlinearity(double alpha);

  virtual double value(double u) const;

  virtual Ord value(Ord u) const;

  protected:
    double alpha;
};

/* Weak forms */

// NOTE: The linear problem in the Picard's method is 
//       solved using the Newton's method.

class CustomWeakFormPicard : public WeakForm
{
public:
  CustomWeakFormPicard(Solution* prev_iter_sln, HermesFunction* lambda, HermesFunction* f);

private:
  class CustomJacobian : public WeakForm::MatrixFormVol
  {
  public:
    CustomJacobian(int i, int j, HermesFunction* lambda) : WeakForm::MatrixFormVol(i, j), lambda(lambda) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
    
    protected:
      HermesFunction* lambda;
  };

  class CustomResidual : public WeakForm::VectorFormVol
  {
  public:
    CustomResidual(int i, HermesFunction* lambda, HermesFunction* f) : WeakForm::VectorFormVol(i), lambda(lambda), f(f) {};

    virtual double value(int n, double *wt, Func<scalar> *u_ext[],
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    private:
      HermesFunction* lambda;
      HermesFunction* f;
  };
};

/* Essential boundary conditions */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition {
public:
  CustomEssentialBCNonConst(std::string marker) 
           : EssentialBoundaryCondition(Hermes::vector<std::string>(marker)) {};

  virtual EssentialBoundaryCondition::EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, 
                       double t_x, double t_y) const;
};
