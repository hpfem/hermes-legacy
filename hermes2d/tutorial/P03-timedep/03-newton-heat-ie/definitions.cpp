#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

/* Nonlinearity lambda(u) = pow(u, alpha) */

class CustomNonlinearity : public HermesFunction
{
public:
  CustomNonlinearity(double alpha): HermesFunction()
  {
    this->is_const = false;
    this->alpha = alpha;
  }

  virtual scalar value(double u) const
  {
    return 1 + pow(u, alpha);
  }

  virtual Ord value(Ord u) const
  {
    // If alpha is not an integer, then the function
    // is non-polynomial. 
    // NOTE: Setting Ord to 10 is safe but costly,
    // one could save here by looking at special cases 
    // of alpha. 
    return Ord(10);
  }

  virtual scalar derivative(double u) const
  {
    return alpha * pow(u, alpha - 1.0);
  }

  virtual Ord derivative(Ord u) const
  {
    // Same comment as above applies.
    return Ord(10);
  }

  protected:
    double alpha;
};

/* Weak forms */

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(HermesFunction* lambda, HermesFunction* f, double tau, Solution* sln_prev_time) 
    : WeakForm(1) 
  {
    // Jacobian.
    add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 0, HERMES_ANY, new HermesFunction(1. / tau)));
    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, HERMES_ANY, lambda));

    // Residual.
    add_vector_form(new WeakFormsH1::DefaultResidualVol(0, HERMES_ANY, new HermesFunction(1. / tau)));
    add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, HERMES_ANY, lambda));
    add_vector_form(new WeakFormsH1::DefaultVectorFormVol(0, HERMES_ANY, f));
    CustomVectorFormVol* vector_form = new CustomVectorFormVol(0, tau);
    vector_form->ext.push_back(sln_prev_time);
    add_vector_form(vector_form);
  };

private:
  // This form (residual) is custom since it contains previous time level solution.
  class CustomVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol(int i, double tau) : WeakForm::VectorFormVol(i), tau(tau) 
    { 
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const 
    {
      scalar result = 0;
      Func<scalar>* u_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++)
        result -= wt[i] * u_prev_time->val[i] * v->val[i] / tau;
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      Ord result = 0;
      Func<Ord>* u_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++)
        result -= wt[i] * u_prev_time->val[i] * v->val[i] / tau;
      return result;
    }

    protected:
      double tau;
  };
};

/* Essential boundary condition */

class EssentialBCNonConst : public EssentialBoundaryCondition {
public:
  EssentialBCNonConst(std::string marker) : EssentialBoundaryCondition(Hermes::vector<std::string>())
  {
    markers.push_back(marker);
  }

  ~EssentialBCNonConst() {};

  inline EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition::BC_FUNCTION; }

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    return (x+10)*(y+10)/100.;
  }
};

/* Initial condition */

class CustomInitialCondition : public ExactSolutionScalar
{
public:
  CustomInitialCondition(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = (y+10)/100.;
    dy = (x+10)/100.;
  };

  virtual scalar value (double x, double y) const {
    return (x+10)*(y+10)/100.;
  }

  virtual Ord ord(Ord x, Ord y) const {
    return (x+10)*(y+10)/100.;
  }
};
