#include "hermes2d.h"

using namespace Hermes;

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(double Le, double alpha, double beta, double kappa, double x1, double tau, bool JFNK, bool PRECOND, Filter* omega, Filter* omega_dt, Filter* omega_dc, Solution* t_prev_time_1, Solution* c_prev_time_1, Solution* t_prev_time_2, Solution* c_prev_time_2);

  ~CustomWeakForm() {};

private:
  class PreconditionerForm_0 : public MatrixFormVol
  {
  public:
    PreconditionerForm_0(double tau, double Le) 
            : MatrixFormVol(0, 0, HERMES_ANY, HERMES_SYM), tau(tau), Le(Le) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    double tau, Le;
  };
  
  class PreconditionerForm_1 : public MatrixFormVol
  {
  public:
    PreconditionerForm_1(double tau, double Le) 
            : MatrixFormVol(1, 1, HERMES_ANY, HERMES_SYM), tau(tau), Le(Le) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    double tau, Le;
  };
    
  class JacobianFormVol_0_0 : public MatrixFormVol
  {
  public:
    JacobianFormVol_0_0(double tau) 
            : MatrixFormVol(0, 0, HERMES_ANY, HERMES_SYM), tau(tau) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    double tau;
  };
  
  class JacobianFormVol_0_1 : public MatrixFormVol
  {
  public:
    JacobianFormVol_0_1(double tau) 
            : MatrixFormVol(0, 1, HERMES_ANY, HERMES_SYM), tau(tau) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    double tau;
  };

  class JacobianFormVol_1_0 : public MatrixFormVol
  {
  public:
    JacobianFormVol_1_0(double tau) 
            : MatrixFormVol(1, 0, HERMES_ANY, HERMES_SYM), tau(tau) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    double tau;
  };

  class JacobianFormVol_1_1 : public MatrixFormVol
  {
  public:
    JacobianFormVol_1_1(double tau, double Le) 
            : MatrixFormVol(1, 1, HERMES_ANY, HERMES_SYM), tau(tau), Le(Le) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    double tau, Le;
  };

  class JacobianFormSurf_0_0 : public MatrixFormSurf
  {
  public:
    JacobianFormSurf_0_0(std::string bnd_marker, double kappa) 
            : MatrixFormSurf(0, 0, bnd_marker), kappa(kappa) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    double kappa;
  };

  class ResidualFormVol_0 : public VectorFormVol
  {
  public:
    ResidualFormVol_0(double tau) 
            : VectorFormVol(0), tau(tau)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double tau;
  };

  class ResidualFormVol_1 : public VectorFormVol
  {
  public:
    ResidualFormVol_1(double tau, double Le) 
            : VectorFormVol(1), tau(tau), Le(Le)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double tau, Le;
  };
  
  class ResidualFormSurf_0 : public VectorFormSurf
  {
  public:
    ResidualFormSurf_0(std::string bnd_marker, double kappa) 
            : VectorFormSurf(0, bnd_marker), kappa(kappa)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double kappa;
  };

  double Le;
  double alpha;
  double beta;
  double kappa;
  double x1;
};

class InitialSolutionTemperature : public ExactSolutionScalar
{
public:
  InitialSolutionTemperature(Mesh* mesh, double x1) : ExactSolutionScalar(mesh), x1(x1) {};

  virtual double value (double x, double y) const {
    return (x <= x1) ? 1.0 : exp(x1 - x);
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = (x <= x1) ? 0.0 : -exp(x1 - x);
    dy = 0.0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return -exp(x1 - x);
  }

  // Value.
  double x1;
};

class InitialSolutionConcentration : public ExactSolutionScalar
{
public:
  InitialSolutionConcentration(Mesh* mesh, double x1, double Le) : ExactSolutionScalar(mesh), x1(x1), Le(Le) {};

  virtual double value (double x, double y) const {
    return (x <= x1) ? 0.0 : 1.0 - exp(Le*(x1 - x));
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = (x <= x1) ? 0.0 : Le * exp(x1 - x);
    dy = 0.0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return exp(Le*(x1 - x));
  }

  // Value.
  double x1, Le;
};

class CustomFilter : public DXDYFilter
{
public:
  CustomFilter(Hermes::vector<MeshFunction*> solutions, double Le, double alpha, double beta, double kappa, double x1, double tau) : DXDYFilter(solutions), Le(Le), alpha(alpha), beta(beta), kappa(kappa), x1(x1), tau(tau)
  {
  }

private:
  virtual void filter_fn (int n, Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);

  double Le;
  double alpha;
  double beta;
  double kappa;
  double x1;
  double tau;
};

class CustomFilterDc : public DXDYFilter
{
public:
  CustomFilterDc(Hermes::vector<MeshFunction*> solutions, double Le, double alpha, double beta, double kappa, double x1, double tau) : DXDYFilter(solutions), Le(Le), alpha(alpha), beta(beta), kappa(kappa), x1(x1), tau(tau)
  {
  }

private:
  virtual void filter_fn (int n, Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);

  double Le;
  double alpha;
  double beta;
  double kappa;
  double x1;
  double tau;
};

class CustomFilterDt : public DXDYFilter
{
public:
  CustomFilterDt(Hermes::vector<MeshFunction*> solutions, double Le, double alpha, double beta, double kappa, double x1, double tau) : DXDYFilter(solutions), Le(Le), alpha(alpha), beta(beta), kappa(kappa), x1(x1), tau(tau)
  {
  }

private:
  virtual void filter_fn (int n, Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);

  double Le;
  double alpha;
  double beta;
  double kappa;
  double x1;
  double tau;
};
