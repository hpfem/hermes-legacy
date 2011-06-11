#include "hermes2d.h"

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(Hermes::vector<std::string> newton_boundaries, double heatcap, 
                 double rho, double tau, double lambda, double alpha, double temp_ext, 
                 Solution* sln_prev_time, bool JFNK = false);

  ~CustomWeakForm() {};

private:
  class JacobianFormVol : public WeakForm::MatrixFormVol
  {
  public:
    JacobianFormVol(int i, int j, double heatcap, double rho, double lambda, double tau) 
            : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), heatcap(heatcap), rho(rho), lambda(lambda), tau(tau) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    double heatcap, rho, lambda, tau;
  };

  class JacobianFormSurf : public WeakForm::MatrixFormSurf
  {
  public:
    JacobianFormSurf(int i, int j, Hermes::vector<std::string> newton_boundaries, double alpha, double lambda) 
            : WeakForm::MatrixFormSurf(i, j, newton_boundaries), alpha(alpha), lambda(lambda) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    double alpha, lambda;
  };

  class ResidualFormVol : public WeakForm::VectorFormVol
  {
  public:
    ResidualFormVol(int i, double heatcap, double rho, double lambda, double tau) 
            : WeakForm::VectorFormVol(i), heatcap(heatcap), rho(rho), lambda(lambda), tau(tau)  {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double heatcap, rho, lambda, tau;
  };
  
  class ResidualFormSurf : public WeakForm::VectorFormSurf
  {
  public:
    ResidualFormSurf(int i, Hermes::vector<std::string> newton_boundaries, double alpha, double lambda, double temp_ext) 
            : WeakForm::VectorFormSurf(i, newton_boundaries), alpha(alpha), lambda(lambda), temp_ext(temp_ext)  {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double alpha, lambda, temp_ext;
  };
};

