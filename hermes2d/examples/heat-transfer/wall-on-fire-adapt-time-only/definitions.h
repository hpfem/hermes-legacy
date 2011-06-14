#include "hermes2d.h"

/* Space-dependent thermal conductivity */

double lambda(double x, double y);

/* Time-dependent fire temperature */

template<typename Real>
Real T_fire_x(Real x);

template<typename Real>
Real T_fire_t(Real t);

/* Weak forms */

class CustomWeakFormHeatRK : public WeakForm
{
public:
  CustomWeakFormHeatRK(std::string bdy_fire, std::string bdy_air,
                       double alpha_fire, double alpha_air, double rho, double heatcap,
                       double temp_ext_air, double temp_init, double* current_time_ptr);

private:
  // This form is custom since it contains space-dependent thermal conductivity.
  class CustomJacobianVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomJacobianVol(int i, int j, double rho, double heatcap)
          : WeakForm::MatrixFormVol(i, j), rho(rho), heatcap(heatcap) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e,
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    // This is needed for the rk_time_step() method.
    virtual WeakForm::MatrixFormVol* clone();

    double rho, heatcap;
  };

  // This form is custom since it contains space-dependent thermal conductivity.
  class CustomFormResidualVol : public WeakForm::VectorFormVol
  {
  public:
    CustomFormResidualVol(int i, double rho, double heatcap)
          : WeakForm::VectorFormVol(i), rho(rho), heatcap(heatcap) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Needed for the rk_time_step() method.
    virtual WeakForm::VectorFormVol* clone(); 

    double rho, heatcap;
  };

  // Custom due to time-dependent exterior temperature.
  class CustomFormResidualSurfFire : public WeakForm::VectorFormSurf
  {
  public:
    CustomFormResidualSurfFire(int i, std::string area, double alpha_fire, double rho,
                               double heatcap, double* current_time_ptr)
          : WeakForm::VectorFormSurf(i, area), alpha_fire(alpha_fire), rho(rho),
                           heatcap(heatcap), current_time_ptr(current_time_ptr) {};

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                            Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Needed for the rk_time_step() method.
    virtual WeakForm::VectorFormSurf* clone();

    // Fire temperature as function of x and t.
    template<typename Real>
    Real T_fire(Real x, Real t) const;

    double alpha_fire, rho, heatcap, *current_time_ptr;
  };
};

