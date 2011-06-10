#include "hermes2d.h"

class CustomWeakFormHeatRK1 : public WeakForm
{
public:
  CustomWeakFormHeatRK1(std::string bdy_air, double alpha, double lambda, double heatcap, double rho,
                        double time_step, double* current_time_ptr, double temp_init, double t_final,
                        Solution* prev_time_sln);

private:
  // This form is custom since it contains previous time-level solution.
  class CustomVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol(int i, double time_step)
      : WeakForm::VectorFormVol(i), time_step(time_step) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    double time_step;
  };

  // This form is custom since it contains time-dependent exterior temperature.
  class CustomVectorFormSurf : public WeakForm::VectorFormSurf
  {
  public:
    CustomVectorFormSurf(int i, std::string area, double alpha, double rho, double heatcap,
                         double* current_time_ptr, double temp_init, double t_final)
          : WeakForm::VectorFormSurf(i, area), alpha(alpha), rho(rho), heatcap(heatcap), current_time_ptr(current_time_ptr),
                                     temp_init(temp_init), t_final(t_final), time_step(time_step) {};

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Time-dependent exterior temperature.
    template<typename Real>
    Real temp_ext(Real t) const;

    double alpha, rho, heatcap, *current_time_ptr, temp_init, t_final, time_step;;
  };
};

