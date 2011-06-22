#include "hermes2d.h"

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(std::string material_1, double eps_1, std::string material_2, double eps_2, bool adapt_eval, int adapt_order_increase, double adapt_rel_error_tol);

private:
  class MatrixFormLaplace : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormLaplace(int i, int j, std::string area, double a = 1) : WeakForm::MatrixFormVol(i, j, area, HERMES_SYM), a(a) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Constant multiplicating the Laplacian.
    const double a;
  };
};
