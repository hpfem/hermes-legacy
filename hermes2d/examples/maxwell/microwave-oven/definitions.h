#include "hermes2d.h"
#include "weakform/weakform.h"
#include "integrals/hcurl.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/weakforms_hcurl.h"

/* Weak forms */

class CustomMatrixForm : public WeakForm::MatrixFormVol
{
public:
  CustomMatrixForm(int i, int j, double e_0, double mu_0, double mu_r, double kappa, double omega, double J, bool align_mesh) 
        : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), e_0(e_0), mu_0(mu_0), 
          mu_r(mu_r), kappa(kappa), omega(omega), J(J), align_mesh(align_mesh) {};

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                       Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                  Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

  // Gamma as a function of x, y.
  double gamma(int marker, double x, double y) const;

  Ord gamma(int marker, Ord x, Ord y) const;

  // Relative permittivity as a function of x, y.
  double er(int marker, double x, double y) const;

  Ord er(int marker, Ord x, Ord y) const;

  // Geometry of the load.
  bool in_load(double x, double y) const;

  private:
  double e_0, mu_0, mu_r, kappa, omega, J;
  bool align_mesh;
};

class CustomResidualForm : public WeakForm::VectorFormVol
{
public:
  CustomResidualForm(int i, double e_0, double mu_0, double mu_r, double kappa, double omega, double J, bool align_mesh) 
        : WeakForm::VectorFormVol(i), e_0(e_0), mu_0(mu_0), 
          mu_r(mu_r), kappa(kappa), omega(omega), J(J), align_mesh(align_mesh) {};

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                     Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                       ExtData<scalar> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                  ExtData<Ord> *ext) const;

  // Gamma as a function of x, y.
  double gamma(int marker, double x, double y) const;

  Ord gamma(int marker, Ord x, Ord y) const;

  // Relative permittivity as a function of x, y.
  double er(int marker, double x, double y) const;

  Ord er(int marker, Ord x, Ord y) const;

  // Geometry of the load.
  bool in_load(double x, double y) const;

  private:
  double e_0, mu_0, mu_r, kappa, omega, J;
  bool align_mesh;
};


class CustomVectorFormSurf : public WeakForm::VectorFormSurf
{
public:
  CustomVectorFormSurf(double omega, double J) 
        : WeakForm::VectorFormSurf(0), omega(omega), J(J) {};

  template<typename Scalar, typename Real>
  Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                          Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], 
                       Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                  Geom<Ord> *e, ExtData<Ord> *ext) const;

  double omega, J;
};

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(double e_0, double mu_0, double mu_r, double kappa, double omega, 
                 double J, bool align_mesh);
};

