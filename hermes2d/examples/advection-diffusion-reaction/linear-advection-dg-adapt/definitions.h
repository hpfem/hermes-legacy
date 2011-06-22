#include "hermes2d.h"
#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(std::string left_bottom_bnd_part);

private:
  class MatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVol(int i, int j) : WeakForm::MatrixFormVol(i, j) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
  };

  class VectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVol(int i) : WeakForm::VectorFormVol(i) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    template<typename Real>
    Real F(Real x, Real y) const;
  };

  class MatrixFormSurface : public WeakForm::MatrixFormSurf
  {
  public:
    MatrixFormSurface(int i, int j) : WeakForm::MatrixFormSurf(i, j, H2D_DG_BOUNDARY_EDGE) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
  };

  class MatrixFormInterface : public WeakForm::MatrixFormSurf
  {
  public:
    MatrixFormInterface(int i, int j) : WeakForm::MatrixFormSurf(i, j, H2D_DG_INNER_EDGE) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
  };

  class VectorFormSurface : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormSurface(int i, std::string left_bottom_bnd_part) : WeakForm::VectorFormSurf(i, H2D_DG_BOUNDARY_EDGE), left_bottom_bnd_part(left_bottom_bnd_part) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    template<typename Real>
    Real F(Real x, Real y) const;

    template<typename Real, typename Scalar>
    Scalar g(std::string ess_bdy_marker, Real x, Real y) const;
    
    // Member.
    std::string left_bottom_bnd_part;
  };
  
  double calculate_a_dot_v(double x, double y, double vx, double vy) const;

  Ord calculate_a_dot_v(Ord x, Ord y, Ord vx, Ord vy) const;

  double upwind_flux(double u_cent, double u_neib, double a_dot_n) const;

  Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n) const;
};
