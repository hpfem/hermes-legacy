#include "definitions.h"

CustomWeakForm::CustomWeakForm() : WeakForm(1) 
{
  add_matrix_form(new MatrixFormVol(0, 0));
  add_matrix_form_surf(new MatrixFormSurface(0, 0));
  add_matrix_form_surf(new MatrixFormInterface(0, 0));

  add_vector_form(new ResidualForm(0));
  add_vector_form_surf(new ResidualFormSurface(0));
  add_vector_form_surf(new ResidualFormInterface(0));

  add_vector_form_surf(new VectorFormSurface(0));
};

CustomWeakForm::MatrixFormVol::MatrixFormVol(int i, int j) : WeakForm::MatrixFormVol(i, j)
{
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::MatrixFormVol::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * u->val[i] * (v->dy[i] - v->dx[i]);
  return result;
}

scalar CustomWeakForm::MatrixFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::MatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);

  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

CustomWeakForm::MatrixFormSurface::MatrixFormSurface(int i, int j) : WeakForm::MatrixFormSurf(i, j, HERMES_ANY) 
{
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::MatrixFormSurface::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    if(e->nx[i] - e->ny[i] > 0)
      result += wt[i] * u->val[i] * v->val[i] * (e->nx[i] - e->ny[i]);
  }

  return result;
}

scalar CustomWeakForm::MatrixFormSurface::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::MatrixFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);

  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

CustomWeakForm::MatrixFormInterface::MatrixFormInterface(int i, int j) : WeakForm::MatrixFormSurf(i, j, H2D_DG_INNER_EDGE) 
{
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::MatrixFormInterface::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;

  for (int i = 0; i < n; i++) {
    if(e->nx[i] - e->ny[i] > 0)
      result += wt[i] * u->get_val_central(i) * (v->get_val_central(i) - v->get_val_neighbor(i)) * (e->nx[i] - e->ny[i]);
    else
      result += wt[i] * u->get_val_neighbor(i) * (v->get_val_central(i) - v->get_val_neighbor(i)) * (e->nx[i] - e->ny[i]);
  }

  return result;
}

scalar CustomWeakForm::MatrixFormInterface::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::MatrixFormInterface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);
}

CustomWeakForm::ResidualForm::ResidualForm(int i) : WeakForm::VectorFormVol(i) 
{
}

scalar CustomWeakForm::ResidualForm::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * u_ext[0]->val[i] * (v->dy[i] - v->dx[i]);
  return result;
}

Ord CustomWeakForm::ResidualForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);

  Ord result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * u_ext[0]->val[i] * (v->dy[i] - v->dx[i]);
  return result;
}

template<typename Real>
Real CustomWeakForm::ResidualForm::F(Real x, Real y) const
{
  return 0;
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::ResidualForm::g(std::string ess_bdy_marker, Real x, Real y) const 
{
  if (ess_bdy_marker == left_bottom_bnd_part) return 1; else return 0;
}

CustomWeakForm::ResidualFormSurface::ResidualFormSurface(int i) : WeakForm::VectorFormSurf(i, HERMES_ANY) 
{
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::ResidualFormSurface::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    if(e->nx[i] - e->ny[i] > 0)
      result += wt[i] * u_ext[0]->val[i] * v->val[i] * (e->nx[i] - e->ny[i]);
  }

  return result;
}

scalar CustomWeakForm::ResidualFormSurface::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakForm::ResidualFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);
}

CustomWeakForm::ResidualFormInterface::ResidualFormInterface(int i) : WeakForm::VectorFormSurf(i, H2D_DG_INNER_EDGE)
{ 
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::ResidualFormInterface::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;

  for (int i = 0; i < n; i++) {
    if(e->nx[i] - e->ny[i] > 0)
      result += wt[i] * u_ext[0]->get_val_central(i) * v->val[i] * (e->nx[i] - e->ny[i]);
    else
      result += wt[i] * u_ext[0]->get_val_neighbor(i) * v->val[i] * (e->nx[i] - e->ny[i]);
  }

  return result;
}

scalar CustomWeakForm::ResidualFormInterface::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakForm::ResidualFormInterface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);
}

CustomWeakForm::VectorFormSurface::VectorFormSurface(int i) : WeakForm::VectorFormSurf(i, HERMES_ANY) 
{
}

scalar CustomWeakForm::VectorFormSurface::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  double result = 0;
  for (int i = 0; i < n; i++) {
    if(e->nx[i] - e->ny[i] <= 0)
      result += wt[i] * (e->x[i] + e->y[i]) * v->val[i] * (e->nx[i] - e->ny[i]);
  }

  return result;
}

Ord CustomWeakForm::VectorFormSurface::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(10);

  Ord result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (e->x[i] + e->y[i]) * 0.5 * v->val[i] * (e->nx[i] - e->ny[i]);
  return result;
}
