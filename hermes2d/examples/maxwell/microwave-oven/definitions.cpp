#include "definitions.h"

template<typename Real, typename Scalar>
Scalar CustomMatrixForm::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  cplx ikappa = cplx(0.0, kappa);
  Scalar result1 = 0;
  Scalar result2 = 0;
  Scalar result3 = 0;
  for (int i = 0; i < n; i++)
    result1 += wt[i] * gamma(e->elem_marker, e->x[i], e->y[i]) * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
  for (int i = 0; i < n; i++)
    result2 += wt[i] * er(e->elem_marker, e->x[i], e->y[i]) * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
  for (int i = 0; i < n; i++)
    result3 += wt[i] * (u->curl[i] * conj(v->curl[i]));
  return 1.0/mu_r * result3 - ikappa * sqrt(mu_0 / e_0) * result1 - sqr(kappa) * result2;
}

scalar CustomMatrixForm::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                               Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomMatrixForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                          Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

double CustomMatrixForm::gamma(int marker, double x, double y) const
{
  if (align_mesh && marker == 1) return 0.03;
  if (!align_mesh && in_load(x,y)) 
  {
    double cx = -0.152994121;  double cy =  0.030598824;
    double r = sqrt(sqr(cx - x) + sqr(cy - y));
    return (0.03 + 1)/2.0 - (0.03 - 1) * atan(10.0*(r -  0.043273273)) / M_PI;
  }
  return 0.0;
}

Ord CustomMatrixForm::gamma(int marker, Ord x, Ord y) const
{  
  return Ord(0.0); 
}

double CustomMatrixForm::er(int marker, double x, double y) const
{
  if (align_mesh && marker == 1) return 7.5;
  if (!align_mesh && in_load(x,y)) 
  {
    double cx = -0.152994121;  double cy =  0.030598824;
    double r = sqrt(sqr(cx - x) + sqr(cy - y));
    return (7.5 + 1)/2.0 - (7.5 - 1) * atan(10.0*(r -  0.043273273)) / M_PI;
  }
  return 1.0;
}

Ord CustomMatrixForm::er(int marker, Ord x, Ord y) const
{  
  return Ord(1.0); 
}

bool CustomMatrixForm::in_load(double x, double y) const
{
  double cx = -0.152994121;
  double cy =  0.030598824;
  double r = 0.043273273;
  if (sqr(cx - x) + sqr(cy - y) < sqr(r)) return true;
  else return false;
}


template<typename Real, typename Scalar>
Scalar CustomResidualForm::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                       Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  cplx ikappa = cplx(0.0, kappa);
  Scalar result1 = 0;
  Scalar result2 = 0;
  Scalar result3 = 0;
  for (int i = 0; i < n; i++)
    result1 += wt[i] * gamma(e->elem_marker, e->x[i], e->y[i]) * (u_ext[0]->val0[i] * conj(v->val0[i]) 
               + u_ext[0]->val1[i] * conj(v->val1[i]));
  for (int i = 0; i < n; i++)
    result2 += wt[i] * er(e->elem_marker, e->x[i], e->y[i]) * (u_ext[0]->val0[i] * conj(v->val0[i]) 
               + u_ext[0]->val1[i] * conj(v->val1[i]));
  for (int i = 0; i < n; i++)
    result3 += wt[i] * (u_ext[0]->curl[i] * conj(v->curl[i]));
  return 1.0/mu_r * result3 - ikappa * sqrt(mu_0 / e_0) * result1 - sqr(kappa) * result2;
}

scalar CustomResidualForm::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                                 Geom<double> *e, ExtData<scalar> *ext) const 
{
  return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
}

Ord CustomResidualForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                            Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

double CustomResidualForm::gamma(int marker, double x, double y) const
{
  if (align_mesh && marker == 1) return 0.03;
  if (!align_mesh && in_load(x,y)) 
  {
    double cx = -0.152994121;  double cy =  0.030598824;
    double r = sqrt(sqr(cx - x) + sqr(cy - y));
    return (0.03 + 1)/2.0 - (0.03 - 1) * atan(10.0*(r -  0.043273273)) / M_PI;
  }
  return 0.0;
}

Ord CustomResidualForm::gamma(int marker, Ord x, Ord y) const
{  
  return Ord(0.0); 
}

double CustomResidualForm::er(int marker, double x, double y) const
{
  if (align_mesh && marker == 1) return 7.5;
  if (!align_mesh && in_load(x,y)) 
  {
    double cx = -0.152994121;  double cy =  0.030598824;
    double r = sqrt(sqr(cx - x) + sqr(cy - y));
    return (7.5 + 1)/2.0 - (7.5 - 1) * atan(10.0*(r -  0.043273273)) / M_PI;
  }
  return 1.0;
}

Ord CustomResidualForm::er(int marker, Ord x, Ord y) const
{  
  return Ord(1.0); 
}

bool CustomResidualForm::in_load(double x, double y) const
{
  double cx = -0.152994121;
  double cy =  0.030598824;
  double r = 0.043273273;
  if (sqr(cx - x) + sqr(cy - y) < sqr(r)) return true;
  else return false;
}

template<typename Scalar, typename Real>
Scalar CustomVectorFormSurf::vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                              Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  cplx ii = cplx(0.0, 1.0);
  return ii * omega * J * int_v1<Real, Scalar>(n, wt, v); // just second component of v, since J = (0, J)
}

scalar CustomVectorFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], 
                                   Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return vector_form_surf<scalar, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomVectorFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                              Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
}


CustomWeakForm::CustomWeakForm(double e_0, double mu_0, double mu_r, double kappa, double omega, 
                               double J, bool align_mesh) : WeakForm(1) 
{
  // Jacobian forms - volumetric.
  add_matrix_form(new CustomMatrixForm(0, 0, e_0, mu_0, mu_r, kappa, omega, J, align_mesh));

  // Residual forms - volumetric.
  add_vector_form(new CustomResidualForm(0, e_0, mu_0, mu_r, kappa, omega, J, align_mesh));

  // Residual forms - surface.
  add_vector_form_surf(new CustomVectorFormSurf(omega, J));
}

