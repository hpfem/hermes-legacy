#include "definitions.h"

void CustomInitialCondition::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  scalar val = exp(-20*(x*x + y*y));
  dx = val * (-40.0 * x);
  dy = val * (-40.0 * y);
}

scalar CustomInitialCondition::value (double x, double y) const 
{
  return exp(-20*(x*x + y*y));
}

Ord CustomInitialCondition::ord(Ord x, Ord y) const 
{
  return exp(-20*(x*x + y*y));
}

CustomWeakFormGPRK::CustomWeakFormGPRK(double h, double m, double g, double omega) : WeakForm(1)
{
  // Jacobian volumetric part.
  add_matrix_form(new CustomFormMatrixFormVol(0, 0, h, m, g, omega));

  // Residual - volumetric.
  add_vector_form(new CustomFormVectorFormVol(0, h, m, g, omega));
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormGPRK::CustomFormMatrixFormVol::matrix_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                                   Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Scalar result = 0;
  Func<Scalar>* psi_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (h*h/(2*m*ii*h) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                     + 2*g/(ii*h)* u->val[i] * psi_prev_newton->val[i] * conj(psi_prev_newton->val[i]) * v->val[i]
                     + (g/ii*h) * psi_prev_newton->val[i] * psi_prev_newton->val[i] * u->val[i] * v->val[i]
                     + .5*m*omega*omega/(ii*h) * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * u->val[i] * v->val[i]);
  return result;
}

scalar CustomWeakFormGPRK::CustomFormMatrixFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                                                          Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return matrix_form_rk<double, scalar>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakFormGPRK::CustomFormMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                                                     Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form_rk<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

WeakForm::MatrixFormVol* CustomWeakFormGPRK::CustomFormMatrixFormVol::clone() 
{
  return new CustomFormMatrixFormVol(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormGPRK::CustomFormVectorFormVol::vector_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                                                   Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  scalar ii = cplx(0.0, 1.0);  // imaginary unit, ii^2 = -1

  Scalar result = 0;
  Func<Scalar>* psi_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (h*h/(2*m*ii*h) * (psi_prev_newton->dx[i] * v->dx[i] + psi_prev_newton->dy[i] * v->dy[i])
                     + 2*g/(ii*h)* psi_prev_newton->val[i] *  psi_prev_newton->val[i] * conj(psi_prev_newton->val[i]) * v->val[i]
                     + (g/ii*h) * psi_prev_newton->val[i] * psi_prev_newton->val[i] * psi_prev_newton->val[i] * v->val[i]
                     + .5*m*omega*omega/(ii*h) * (e->x[i] * e->x[i] + e->y[i] * e->y[i]) * psi_prev_newton->val[i] * v->val[i]);

  return result;
}

scalar CustomWeakFormGPRK::CustomFormVectorFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e,
                                                           ExtData<scalar> *ext) const 
{
  return vector_form_rk<double, scalar>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormGPRK::CustomFormVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form_rk<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

WeakForm::VectorFormVol* CustomWeakFormGPRK::CustomFormVectorFormVol::clone() 
{
  return new CustomFormVectorFormVol(*this);
}

