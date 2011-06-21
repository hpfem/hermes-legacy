#include "definitions.h"

using namespace WeakFormsHcurl;

scalar2 CustomInitialConditionWave::value (double x, double y) const 
{
  return scalar2(sin(x) * cos(y), -cos(x) * sin(y));
}

void CustomInitialConditionWave::derivatives (double x, double y, scalar2& dx, scalar2& dy) const 
{
    dx[0] = cos(x) * cos(y);
    dx[1] = sin(x) * sin(y);
    dy[0] = -sin(x) * sin(y);
    dy[1] = -cos(x) * cos(y);
}

Ord CustomInitialConditionWave::ord(Ord x, Ord y) const
{
  return Ord(10);
}

CustomWeakFormWave::CustomWeakFormWave(double c_squared) : WeakForm(2) 
{
  // Jacobian.
  add_matrix_form(new MatrixFormVolWave_0_1(c_squared));
  add_matrix_form(new MatrixFormVolWave_1_0);

  // Residual.
  add_vector_form(new VectorFormVolWave_0(c_squared));
  add_vector_form(new VectorFormVolWave_1());
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::MatrixFormVolWave_0_1::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  for (int i=0; i < n; i++) 
  {
    result += wt[i] * (u->dy[i] * v->val0[i] - u->dx[i] * v->val1[i]);
  }
  return c_squared * result;
}

scalar CustomWeakFormWave::MatrixFormVolWave_0_1::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                                                        Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakFormWave::MatrixFormVolWave_0_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

WeakForm::MatrixFormVol* CustomWeakFormWave::MatrixFormVolWave_0_1::clone() 
{
  return new MatrixFormVolWave_0_1(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::MatrixFormVolWave_1_0::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  for (int i=0; i < n; i++) 
  {
    result -= wt[i] * u->curl[i] * v->val[i];
  }
  return result;
}

scalar CustomWeakFormWave::MatrixFormVolWave_1_0::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                                                        Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakFormWave::MatrixFormVolWave_1_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

WeakForm::MatrixFormVol* CustomWeakFormWave::MatrixFormVolWave_1_0::clone() 
{
  return new MatrixFormVolWave_1_0(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::VectorFormVolWave_0::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                                            Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  Func<Scalar>* sln_prev = u_ext[1];

  for (int i = 0; i < n; i++) 
  {
    result += wt[i] * (sln_prev->dy[i] * v->val0[i] - sln_prev->dx[i] * v->val1[i]);
  }
  return c_squared * result;
}

scalar CustomWeakFormWave::VectorFormVolWave_0::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                                                      Geom<double> *e, ExtData<scalar> *ext) const 
{
  return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormWave::VectorFormVolWave_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                                                 ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

WeakForm::VectorFormVol* CustomWeakFormWave::VectorFormVolWave_0::clone() 
{
  return new VectorFormVolWave_0(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::VectorFormVolWave_1::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                                            Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  Func<Scalar>* sln_prev = u_ext[0];
      
  for (int i = 0; i < n; i++) 
  {
    result -= wt[i] * sln_prev->curl[i] * v->val[i];
  }
  return result;
}

scalar CustomWeakFormWave::VectorFormVolWave_1::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                                                      Geom<double> *e, ExtData<scalar> *ext) const 
{
  return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormWave::VectorFormVolWave_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

WeakForm::VectorFormVol* CustomWeakFormWave::VectorFormVolWave_1::clone() 
{
  return new VectorFormVolWave_1(*this);
}

