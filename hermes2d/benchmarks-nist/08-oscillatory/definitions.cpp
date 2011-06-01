#include "definitions.h"

double CustomFunction::value(double x, double y) const 
{
  return (-sin(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),4)
           + 2*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),2)*pow((pow(x,2) + pow(y,2)),(1.0/2.0)))
           + pow(x,2)*sin(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),4)*(pow(x,2) + pow(y,2)))
           + pow(y,2)*sin(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),4)*(pow(x,2) + pow(y,2)))
           - pow(x,2)*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),2)*pow((pow(x,2) + pow(y,2)),(3.0/2.0)))
           - pow(y,2)*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),2)*pow((pow(x,2) + pow(y,2)),(3.0/2.0)))
           - 2*pow(x,2)*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),3)*(pow(x,2) + pow(y,2)))
           - 2*pow(y,2)*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha
           + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),3)*(pow(x,2) + pow(y,2))));
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord(10);
}


double CustomExactSolution::value(double x, double y) const 
{
  double r = sqrt(x*x + y*y);
  return sin(1/(alpha + r));
}

void CustomExactSolution::derivatives(double x, double y, scalar& dx, scalar& dy) const 
{
  double r = sqrt(x*x + y*y);
  double h = 1/(alpha + r);
  dx = -cos(h) * h * h * x / r;
  dy = -cos(h) * h * h * y / r;
}

Ord CustomExactSolution::ord(Ord x, Ord y) const 
{
  return Ord(10);
}


CustomWeakForm::CustomWeakForm(CustomFunction* f) : WeakForm(1) 
{
  // Jacobian.
  add_matrix_form(new CustomMatrixFormVol(0, 0, f->alpha));

  // Residual.
  add_vector_form(new CustomVectorFormVol(0, f));
}


template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormVol::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar val = 0;
  for (int i=0; i < n; i++) 
  {
    Scalar x = e->x[i];
    Scalar y = e->y[i];
    Scalar r = sqrt(x*x + y*y);
    Scalar h = 1/(alpha + r);
    Scalar grad_u_grad_v = u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i];
    val += wt[i] * (grad_u_grad_v - pow(h, 4) * u->val[i] * v->val[i]);
  }

  return val;
}

scalar CustomWeakForm::CustomMatrixFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                                  Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::CustomMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                             Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}


scalar CustomWeakForm::CustomVectorFormVol::value(int n, double *wt, Func<scalar> *u_ext[],
                                                  Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  scalar val = 0;
  for (int i=0; i < n; i++) 
  {
    scalar x = e->x[i];
    scalar y = e->y[i];
    scalar r = sqrt(x*x + y*y);
    scalar h = 1/(f->alpha + r);
    scalar grad_u_grad_v = u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i];
    val += wt[i] * (grad_u_grad_v - pow(h, 4) * u_ext[0]->val[i] * v->val[i]);
    val -= wt[i] * f->value(e->x[i], e->y[i]) * v->val[i]; 
  }

  return val;
}

Ord CustomWeakForm::CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                             Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  Ord val = 0;
  for (int i=0; i < n; i++) 
  {
    Ord x = e->x[i];
    Ord y = e->y[i];
    Ord r = sqrt(x*x + y*y);
    Ord h = 1/(f->alpha + r);
    Ord grad_u_grad_v = u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i];
    val += wt[i] * (grad_u_grad_v - pow(h, 4) * u_ext[0]->val[i] * v->val[i]);
    val -= wt[i] * f->value(e->x[i], e->y[i]) * v->val[i]; 
  }

  return val;
}
