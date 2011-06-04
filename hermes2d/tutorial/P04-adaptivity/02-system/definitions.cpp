#include "definitions.h"

double CustomExactFunction1::val(double x) 
{
  return cos(M_PI*x/2);
}
  
double CustomExactFunction1::dx(double x) 
{
  return -sin(M_PI*x/2)*(M_PI/2.);
}
  
double CustomExactFunction1::ddxx(double x) 
{
  return -cos(M_PI*x/2)*(M_PI/2.)*(M_PI/2.);
}


double CustomExactFunction2::val(double x) 
{
  return 1. - (exp(K*x) + exp(-K*x))/(exp(K) + exp(-K));
}
  
double CustomExactFunction2::dx(double x) 
{
  return -K*(exp(K*x) - exp(-K*x))/(exp(K) + exp(-K));
}
  
double CustomExactFunction2::ddxx(double x) 
{
  return -K*K*(exp(K*x) + exp(-K*x))/(exp(K) + exp(-K));
}


CustomRightHandSide1::CustomRightHandSide1(double K, double d_u, double sigma)
      : HermesFunction(), d_u(d_u), sigma(sigma) 
{
  cef1 = new CustomExactFunction1();
  cef2 = new CustomExactFunction2(K);
}

scalar CustomRightHandSide1::value(double x, double y) const 
{
  double Laplace_u = cef1->ddxx(x) * cef1->val(y)
                     + cef1->val(x) * cef1->ddxx(y);
  double u = cef1->val(x) * cef1->val(y);
  double v = cef2->val(x) * cef2->val(y);
  return -d_u * d_u * Laplace_u - u + sigma * v;
}

Ord CustomRightHandSide1::ord(Ord x, Ord y) const 
{
  return Ord(10);
}

CustomRightHandSide1::~CustomRightHandSide1() 
{ 
  delete cef1; 
  delete cef2;
}

CustomRightHandSide2::CustomRightHandSide2(double K, double d_v)
      : HermesFunction(), d_v(d_v) 
{
  cef1 = new CustomExactFunction1();
  cef2 = new CustomExactFunction2(K);
}

scalar CustomRightHandSide2::value(double x, double y) const 
{
  double Laplace_v = cef2->ddxx(x) * cef2->val(y)
                     + cef2->val(x) * cef2->ddxx(y);
  double u = cef1->val(x) * cef1->val(y);
  double v = cef2->val(x) * cef2->val(y);
  return -d_v*d_v * Laplace_v - u + v;
}

Ord CustomRightHandSide2::ord(Ord x, Ord y) const 
{
  return Ord(10);
}

CustomRightHandSide2::~CustomRightHandSide2() 
{ 
  delete cef1; 
  delete cef2;
}


ExactSolutionFitzHughNagumo1::ExactSolutionFitzHughNagumo1(Mesh* mesh)
     : ExactSolutionScalar(mesh) 
{
  cef1 = new CustomExactFunction1();
}

scalar ExactSolutionFitzHughNagumo1::value (double x, double y) const 
{
  return cef1->val(x)*cef1->val(y);
}

void ExactSolutionFitzHughNagumo1::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  dx = cef1->dx(x)*cef1->val(y);
  dy = cef1->val(x)*cef1->ddxx(y);
}

Ord ExactSolutionFitzHughNagumo1::ord(Ord x, Ord y) const 
{
  return Ord(10);
}

ExactSolutionFitzHughNagumo1::~ExactSolutionFitzHughNagumo1() 
{
  delete cef1;
}

ExactSolutionFitzHughNagumo2::ExactSolutionFitzHughNagumo2(Mesh* mesh, double K)
     : ExactSolutionScalar(mesh) 
{
  cef2 = new CustomExactFunction2(K);
}

scalar ExactSolutionFitzHughNagumo2::value (double x, double y) const 
{
  return cef2->val(x)*cef2->val(y);
}

void ExactSolutionFitzHughNagumo2::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  dx = cef2->dx(x)*cef2->val(y);
  dy = cef2->val(x)*cef2->dx(y);
}

Ord ExactSolutionFitzHughNagumo2::ord(Ord x, Ord y) const 
{
  return Ord(10);
}

ExactSolutionFitzHughNagumo2::~ExactSolutionFitzHughNagumo2() 
{
  delete cef2;
}

scalar CustomResidual1::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                              Geom<double> *e, ExtData<scalar> *ext) const
{
   scalar result = 0;
   for (int i = 0; i < n; i++) 
   {
     result += wt[i] * (    sqr(d_u) * (u_ext[0]->dx[i]*v->dx[i] + u_ext[0]->dy[i]*v->dy[i]) 
                          - u_ext[0]->val[i]*v->val[i] 
                          + sigma*u_ext[1]->val[i]*v->val[i]
                          - g1->value(e->x[i], e->y[i])*v->val[i]
                       );
   }
 
   return result;
}

Ord CustomResidual1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                         Geom<Ord> *e, ExtData<Ord> *ext) const 
{
   Ord result = 0;
   for (int i = 0; i < n; i++) 
   {
     result += wt[i] * (    sqr(d_u) * (u_ext[0]->dx[i]*v->dx[i] + u_ext[0]->dy[i]*v->dy[i]) 
                          - u_ext[0]->val[i]*v->val[i] 
                          + sigma*u_ext[1]->val[i]*v->val[i]
                          - g1->ord(e->x[i], e->y[i])*v->val[i]
                        );
   }

   return result;
}

WeakForm::VectorFormVol* CustomResidual1::clone() 
{
  return new CustomResidual1(*this);
}

scalar CustomResidual2::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                              Geom<double> *e, ExtData<scalar> *ext) const
{
   scalar result = 0;
   for (int i = 0; i < n; i++) 
   {
     result += wt[i] * (    sqr(d_v) * (u_ext[1]->dx[i]*v->dx[i] + u_ext[1]->dy[i]*v->dy[i]) 
                          - u_ext[0]->val[i]*v->val[i] 
                          + u_ext[1]->val[i]*v->val[i]
                          - g2->value(e->x[i], e->y[i])*v->val[i]
                       );
   }
 
   return result;
  }

Ord CustomResidual2::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                         Geom<Ord> *e, ExtData<Ord> *ext) const 
{
   Ord result = 0;
   for (int i = 0; i < n; i++) 
   {
     result += wt[i] * (    sqr(d_v) * (u_ext[1]->dx[i]*v->dx[i] + u_ext[1]->dy[i]*v->dy[i]) 
                          - u_ext[0]->val[i]*v->val[i] 
                          + u_ext[1]->val[i]*v->val[i]
                          - g2->ord(e->x[i], e->y[i])*v->val[i]
                        );
   }

  return result;
} 

WeakForm::VectorFormVol* CustomResidual2::clone() 
{
  return new CustomResidual2(*this);
}

CustomWeakForm::CustomWeakForm(CustomRightHandSide1* g1, CustomRightHandSide2* g2) : WeakForm(2) 
{
  // Jacobian.
	add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, HERMES_ANY, new HermesFunction(g1->d_u * g1->d_u)));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 0, HERMES_ANY, new HermesFunction(-1.0)));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 1, HERMES_ANY, new HermesFunction(g1->sigma), HERMES_NONSYM));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(1, 0, HERMES_ANY, new HermesFunction(-1.0), HERMES_NONSYM));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(1, 1, HERMES_ANY, new HermesFunction(g2->d_v * g2->d_v)));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(1, 1, HERMES_ANY, new HermesFunction(1.0)));

  // Residual.
  add_vector_form(new CustomResidual1(g1->d_u, g1->sigma, g1));
  add_vector_form(new CustomResidual2(g2->d_v, g2));
}
