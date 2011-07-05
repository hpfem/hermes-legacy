#include "definitions.h"

double CustomExactSolution::value (double x, double y) const 
{
  return 0;
}

void CustomExactSolution::derivatives(double x, double y, scalar& dx, scalar& dy) const 
{
  dx = 0; 
  dy = 0;
}

Ord CustomExactSolution::ord(Ord x, Ord y) const 
{
  return Ord(20);
}

double CustomFunction::value(double x, double y, double t) const 
{
  double S = s;
  double C = c;
  double f = -2*pow(S, 3)*pow(x, 2)*(-t + sqrt(pow(x, 2) + pow(y, 2)))*(x - x0)*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*(pow(x, 2) + pow(y, 2))*pow(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1, 2)) - 2*pow(S, 3)*pow(y, 2)*(-t + sqrt(pow(x, 2) + pow(y, 2)))*(x - x0)*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*(pow(x, 2) + pow(y, 2))*pow(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1, 2)) - S*pow(x, 2)*(x - x0)*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*pow(pow(x, 2) + pow(y, 2), 3.0/2.0)*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) + 2*S*x*(x - x0)*(y - y0)*(y - y1)*atan(t)/(C*sqrt(pow(x, 2) + pow(y, 2))*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) + 2*S*x*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*sqrt(pow(x, 2) + pow(y, 2))*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) - S*pow(y, 2)*(x - x0)*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*pow(pow(x, 2) + pow(y, 2), 3.0/2.0)*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) + 2*S*y*(x - x0)*(x - x1)*(y - y0)*atan(t)/(C*sqrt(pow(x, 2) + pow(y, 2))*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) + 2*S*y*(x - x0)*(x - x1)*(y - y1)*atan(t)/(C*sqrt(pow(x, 2) + pow(y, 2))*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) + S*(x - x0)*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) + 2*S*(x - x0)*(x - x1)*(y - y0)*(y - y1)*atan(t)/(C*sqrt(pow(x, 2) + pow(y, 2))*(pow(S, 2)*pow(-t + sqrt(pow(x, 2) + pow(y, 2)), 2) + 1)) - 2*(x - x0)*(x - x1)*(-atan(S*(-t + sqrt(pow(x, 2) + pow(y, 2)))) + M_PI/2)*atan(t)/C - 2*(y - y0)*(y - y1)*(-atan(S*(-t + sqrt(pow(x, 2) + pow(y, 2)))) + M_PI/2)*atan(t)/C + (x - x0)*(x - x1)*(y - y0)*(y - y1)*(-atan(S*(-t + sqrt(pow(x, 2) + pow(y, 2)))) + M_PI/2)/(C*(pow(t, 2) + 1));

  return f;
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord(20);
}

CustomVectorFormVol::CustomVectorFormVol(int i, std::string area,
  HermesFunction* coeff,
  GeomType gt)
  : WeakForm::VectorFormVol(i, area), coeff(coeff), gt(gt)
{
  // If coeff is HERMES_ONE, initialize it to be constant 1.0.
  if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
}

CustomVectorFormVol::CustomVectorFormVol(int i, Hermes::vector<std::string> areas,
  HermesFunction* coeff,
  GeomType gt)
  : WeakForm::VectorFormVol(i, areas), coeff(coeff), gt(gt)
{
  // If coeff is HERMES_ONE, initialize it to be constant 1.0.
  if (coeff == HERMES_ONE) this->coeff = new HermesFunction(1.0);
}

CustomVectorFormVol::~CustomVectorFormVol() 
{
  // FIXME: Should be deleted here only if it was created here.
  //if (coeff != HERMES_ONE) delete coeff;
};

scalar CustomVectorFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
  Geom<double> *e, ExtData<scalar> *ext) const 
{
  scalar result = 0;
  if (gt == HERMES_PLANAR) {
    for (int i = 0; i < n; i++) {
      result += wt[i] * static_cast<CustomFunction*>(coeff)->value(e->x[i], e->y[i], this->get_current_stage_time()) * v->val[i];
    }
  }
  else {
    if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * e->y[i] * static_cast<CustomFunction*>(coeff)->value(e->x[i], e->y[i], this->get_current_stage_time()) * v->val[i];
      }
    }
    else {
      for (int i = 0; i < n; i++) {
        result += wt[i] * e->x[i] * static_cast<CustomFunction*>(coeff)->value(e->x[i], e->y[i], this->get_current_stage_time()) * v->val[i];
      }
    }
  }
  return result;
}

Ord CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
  Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  Ord result = 0;
  if (gt == HERMES_PLANAR) {
    for (int i = 0; i < n; i++) {
      result += wt[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
    }
  }
  else {
    if (gt == HERMES_AXISYM_X) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
      }
    }
    else {
      for (int i = 0; i < n; i++) {
        result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
      }
    }
  }

  return result;
}

WeakForm::VectorFormVol* CustomVectorFormVol::clone() 
{
  return new CustomVectorFormVol(*this);
}

CustomWeakFormPoisson::CustomWeakFormPoisson(std::string area,
  HermesFunction* coeff,
  HermesFunction* f,
  GeomType gt) : DefaultWeakFormPoisson()
{
  // Jacobian.
  // NOTE: The flag HERMES_NONSYM is important here.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, area, coeff, HERMES_NONSYM, gt));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, area, coeff, gt));
  add_vector_form(new CustomVectorFormVol(0, area, f, gt));
};