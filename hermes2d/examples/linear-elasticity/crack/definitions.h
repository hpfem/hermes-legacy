#include "hermes2d.h"

//#define USE_MULTICOMPONENT_FORMS

class CustomWeakFormLinearElasticity : public WeakForm
{
public:
  CustomWeakFormLinearElasticity(double E, double nu, double rho_g);
};
