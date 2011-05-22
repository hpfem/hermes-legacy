#include "hermes2d.h"

/* Weak form */

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(double K_squared);
};


